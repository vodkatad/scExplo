
from __future__ import annotations

import os
import re
import itertools
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from sklearn.metrics import silhouette_score

# Optional: external integrations
import bbknn
import scanpy.external as sce

# -------------------- Constants --------------------

DEFAULT_BATCH_KEY = "sample"
DEFAULT_UID_SEP = "|"

# Cell-cycle gene sets (human; Seurat/Tirosh)
S_GENES = [
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7",
    "DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76",
    "SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6",
    "EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1",
    "E2F8"
]
G2M_GENES = [
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2",
    "CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","CKAP2L","CKAP2","AURKB",
    "BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","CDC20","TTK",
    "CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23",
    "HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3",
    "GAS2L3","CBX5","CENPA"
]

# -------------------- File / naming helpers --------------------

def safe_dirname(s: str) -> str:
    """Make a filesystem-safe name out of arbitrary text."""
    return re.sub(r"[^A-Za-z0-9_=.,+-]", "_", str(s))

def outdir_for(root: str, method: str, params: Dict) -> str:
    """results_grid/<method>/<k=v__k=v__...>/"""
    parts = [f"{k}={params[k]}" for k in sorted(params)]
    name = "__".join(safe_dirname(p) for p in parts)
    path = os.path.join(root, method, name)
    os.makedirs(path, exist_ok=True)
    return path

# -------------------- Data preparation --------------------

def strip_prefix_from_genes(
    df: pd.DataFrame,
    meta_cols: Tuple[str, ...] = ("cell_id","sample"),
    sep: str = ":",
    on_duplicate: str = "first",  # "error" | "first" | "mean" | "sum" | "suffix"
) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """
    Rename gene columns by removing the prefix before `sep`.
    Handle duplicates according to `on_duplicate` policy.
    Returns: (df_out, dup_report) where dup_report maps gene -> original columns
    """
    df = df.copy()
    meta_cols = list(meta_cols)
    gene_cols = [c for c in df.columns if c not in meta_cols]

    # map old gene col -> symbol after sep
    new_names = {c: c.split(sep)[-1] for c in gene_cols}
    df_ren = df.rename(columns=new_names)

    # build duplicate report
    dup_report: Dict[str, List[str]] = {}
    counts = pd.Series([new_names[c] for c in gene_cols]).value_counts()
    dups = counts[counts > 1].index.tolist()
    if dups:
        for g in dups:
            dup_report[g] = [c for c in gene_cols if new_names[c] == g]

        if on_duplicate == "error":
            raise ValueError(f"Duplicate genes after renaming: {dup_report}")

        elif on_duplicate in ("mean","sum"):
            agg_func = np.nanmean if on_duplicate == "mean" else np.nansum
            gdf = df_ren.drop(columns=meta_cols, errors="ignore")
            # collapse duplicates by aggregating columns with same name
            collapsed = {}
            for g, _cols in dup_report.items():
                mask = gdf.columns == g
                if mask.sum() > 1:
                    collapsed[g] = agg_func(gdf.loc[:, mask].to_numpy(), axis=1)
                else:
                    collapsed[g] = gdf.loc[:, g].to_numpy()
            keep_mask = ~gdf.columns.duplicated(keep=False)
            gdf_unique = gdf.loc[:, keep_mask].copy()
            for g, vec in collapsed.items():
                gdf_unique[g] = vec
            df_out = pd.concat([df_ren.loc[:, meta_cols], gdf_unique], axis=1)

        elif on_duplicate == "first":
            cols_meta = list(meta_cols)
            cols_gene_unique = ~df_ren.columns.duplicated(keep="first")
            df_out = pd.concat([df_ren.loc[:, cols_meta],
                                df_ren.loc[:, cols_gene_unique & ~df_ren.columns.isin(cols_meta)]],
                               axis=1)
        elif on_duplicate == "suffix":
            seen = {}
            newcols = []
            for c in df_ren.columns:
                if c in meta_cols:
                    newcols.append(c); continue
                name = c
                if name not in seen:
                    seen[name] = 1
                    newcols.append(name)
                else:
                    newcols.append(f"{name}.{seen[name]}")
                    seen[name] += 1
            df_out = df_ren.copy()
            df_out.columns = newcols
        else:
            raise ValueError(f"Unsupported on_duplicate={on_duplicate}")
    else:
        df_out = df_ren

    # reorder: meta first
    gene_cols_new = [c for c in df_out.columns if c not in meta_cols]
    df_out = pd.concat([df_out.loc[:, meta_cols], df_out.loc[:, gene_cols_new]], axis=1)

    return df_out, dup_report

def make_anndata_from_df(
    df: pd.DataFrame,
    cell_col: str = "cell_id",
    sample_col: str = "sample",
    uid_sep: str = DEFAULT_UID_SEP,
    set_raw: bool = True,
) -> ad.AnnData:
    """Create AnnData from a wide table with row=cell, gene columns + meta (cell_id, sample)."""
    meta_cols = [cell_col, sample_col]
    gene_cols = [c for c in df.columns if c not in meta_cols]
    uid = df[sample_col].astype(str) + uid_sep + df[cell_col].astype(str)

    adata = ad.AnnData(
        X = df[gene_cols].to_numpy(),
        obs = pd.DataFrame({
            "cell_id": df[cell_col].astype(str).values,
            "sample": df[sample_col].astype("category").values,
            "uid": uid.values
        }, index=uid),
        var = pd.DataFrame(index=pd.Index(gene_cols, name="gene"))
    )
    adata.obs_names = adata.obs["uid"].astype(str)
    if set_raw:
        adata.raw = adata
    return adata

# -------------------- Cell cycle --------------------

def _match_genes_names(var_names, genes):
    vset = set(var_names)
    present = [g for g in genes if g in vset]
    if present:
        return present
    upmap = {g.upper(): g for g in var_names}
    return [upmap[g.upper()] for g in genes if g.upper() in upmap]

def score_cell_cycle(adata_full: ad.AnnData,
                     s_genes: Optional[List[str]] = None,
                     g2m_genes: Optional[List[str]] = None) -> None:
    """
    Score cell-cycle on the full object (before HVG). Leaves results in adata_full.obs:
    'S_score', 'G2M_score', 'phase'.
    """
    s_list = s_genes if s_genes is not None else S_GENES
    g_list = g2m_genes if g2m_genes is not None else G2M_GENES
    s_present = _match_genes_names(adata_full.var_names, s_list)
    g_present = _match_genes_names(adata_full.var_names, g_list)
    if len(s_present) < 5 or len(g_present) < 5:
        print(f"⚠️ Cell-cycle: small overlap (S={len(s_present)}, G2M={len(g_present)}).")
    sc.tl.score_genes_cell_cycle(adata_full, s_genes=s_present, g2m_genes=g_present)

# -------------------- HVG / PCA --------------------

def select_hvg_cell_ranger(adata: ad.AnnData, n_top_genes: int = 1000,
                           batch_key: str = DEFAULT_BATCH_KEY, subset: bool = True) -> None:
    """Select HVG with flavor=cell_ranger (batch-aware if batch_key provided)."""
    sc.pp.highly_variable_genes(adata, flavor="cell_ranger",
                                n_top_genes=n_top_genes, batch_key=batch_key, subset=subset)

def scale_and_pca(adata: ad.AnnData, n_comps: int = 120,
                  max_value: Optional[float] = 10, random_state: int = 42) -> None:
    """Scale and run PCA (store in adata.uns['pca'])."""
    sc.pp.scale(adata, max_value=max_value)
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=random_state)

def pcs_needed_for_target(adata: ad.AnnData, target: float = 0.90) -> Optional[int]:
    """Return #PCs to reach target cumulative variance if possible (uses last computed PCA)."""
    evr = adata.uns.get("pca", {}).get("variance_ratio", None)
    if evr is None:
        return None
    cum = np.cumsum(evr)
    return int(np.searchsorted(cum, target) + 1) if cum[-1] >= target else None

# -------------------- Integrations --------------------

def integrate_bbknn(adata_base: ad.AnnData, n_pcs: int, umap_n_neighbors: int = 15,
                    neighbors_within_batch: int = 5, trim: Optional[int] = None,
                    batch_key: str = DEFAULT_BATCH_KEY, random_state: int = 42) -> ad.AnnData:
    """Return a copy integrated with BBKNN; keeps .raw; computes neighbors+UMAP."""
    adx = adata_base.copy()
    sc.pp.neighbors(adx, use_rep="X_pca", n_neighbors=umap_n_neighbors, n_pcs=n_pcs, random_state=random_state)
    bbknn.bbknn(adx, batch_key=batch_key, neighbors_within_batch=neighbors_within_batch, trim=trim, n_pcs=n_pcs)
    sc.tl.umap(adx, min_dist=0.3, random_state=random_state)
    return adx

def integrate_harmony(base, n_pcs, umap_n_neighbors, theta, lam,
                      batch_key="sample", random_state=42):
    import scanpy as sc
    import scanpy.external as sce

    adx = base.copy()
    # usa SEMPRE 'lambda' (dict expansion, perché 'lambda' è parola riservata in Python)
    sce.pp.harmony_integrate(
        adx, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony",
        theta=theta, **{"lambda": lam}, max_iter_harmony=20
    )

    sc.pp.neighbors(
        adx, use_rep="X_pca_harmony", n_neighbors=umap_n_neighbors, n_pcs=n_pcs,
        random_state=random_state
    )
    sc.tl.umap(adx, min_dist=0.3, random_state=random_state)  # niente n_neighbors qui
    return adx


# -------------------- Metrics --------------------

def batch_entropy_from_graph(adata: ad.AnnData, batch_key: str = DEFAULT_BATCH_KEY) -> float:
    """Normalized entropy of batch labels in the neighbor graph (↑ better)."""
    B = adata.obs[batch_key].astype(str).values
    batches = np.unique(B)
    if len(batches) < 2:
        return np.nan
    b2i = {b:i for i,b in enumerate(batches)}
    G = adata.obsp.get("connectivities")
    if G is None:
        return np.nan
    if not sparse.isspmatrix_csr(G):
        G = sparse.csr_matrix(G)
    log_nb = np.log(len(batches))
    ent = np.zeros(adata.n_obs, float)
    for i in range(adata.n_obs):
        row = G.getrow(i); idx, w = row.indices, row.data
        if w.size == 0: continue
        sums = np.bincount([b2i[B[j]] for j in idx], weights=w, minlength=len(batches))
        p = sums / (sums.sum() + 1e-12)
        h = -(p[p>0]*np.log(p[p>0])).sum()
        ent[i] = h / (log_nb if log_nb>0 else 1.0)
    return float(np.nanmean(ent))

def safe_silhouette(X: np.ndarray, labels: np.ndarray) -> float:
    """Return silhouette score or NaN if not computable."""
    labels = np.asarray(labels)
    n_labels = len(np.unique(labels))
    if X is None or X.shape[0] < 3 or n_labels < 2 or n_labels >= X.shape[0]:
        return np.nan
    try:
        return float(silhouette_score(X, labels))
    except Exception:
        return np.nan

# -------------------- Plots / Save --------------------

def save_umap_gene_panel(adx: ad.AnnData, outdir: str, genes: List[str]) -> None:
    """Save a grid of UMAPs for given genes; falls back to raw if genes not in .var."""
    # if genes not in var, try raw (case-sensitive first, then case-insensitive)
    use_raw = False
    present = [g for g in genes if g in adx.var_names]
    if not present and adx.raw is not None:
        rv = list(adx.raw.var_names)
        present = [g for g in genes if g in rv]
        if not present:
            up = {x.upper(): x for x in rv}
            present = [up[g.upper()] for g in genes if g.upper() in up]
        if present:
            use_raw = True
    if not present:
        print("⚠️ save_umap_gene_panel: none of the requested genes are present (var/raw).")
        return
    sc.settings.figdir = outdir
    sc.pl.umap(adx, color=present, ncols=4, wspace=0.3, frameon=False,
               use_raw=use_raw, show=False, save="_genes_panel.png")

def save_umap_cellcycle(adx: ad.AnnData, outdir: str) -> None:
    """Save UMAP colored by 'phase' if present."""
    if "phase" not in adx.obs.columns:
        print("⚠️ save_umap_cellcycle: 'phase' not present in obs.")
        return
    sc.settings.figdir = outdir
    sc.pl.umap(adx, color=["phase"], frameon=False, show=False, save="_cellcycle.png")

def save_umap_clusters(adx: ad.AnnData, outdir: str, cluster_key: str) -> None:
    """Save UMAP colored by the given cluster labels key (e.g., 'leiden_r0.8')."""
    sc.settings.figdir = outdir
    sc.pl.umap(adx, color=[cluster_key], frameon=False, show=False, save=f"_{cluster_key}.png")

def save_markers(adx: ad.AnnData, outdir: str, cluster_key: str, n_top: int = 25) -> None:
    """Run rank_genes_groups and save both PNG and CSV of top markers per cluster."""
    use_raw = adx.raw is not None
    sc.tl.rank_genes_groups(adx, groupby=cluster_key, method="wilcoxon", use_raw=use_raw)
    sc.settings.figdir = outdir
    sc.pl.rank_genes_groups(adx, n_genes=n_top, sharey=False, show=False,
                            save=f"_{cluster_key}_top{n_top}.png")
    df = sc.get.rank_genes_groups_df(adx, group=None)
    top = df.groupby("group", sort=False).head(n_top)
    top.to_csv(os.path.join(outdir, f"markers_{cluster_key}_top{n_top}.csv"), index=False)

# -------------------- Cluster selection --------------------

def pick_best_cluster(rows_local: List[Dict]) -> Tuple[Optional[str], Optional[Dict]]:
    """Pick the 'best' clustering among rows_local by metrics (see docstring)."""
    if not rows_local:
        return None, None
    rl = [r for r in rows_local if not np.isnan(r.get("silhouette_cluster", np.nan))]
    if not rl:
        rl = rows_local[:]
    def algo_priority(key: str) -> int:
        return 1 if key.startswith("leiden_") else 0
    rl.sort(key=lambda r: (
        r.get("silhouette_cluster", -np.inf),
        r.get("batch_entropy",    -np.inf),
        -r.get("silhouette_batch", np.inf),
        algo_priority(r["clustering"])
    ), reverse=True)
    best = rl[0]
    return best["clustering"], best

def save_summary_umap(adx: ad.AnnData, outdir: str, best_key: Optional[str]) -> None:
    """Save a single PNG with UMAP colored by sample, phase, and best cluster key."""
    if best_key is None:
        return
    sc.settings.figdir = outdir
    sc.pl.umap(adx, color=["sample", "phase", best_key], ncols=3,
               frameon=False, wspace=0.3, show=False, save="_summary.png")
    with open(os.path.join(outdir, "best_cluster.txt"), "w") as f:
        f.write(best_key + "\n")

# -------------------- Grid search (no KMeans) --------------------

def run_grid(
    adata_base: ad.AnnData,
    grid: Dict,
    out_root: str = "./results_grid",
    batch_key: str = DEFAULT_BATCH_KEY,
    gene_panel: Optional[List[str]] = None,
    random_state: int = 42
) -> pd.DataFrame:
    """
    Run integrations (BBKNN/Harmony) across hyperparameters, do Leiden/Louvain,
    compute metrics, and save UMAPs (gene panel, cell cycle, clusters) + markers.
    grid example:
        {
          "n_pcs": [80, 100],
          "umap_n_neighbors": 15,
          "bbknn_neighbors_within_batch": [3,5],
          "bbknn_trim": [None, 100],
          "harmony_theta": [2.0, 4.0],
          "harmony_lambda": [2.0],
          "leiden_res": [0.8, 1.2],
          "louvain_res": [0.8, 1.2],
        }
    """
    os.makedirs(out_root, exist_ok=True)
    rows = []
    NN = grid["umap_n_neighbors"]
    panel = gene_panel or []

    # ------------ BBKNN -------------
    for n_pcs in grid["n_pcs"]:
        for nwb in grid["bbknn_neighbors_within_batch"]:
            for trim in grid["bbknn_trim"]:
                adx = integrate_bbknn(adata_base, n_pcs=n_pcs, umap_n_neighbors=NN,
                                      neighbors_within_batch=nwb, trim=trim,
                                      batch_key=batch_key, random_state=random_state)
                params = {"npcs": n_pcs, "umap_nn": NN, "nwb": nwb, "trim": trim}
                outdir = outdir_for(out_root, "bbknn", params)

                save_umap_gene_panel(adx, outdir, panel)
                save_umap_cellcycle(adx, outdir)
                save_umap_sample(adx, outdir, sample_key=batch_key)

                # clustering + metrics
                rows_local = []
                for res in grid["leiden_res"]:
                    key = f"leiden_r{res}"
                    sc.tl.leiden(adx, resolution=res, key_added=key, random_state=random_state)
                    labels = adx.obs[key].astype(str).values
                    row = {
                        "method": "bbknn", **params, "clustering": key,
                        "n_clusters": int(len(np.unique(labels))),
                        "silhouette_cluster": safe_silhouette(adx.obsm["X_umap"], labels),
                        "silhouette_batch":  safe_silhouette(adx.obsm["X_umap"], adx.obs[batch_key].astype(str).values),
                        "batch_entropy":     batch_entropy_from_graph(adx, batch_key=batch_key),
                    }
                    rows.append(row); rows_local.append(row)
                    # per-clustering saves
                    subdir = os.path.join(outdir, key); os.makedirs(subdir, exist_ok=True)
                    pd.DataFrame([row]).to_csv(os.path.join(subdir, "metrics.csv"), index=False)
                    save_umap_clusters(adx, outdir, key)
                    save_markers(adx, outdir, key, n_top=25)

                for res in grid["louvain_res"]:
                    key = f"louvain_r{res}"
                    sc.tl.louvain(adx, resolution=res, key_added=key, random_state=random_state)
                    labels = adx.obs[key].astype(str).values
                    row = {
                        "method": "bbknn", **params, "clustering": key,
                        "n_clusters": int(len(np.unique(labels))),
                        "silhouette_cluster": safe_silhouette(adx.obsm["X_umap"], labels),
                        "silhouette_batch":  safe_silhouette(adx.obsm["X_umap"], adx.obs[batch_key].astype(str).values),
                        "batch_entropy":     batch_entropy_from_graph(adx, batch_key=batch_key),
                    }
                    rows.append(row); rows_local.append(row)
                    subdir = os.path.join(outdir, key); os.makedirs(subdir, exist_ok=True)
                    pd.DataFrame([row]).to_csv(os.path.join(subdir, "metrics.csv"), index=False)
                    save_umap_clusters(adx, outdir, key)
                    save_markers(adx, outdir, key, n_top=25)

                best_key, best_row = pick_best_cluster(rows_local)
                save_summary_umap(adx, outdir, best_key)
                
                # save integrated object
                adx.write(os.path.join(outdir, "bbknn_integrated.h5ad"), compression="lzf")

    # ------------ HARMONY -------------
    for n_pcs in grid["n_pcs"]:
        for theta in grid["harmony_theta"]:
            for lam in grid["harmony_lambda"]:
                adx = integrate_harmony(adata_base, n_pcs=n_pcs, umap_n_neighbors=NN,
                        theta=theta, lam=lam, batch_key=batch_key, random_state=random_state)
                params = {"npcs": n_pcs, "umap_nn": NN, "theta": theta, "lambda": lam}
                outdir = outdir_for(out_root, "harmony", params)

                save_umap_gene_panel(adx, outdir, panel)
                save_umap_cellcycle(adx, outdir)
                save_umap_sample(adx, outdir, sample_key=batch_key)
                rows_local = []
                for res in grid["leiden_res"]:
                    key = f"leiden_r{res}"
                    sc.tl.leiden(adx, resolution=res, key_added=key, random_state=random_state)
                    labels = adx.obs[key].astype(str).values
                    row = {
                        "method": "harmony", **params, "clustering": key,
                        "n_clusters": int(len(np.unique(labels))),
                        "silhouette_cluster": safe_silhouette(adx.obsm["X_umap"], labels),
                        "silhouette_batch":  safe_silhouette(adx.obsm["X_umap"], adx.obs[batch_key].astype(str).values),
                        "batch_entropy":     batch_entropy_from_graph(adx, batch_key=batch_key),
                    }
                    rows.append(row); rows_local.append(row)
                    subdir = os.path.join(outdir, key); os.makedirs(subdir, exist_ok=True)
                    pd.DataFrame([row]).to_csv(os.path.join(subdir, "metrics.csv"), index=False)
                    save_umap_clusters(adx, outdir, key)
                    save_markers(adx, outdir, key, n_top=25)

                for res in grid["louvain_res"]:
                    key = f"louvain_r{res}"
                    sc.tl.louvain(adx, resolution=res, key_added=key, random_state=random_state)
                    labels = adx.obs[key].astype(str).values
                    row = {
                        "method": "harmony", **params, "clustering": key,
                        "n_clusters": int(len(np.unique(labels))),
                        "silhouette_cluster": safe_silhouette(adx.obsm["X_umap"], labels),
                        "silhouette_batch":  safe_silhouette(adx.obsm["X_umap"], adx.obs[batch_key].astype(str).values),
                        "batch_entropy":     batch_entropy_from_graph(adx, batch_key=batch_key),
                    }
                    rows.append(row); rows_local.append(row)
                    subdir = os.path.join(outdir, key); os.makedirs(subdir, exist_ok=True)
                    pd.DataFrame([row]).to_csv(os.path.join(subdir, "metrics.csv"), index=False)
                    save_umap_clusters(adx, outdir, key)
                    save_markers(adx, outdir, key, n_top=25)

                best_key, best_row = pick_best_cluster(rows_local)
                save_summary_umap(adx, outdir, best_key)
                adx.write(os.path.join(outdir, "harmony_integrated.h5ad"), compression="lzf")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out_root, "summary_metrics.csv"), index=False)
    return df

# -------------------- Extra: attach external UMAP coords by uid --------------------

def add_umap_from_df_uid(adata: ad.AnnData, df_umap: pd.DataFrame,
                         x_col: str = "X", y_col: str = "Y",
                         uid_col: str = "uid",
                         id_col: str = "cell_id", sample_col: str = "sample",
                         obsm_key: str = "X_umap") -> None:
    """
    Attach external UMAP coordinates from a DataFrame that has either 'uid' or (sample, cell_id).
    """
    df = df_umap.copy()
    if uid_col not in df.columns:
        need = {id_col, sample_col}
        if not need.issubset(df.columns):
            raise ValueError(f"Need '{uid_col}' or both {need} in df_umap.")
        df[uid_col] = df[sample_col].astype(str) + DEFAULT_UID_SEP + df[id_col].astype(str)
    df = df[[uid_col, x_col, y_col]].drop_duplicates(uid_col).set_index(uid_col)
    missing = set(adata.obs_names) - set(df.index)
    if missing:
        ex = list(missing)[:5]
        raise ValueError(f"Missing UMAP coords for {len(missing)} cells; examples: {ex}")
    coords = df.reindex(adata.obs_names)[[x_col, y_col]].to_numpy()
    if np.isnan(coords).any():
        raise ValueError("NaNs in aligned UMAP coords.")
    adata.obsm[obsm_key] = coords
def save_umap_sample(adx, outdir, sample_key="sample"):
    """Salva una UMAP colorata per sample."""
    if sample_key not in adx.obs.columns:
        print(f"⚠️ '{sample_key}' non presente in obs; salto il plot.")
        return
    sc.settings.figdir = outdir
    sc.pl.umap(
        adx,
        color=[sample_key],
        frameon=False,
        legend_loc="right", 
        show=False,
        save="_sample.png",    # -> file: umap_sample.png
    )


def make_base_cc_regressed(adata_full, n_hvg=1000, batch_key="sample", n_pcs=120, rng=42):
    ad = adata_full.copy()

    sc.pp.highly_variable_genes(ad, flavor="cell_ranger", n_top_genes=n_hvg,
                                batch_key=batch_key, subset=True)

    sc.pp.regress_out(ad, keys=["S_score", "G2M_score"])

    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=n_pcs, svd_solver="arpack", random_state=rng)
    return ad



