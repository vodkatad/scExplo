#!/usr/bin/env python3

import glob
import os
import sys

link = glob.glob(sys.argv[1] +'_seurat_clusters_*.csv')
if len(link) == 1:
    linked = os.readlink(link[0])
    print(os.path.realpath(linked))
else:
    raise(Exception('step_04 did not generate what is needed :(' + sys.argv[1] + '_seurat_clusters_*.csv'))