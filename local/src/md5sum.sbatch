#!/bin/bash
#SBATCH --job-name=md5sum
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elena.grassi@ircc.it
#SBATCH --partition=global
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=job_name_%j.log
#SBATCH --error=job_name_%j.error
#SBATCH --mem-per-cpu=1024M
#SBATCH --workdir=/home/egrassi/ba/local/share/data
md5sum /home/egrassi/ba/local/share/data/hpcnfs/techunits/genomics/PublicData/Pasini/dmangana/FASTQ/180717_A00302_0025_AHFKJTDMXX/*/* > md5sum_pasini
