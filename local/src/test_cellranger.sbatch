#!/bin/bash
#SBATCH --job-name=test_cellranger_tiny
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elena.grassi@ircc.it
#SBATCH --partition=global
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=job_name_%j.log
#SBATCH --error=job_name_%j.error
#SBATCH --mem-per-cpu=1024M
#SBATCH --workdir=/home/egrassi/ba/local/src/cellranger-2.2.0/
./cellranger testrun --id=tiny
