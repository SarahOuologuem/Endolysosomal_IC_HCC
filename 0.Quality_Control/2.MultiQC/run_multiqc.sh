#!/bin/bash
#SBATCH --job-name=runmultiqc
#SBATCH --output=/dss/dsshome1/05/ru62tig2/0.Quality_Control/2.MultiQC/logs/%j.log
#SBATCH --error=/dss/dsshome1/05/ru62tig2/0.Quality_Control/2.MultiQC/logs/%j.log
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --time=8-00:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --partition=teramem_inter
#SBATCH --mail-user=sarah.ouologuem@tum.de
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

# debugging information
echo "Job started on $(date)"

source /dss/dsshome1/05/ru62tig2/miniconda3/etc/profile.d/conda.sh

conda activate multiqc

multiqc /dss/dsshome1/05/ru62tig2/0.Quality_Control/1.FastQC/outputs -o /dss/dsshome1/05/ru62tig2/0.Quality_Control/2.MultiQC/outputs

