#!/bin/bash
#SBATCH --job-name=runfastqc
#SBATCH --output=/dss/dsshome1/05/ru62tig2/0.Quality_Control/1.FastQC/logs/%j.log
#SBATCH --error=/dss/dsshome1/05/ru62tig2/0.Quality_Control/1.FastQC/logs/%j.log
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

conda activate fastqc

cd /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/X208SC24103610-Z01-F001/01.RawData/

for file in ./*/*.fq.gz; do fastqc $file -o /dss/dsshome1/05/ru62tig2/0.Quality_Control/1.FastQC/outputs/; done

