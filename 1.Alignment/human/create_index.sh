#!/bin/bash
#SBATCH --job-name=createindex_human
#SBATCH --output=/dss/dsshome1/05/ru62tig2/1.Alignment/logs/%j.log
#SBATCH --error=/dss/dsshome1/05/ru62tig2/1.Alignment/logs/%j.log
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --time=1-00:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=18
#SBATCH --partition=teramem_inter
#SBATCH --mail-user=sarah.ouologuem@tum.de
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

# debugging information
echo "Job started on $(date)"

source /dss/dsshome1/05/ru62tig2/miniconda3/etc/profile.d/conda.sh

conda activate STAR

# downloaded gtf from: https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
# downloaded fasta from: https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

STAR --runThreadN 18 \
--runMode genomeGenerate \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--genomeFastaFiles /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/Homo_sapiens.GRCh38.113.gtf \
--sjdbOverhang 149
