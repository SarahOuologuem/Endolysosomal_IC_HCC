#!/bin/bash
#SBATCH --job-name=run_featureCounts
#SBATCH --output=/dss/dsshome1/05/ru62tig2/2.Read_Counting/logs/%j.log
#SBATCH --error=/dss/dsshome1/05/ru62tig2/2.Read_Counting/logs/%j.log
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --time=0-00:30:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --partition=teramem_inter
#SBATCH --mail-user=sarah.ouologuem@tum.de
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

# debugging information
echo "Job started on $(date)"

source /dss/dsshome1/05/ru62tig2/miniconda3/etc/profile.d/conda.sh

conda activate FeatureCounts


cd /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results/
for file in ./*.sam; do featureCounts -p -T 8 -a ../Homo_sapiens.GRCh38.113.gtf  -o "/dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/2.Read_Counting/results/human/${file%.*}.txt" $file; done


cd /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/2.Read_Counting/human/results
for file in ./*.txt; do
    if [ -e tempfile ]; then
        paste  tempfile <(cut -f 7 $file) >tempfile2
        mv tempfile2 tempfile
    else
        cut -f 1,7 $file > tempfile
    fi
done
mv tempfile ../H3B_matrix.tsv


cd /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/
for file in ./*.sam; do featureCounts -p -T 8 -a ../Mus_musculus.GRCm39.113.gtf  -o "/dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/2.Read_Counting/mouse/${file%.*}.txt" $file; done


cd /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/2.Read_Counting/mouse
for file in ./*.txt; do
    if [ -e tempfile ]; then
        paste  tempfile <(cut -f 7 $file) >tempfile2
        mv tempfile2 tempfile
    else
        cut -f 1,7 $file > tempfile
    fi
done
mv tempfile ../RIL_matrix.tsv


