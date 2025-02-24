#!/bin/bash
#SBATCH --job-name=align_mouse
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

cd /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/X208SC24103610-Z01-F001/01.RawData/


STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_dKO_1/RIL_dKO_1_MKRN240023756-1A_22HGLTLT4_L1_1.fq.gz RIL_dKO_1/RIL_dKO_1_MKRN240023756-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_dKO_1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_M_1/RIL_M_1_MKRN240023753-1A_22HGLTLT4_L1_1.fq.gz RIL_M_1/RIL_M_1_MKRN240023753-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_M_1_L1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_M_1/RIL_M_1_MKRN240023753-1A_22HGGMLT4_L7_1.fq.gz RIL_M_1/RIL_M_1_MKRN240023753-1A_22HGGMLT4_L7_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_M_1_L7_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_T_1/RIL_T_1_MKRN240023750-1A_22HGLTLT4_L3_1.fq.gz RIL_T_1/RIL_T_1_MKRN240023750-1A_22HGLTLT4_L3_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_T_1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_WT_1/RIL_WT_1_MKRN240023747-1A_22HGLTLT4_L1_1.fq.gz RIL_WT_1/RIL_WT_1_MKRN240023747-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_WT_1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_dKO_2/RIL_dKO_2_MKRN240023757-1A_22HGLTLT4_L1_1.fq.gz RIL_dKO_2/RIL_dKO_2_MKRN240023757-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_dKO_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_M_2/RIL_M_2_MKRN240023754-1A_22HGLTLT4_L2_1.fq.gz RIL_M_2/RIL_M_2_MKRN240023754-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_M_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_T_2/RIL_T_2_MKRN240023751-1A_22HGLTLT4_L1_1.fq.gz RIL_T_2/RIL_T_2_MKRN240023751-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_T_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_WT_2/RIL_WT_2_MKRN240023748-1A_22HGLTLT4_L1_1.fq.gz RIL_WT_2/RIL_WT_2_MKRN240023748-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_WT_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_dKO_3/RIL_dKO_3_MKRN240023758-1A_22HGLTLT4_L2_1.fq.gz RIL_dKO_3/RIL_dKO_3_MKRN240023758-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_dKO_3_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_M_3/RIL_M_3_MKRN240023755-1A_22HGLTLT4_L5_1.fq.gz RIL_M_3/RIL_M_3_MKRN240023755-1A_22HGLTLT4_L5_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_M_3_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_T_3/RIL_T_3_MKRN240023752-1A_22HGLTLT4_L1_1.fq.gz RIL_T_3/RIL_T_3_MKRN240023752-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_T_3_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/index \
--readFilesIn RIL_WT_3/RIL_WT_3_MKRN240023749-1A_22HGLTLT4_L1_1.fq.gz RIL_WT_3/RIL_WT_3_MKRN240023749-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/mouse/results/RIL_WT_3_

