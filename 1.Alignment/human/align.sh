#!/bin/bash
#SBATCH --job-name=align_human
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
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_dKO_1/H3B_dKO_1_MKRN240023768-1A_22HGLTLT4_L2_1.fq.gz H3B_dKO_1/H3B_dKO_1_MKRN240023768-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_dKO_1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_M_1/H3B_M_1_MKRN240023765-1A_22HGLTLT4_L2_1.fq.gz H3B_M_1/H3B_M_1_MKRN240023765-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_M_1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_T_1/H3B_T_1_MKRN240023762-1A_22HGLTLT4_L2_1.fq.gz H3B_T_1/H3B_T_1_MKRN240023762-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_T_1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_WT_1/H3B_WT_1_MKRN240023759-1A_22HGLTLT4_L2_1.fq.gz H3B_WT_1/H3B_WT_1_MKRN240023759-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_WT_1_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_dKO_2/H3B_dKO_2_MKRN240023769-1A_22HGLTLT4_L2_1.fq.gz H3B_dKO_2/H3B_dKO_2_MKRN240023769-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_dKO_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_M_2/H3B_M_2_MKRN240023766-1A_22HGLTLT4_L1_1.fq.gz H3B_M_2/H3B_M_2_MKRN240023766-1A_22HGLTLT4_L1_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_M_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_T_2/H3B_T_2_MKRN240023763-1A_22HGLTLT4_L4_1.fq.gz H3B_T_2/H3B_T_2_MKRN240023763-1A_22HGLTLT4_L4_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_T_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_WT_2/H3B_WT_2_MKRN240023760-1A_22HGLTLT4_L2_1.fq.gz H3B_WT_2/H3B_WT_2_MKRN240023760-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_WT_2_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_dKO_3/H3B_dKO_3_MKRN240023770-1A_22HGLTLT4_L2_1.fq.gz H3B_dKO_3/H3B_dKO_3_MKRN240023770-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_dKO_3_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_M_3/H3B_M_3_MKRN240023767-1A_22HGLTLT4_L4_1.fq.gz H3B_M_3/H3B_M_3_MKRN240023767-1A_22HGLTLT4_L4_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_M_3_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_T_3/H3B_T_3_MKRN240023764-1A_22HGLTLT4_L4_1.fq.gz H3B_T_3/H3B_T_3_MKRN240023764-1A_22HGLTLT4_L4_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_T_3_

STAR --runThreadN 18 \
--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
--readFilesIn H3B_WT_3/H3B_WT_3_MKRN240023761-1A_22HGLTLT4_L2_1.fq.gz H3B_WT_3/H3B_WT_3_MKRN240023761-1A_22HGLTLT4_L2_2.fq.gz \
--readFilesCommand gunzip -c \
--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results2/H3B_WT_3_

#STAR --runThreadN 18 \
#--genomeDir /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/index \
#--readFilesManifest /dss/dsshome1/05/ru62tig2/1.Alignment/human/samples.tsv \
#--readFilesCommand gunzip -c \
#--outFileNamePrefix /dss/dssfs02/lwp-dss-0001/ui521/ui521-dss-0000/ru62tig2/1.Alignment/human/results/H3B_