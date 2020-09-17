#!/bin/bash
indir=/datb1/wangrui/Project/Human_Gonad/fastq_path
sample=$2
refpath=/datb1/wangrui/bin/10X_pipeline/refdata-cellranger-hg19-1.2.0/
echo "cellranger count --id=$sample --fastqs=$indir/$sample --transcriptome=$refpath --sample=$sample --cells=5000" > $sample.00.cellranger.tmp.sh
qsub -cwd -l vf=100g,p=20 -V $sample.00.cellranger.tmp.sh
