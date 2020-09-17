#!/bin/bash
infile=$1
sample=$2
outdir=$3

echo "samtools idxstats $infile |grep -v 'GL00' |grep -v '*' |sort -k 1.4n > $outdir/${sample}.chr_cor.txt" >${sample}.tmp.sh
qsub -cwd -l vf=1g,h=tanglab2,p=2 -V ${sample}.tmp.sh
