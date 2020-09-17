#!/bin/bash
sample=$1
script=/scriptpath/s03.0_SNPCalling_Pipeline.sh
inpath=/path/to/bam/file/01.bam
outpath=/path/to/printout/result/02.SNP
mkdir -p $outpath
cat $sample|while read sp
do
	echo "bash $script $inpath $sp $outpath" >$sp.SNP.tmp.sh
	qsub -cwd -l vf=15g,p=4 -V  $sp.SNP.tmp.sh
done
