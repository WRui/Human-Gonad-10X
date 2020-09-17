#!/bin/bash
indir=$1  ## Path to mapped bam Result
sample=$2  ## one sample name per line
outdir=$3  ## StatInfo/02.SNP_Calling

ref=/datb1/wangrui/bin/10X_pipeline/refdata-cellranger-hg19-1.2.0/fasta/sep_chr/genome_new.fa
PICARD="java -Djava.io.tmpdir=/datg/wangrui/tmp -jar /datc/wangrui/software/picard-tools-1.130/picard.jar"
GATK="java -Djava.io.tmpdir=/datg/wangrui/tmp -jar /datc/wangrui/software/GenomeAnalysisTK/GenomeAnalysisTK.jar"
knownSites_1=/datc/wangrui/Database/dbSNP/dbsnp_135.hg19.vcf
knownSites_2=/datc/wangrui/Database/indel_annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf
knownSites_3=/datc/wangrui/Database/dbSNP/1000G_phase1.indels.hg19.vcf
vcf=/datc/wangrui/Database/indel_annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf


function do_Remove_ERCC {
inpath=$1
#spike_in=/Share/BP/wangrui/Database/GATK/ERCC_RGC.list
spike_in=/WPSnew/wangrui/Database/Share_Database/Database/GATK/ERCC_RGC.list
sp=$2
outpath=$3

mkdir -p $outpath/$sp
cat $spike_in| grep -f - -v  <(samtools view -h $inpath/$sp/accepted_hits.bam) >  $outpath/$sp/$sp.NPI.sam 

}

function do_ReorderSam {
inpath=$1
PICARD=$PICARD
sp=$2
ref=$ref
outpath=$inpath

$PICARD ReorderSam I=$inpath/$sp/$sp.NPI.sam O=$outpath/$sp/$sp.reorder.sam REFERENCE=$ref 

}


function do_samTobam {
inpath=$1
outpath=$inpath
sp=$2
samtools view -bS $inpath/$sp/$sp.reorder.sam -o $outpath/$sp/$sp.reorder.bam
}

function do_10X_RG_added_sorted {

inpath=$1
PICARD=$PICARD
sp=$2
outpath=$3
mkdir -p $outpath/$sp
$PICARD AddOrReplaceReadGroups I=$inpath/$sp/$sp.reorder.bam O=$outpath/$sp/$sp.rg_added_sorted.bam SO=coordinate RGID=$sp RGLB=$sp RGPL=illumina RGPU=$sp RGSM=$sp

}

function do_Extract_chromosome {

inpath=$1
sp=$2
outpath=$3
mkdir -p $outpath/$sp

samtools view -b $inpath/$sp/possorted_genome_bam.bam 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT >$outpath/$sp/$sp.new.bam

}

function do_reheader_10X {
inpath=$1
PICARD=$PICARD
sp=$2
ref=$ref
outpath=$3

	samtools view -H $outpath/$sp/$sp.new.bam| grep -v '@HD' | grep -v '@SQ' > $outpath/$sp/$sp.tmp_header
	cat /datb1/wangrui/Project/Human_Gonad/StatInfo/s03.2_SNPCalling_GATK/02.SNP/new_header $outpath/$sp/$sp.tmp_header > $outpath/$sp/$sp.new_header
	samtools reheader $outpath/$sp/$sp.new_header $outpath/$sp/$sp.new.bam >$inpath/$sp/$sp.new_reheader.bam
}

function do_ReorderSam_10X {
inpath=$1
PICARD=$PICARD
sp=$2
ref=$ref
outpath=$3

#$PICARD ReorderSam I=$inpath/$sp/$sp.new.bam O=$outpath/$sp/$sp.reorder.bam SORT_ORDER=coordinate #REFERENCE=$ref
$PICARD ReorderSam I=$inpath/$sp/$sp.new_reheader.bam O=$outpath/$sp/$sp.reorder.bam REFERENCE=$ref

}

function do_ReorderSam_10X_2 {
inpath=$1
PICARD=$PICARD
sp=$2
ref=$ref
outpath=$inpath
$PICARD ReorderSam I=$outpath/$sp/$sp.dedupped.bam O=$outpath/$sp/$sp.reorder_V2.bam SORT_ORDER=coordinate
#REFERENCE=$ref

}


function do_RG_added_sorted {
inpath=$1
#PICARD=$PICARD
PICARD="java -Xmx20g -Djava.io.tmpdir=/datg/wangrui/tmp -jar /datc/wangrui/software/picard-tools-1.130/picard.jar"

sp=$2
outpath=$inpath

$PICARD AddOrReplaceReadGroups I=$inpath/$sp/$sp.reorder.bam O=$outpath/$sp/$sp.rg_added_sorted.bam SO=coordinate RGID=$sp RGLB=$sp RGPL=illumina RGPU=$sp RGSM=$sp

}

function do_mark_dup {
path=$1
PICARD=$PICARD
sp=$2
		$PICARD MarkDuplicates I=$path/$sp/$sp.rg_added_sorted.bam O=$path/$sp/$sp.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$path/$sp/$sp.metrics
}

function do_splitNTrim {
path=$1
GATK=$GATK
sp=$2
ref=$ref
		$GATK -T SplitNCigarReads -R $ref -I $path/$sp/$sp.dedupped.bam -o $path/$sp/$sp.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

}

function do_RealignerTargetCreator {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2

		$GATK -T RealignerTargetCreator -R $ref -I $inpath/$sp/$sp.split.bam -o $outpath/$sp/$sp.forIndelRealigner.intervals --known $vcf
}


function do_IndelRealigner {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2


		$GATK -T IndelRealigner -R $ref -I $inpath/$sp/$sp.split.bam -targetIntervals $outpath/$sp/$sp.forIndelRealigner.intervals -o $outpath/$sp/$sp.realignedBam.bam -known $vcf

}

function do_BaseRecalibrator {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
knownSites_1=$knownSites_1
knownSites_2=$knownSites_2
knownSites_3=$knownSites_3
		$GATK -T BaseRecalibrator -R $ref -I $inpath/$sp/$sp.realignedBam.bam -knownSites $knownSites_1 -knownSites $knownSites_2 -knownSites $knownSites_3 -o $outpath/$sp/$sp.recal.grp
}

function do_BaseRecalibrator_BQSR {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
knownSites_1=$knownSites_1
knownSites_2=$knownSites_2
knownSites_3=$knownSites_3
		$GATK -T BaseRecalibrator -R $ref -I $inpath/$sp/$sp.realignedBam.bam -knownSites $knownSites_1 -knownSites $knownSites_2 -knownSites $knownSites_3 -BQSR $inpath/$sp/$sp.recal.grp -o $outpath/$sp/$sp.pos_recal.grp 
}

function do_PrintReads {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
knownSites_1=$knownSites_1
knownSites_2=$knownSites_2
knownSites_3=$knownSites_3
		$GATK -T PrintReads -R $ref -I $inpath/$sp/$sp.realignedBam.bam -BQSR $inpath/$sp/$sp.recal.grp -o $outpath/$sp/$sp.realn_Recal.bam

}


function do_HaplotypeCaller {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
		$GATK -T HaplotypeCaller -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $outpath/$sp/$sp.variants_result.vcf
}

function do_HaplotypeCaller_GVCF {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
        $GATK -T HaplotypeCaller -R $ref -I $inpath/$sp/$sp.realn_Recal.bam -G StandardAnnotation --emitRefConfidence GVCF -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -variant_index_type LINEAR -variant_index_parameter 128000 -o $outpath/$sp/$sp.variants_result.AS.g.vcf
}

function do_VarintFilter {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
ref=$ref
		$GATK -T VariantFiltration -R $ref -V $inpath/$sp/$sp.variants_result.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $outpath/$sp/$sp.filtered_variants_result.vcf
}


function do_SnpEFF_annatation {
inpath=$1
outpath=$inpath
GATK=$GATK
sp=$2
snpEff=/datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.jar
		java -Xmx10g -jar $snpEff -c /datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.filtered_variants_result.vcf >$outpath/${sp}/${sp}.ann.vcf

}

function do_SnpEff_annotation_gVCF {
inpath=$1
outpath=$inpath
sp=$2
snpEff=/datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.jar

java -Xmx10g -jar $snpEff -c /datg/wangrui/Software_NewVersion/snpEff_v4.3/snpEff/snpEff.config hg19 $inpath/${sp}/${sp}.variants_result.AS.g.vcf > $outpath/${sp}/${sp}.ann.g.vcf
}

#do_Remove_ERCC $indir  $sample $outdir #step 00
#do_ReorderSam $outdir $sample          #step 01
#do_samTobam $outdir $sample            #step 02
#do_RG_added_sorted $outdir $sample     #step 03
do_Extract_chromosome $indir $sample $outdir #step add for 10X
#do_reheader_10X $outdir $sample $outdir
do_ReorderSam_10X $outdir $sample $outdir #step add for 10X
#do_10X_RG_added_sorted $outdir $sample $outdir   #step 03
do_mark_dup $outdir $sample            #step 04
#do_ReorderSam_10X_2 $outdir $sample
do_splitNTrim $outdir $sample          #step 05
do_RealignerTargetCreator $outdir $sample #step 06
do_IndelRealigner $outdir $sample      #step 07
do_BaseRecalibrator $outdir $sample    #step 08
do_BaseRecalibrator_BQSR $outdir $sample #step 09
do_PrintReads $outdir $sample          #step 10
do_HaplotypeCaller $outdir $sample     #step 11
do_VarintFilter $outdir $sample        #step 12
do_SnpEFF_annatation $outdir $sample   #step 13

#Option
do_HaplotypeCaller_GVCF $outdir $sample
do_SnpEff_annotation_gVCF $outdir $sample
