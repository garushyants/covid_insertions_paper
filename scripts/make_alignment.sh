#!/bin/bash

conda activate lofreq

SRA_IDs_list='SRA_IDs_list' #file that contains column of SRA IDs that have to be verified

IFS=' ' read -r -a alltasks <<< $(less $SRA_IDs_list)
PREFIX=${alltasks[$SGE_TASK_ID-1]}
BAM=${PREFIX}_sorted.bam

echo "$PREFIX"
if [ -f "${PREFIX}_1.fastq" ]; then
 bowtie2 -x SARS-CoV-2_ref -1 ${PREFIX}_1.fastq -2 ${PREFIX}_2.fastq -S ${PREFIX}.sam
else
 bowtie2 -x SARS-CoV-2_ref -U ${PREFIX}.fastq -S ${PREFIX}.sam
fi
samtools view -b -S -o ${PREFIX}.bam ${PREFIX}.sam
samtools sort -o $BAM ${PREFIX}.bam
samtools index $BAM

###Start LoFreq
LOFREQ=lofreq
REF=SARS-CoV2_reference_NC_045512.fasta
VITERBAM=${PREFIX}_viterbi.bam
SORTED=${PREFIX}_viterbi_sorted.bam
IQBAM=${PREFIX}_viterbi_iq.bam
VCF=${PREFIX}.vcf
SNVS=${PREFIX}.dp4.SNVs.vcf
INDELS=${PREFIX}.dp10.indels.vcf

$LOFREQ viterbi -f $REF -o $VITERBAM $BAM
samtools sort -o $SORTED $VITERBAM
$LOFREQ indelqual -f $REF -o $IQBAM --dindel $SORTED
samtools index $IQBAM

#call variants using  min depth 2
$LOFREQ call -f $REF -o $VCF  -C 2 --call-indels --no-default-filter --force-overwrite --use-orphan $IQBAM

#call indels
$LOFREQ filter  -Q 20 -K 20 --no-defaults  -v 4 -V 0 -a 0.500001 -A 0 --only-snvs -i $VCF -o $SNVS

#do indels with DP >= 20 and AF >= 0.5
$LOFREQ filter  -Q 20 -K 20 --no-defaults  -v 10 -V 0 -a 0.50000 -A 0 --only-indels -i  $VCF -o $INDELS
