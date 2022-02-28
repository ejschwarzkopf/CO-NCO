#!/bin/bash

# Load programs:

module load bwa/0.7.17

module load samtools/1.12

module load picard/2.25.6

module load gatk/4.2.0.0

module load java/openjdk-11

SEQID=Yeast
JAVA=/home2/cheil/programs/jre1.8.0_231/bin #GATK only works with Java v1.8
PARENT=
DIR=
OUTPUTDIR=
REFERENCE=

while [ $# -gt 0 ] ; do
	case $1 in
		--parent) PARENT=$2 ;;
		--directory) DIR=$2 ;;
		--reference-genome) REFERENCE=$2 ;;
		--output-directory) OUTPUTDIR=$2 ;;
	esac

	shift
done

VCFDIR=${OUTPUTDIR}/vcf

cd ${OUTPUTDIR}

mkdir ${VCFDIR}

# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REFERENCE} ${DIR}/${PARENT}.1_1.fastq.gz ${DIR}/${PARENT}.1_2.fastq.gz > ${PARENT}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${PARENT}_R1R2.sam -o ${PARENT}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${PARENT}_R1R2.bam -o ${PARENT}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${PARENT}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${OUTPUTDIR}/${PARENT}_R1R2_sort.bam

# Remove intermediate files
rm ${PARENT}_R1R2.sam
rm ${PARENT}_R1R2.bam



mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${PARENT}_R1R2_sort.bam \
        OUTPUT=${PARENT}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${PARENT}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${PARENT}_comb_R1R2.MD.bam \
        O=${PARENT}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${PARENT} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${PARENT}_comb_R1R2.RG.MD.bam \
        -o ${PARENT}_comb_R1R2.RG.MD.sort.bam
samtools index ${PARENT}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REFERENCE} \
        -I ${OUTPUTDIR}/${PARENT}_comb_R1R2.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${PARENT}_GATK_HaplotypeCaller.vcf
