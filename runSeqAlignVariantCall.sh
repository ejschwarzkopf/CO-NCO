#!/bin/bash

#SBATCH -J seqalign
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out

## Alignment and variant calling for S. uvarum population sequencing 
## To run on NCSU BRC with SLURM
## Chris Large and Caiti S. Heil
## Uses the recommended SNP calling pipeline from Samtools
## Then filters based on the ancestral sequence


SAMPLE=$1 # Passed sample prefix (ex: Sample-01)
DIR=/home2/cheil
WORKDIR=/home2/cheil/s_uvarum # Where files will be created
SEQDIR=/home2/cheil/s_uvarum/fastq/s_uvarum_seq # Location of Fastqs
MODDIR=/usr/local/bin
JAVA=/home2/cheil/programs/jre1.8.0_231/bin #GATK only works with Java v1.8
SEQID=Suva_popseq
REF=/home2/cheil/reference_seq/Sbay.ultrascaf.fasta # Reference genome
#ANNOTATE=/net/dunham/vol2/Archives/Archived_UW_people/Anna_Sunshine/dunham/sequencing_analysis/SNPs/SNVs_paired_end # Location of custom annotation scripts
SCRIPTS=/home2/cheil/scripts # Location of custom scripts
VCFDIR=${WORKDIR}/vcf

cd ${WORKDIR}

#
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${SAMPLE}_1.fastq > ${SAMPLE}_R1.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${SAMPLE}_R1.sam -o ${SAMPLE}_R1.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE}_R1.bam -o ${SAMPLE}_R1_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE}_R1_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${SAMPLE}_R1_sort.bam

# Remove intermediate files
rm ${SAMPLE}_R1.sam
rm ${SAMPLE}_R1.bam



mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
java -Xmx2g -jar ${MODDIR}/picard.jar MarkDuplicates \
        INPUT=${SAMPLE}_R1_sort.bam \
        OUTPUT=${SAMPLE}_comb_R1.MD.bam \
        METRICS_FILE=dup_metrics/${SAMPLE}_comb_R1.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${MODDIR}/picard.jar AddOrReplaceReadGroups \
        I=${SAMPLE}_comb_R1.MD.bam \
        O=${SAMPLE}_comb_R1.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${SAMPLE} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${SAMPLE}_comb_R1.RG.MD.bam \
        -o ${SAMPLE}_comb_R1.RG.MD.sort.bam
samtools index ${SAMPLE}_comb_R1.RG.MD.sort.bam

#(>&2 echo ***GATK - RealingerTargetCreator***)
##GATK Realinger
${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${REF} \
        -I ${SAMPLE}_comb_R1.RG.MD.sort.bam \
        -o ${SAMPLE}_comb_R1.bam.intervals
${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ${REF} \
        -I ${SAMPLE}_comb_R1.RG.MD.sort.bam \
        -targetIntervals ${SAMPLE}_comb_R1.bam.intervals \
        -o ${SAMPLE}_comb_R1.RG.MD.realign.bam
#
samtools sort ${SAMPLE}_comb_R1.RG.MD.realign.bam \
        -o ${SAMPLE}_comb_R1.RG.MD.realign.sort.bam
samtools index ${SAMPLE}_comb_R1.RG.MD.realign.sort.bam
#
${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
	-R ${REF} \
	-T HaplotypeCaller \
	-I ${WORKDIR}/${SAMPLE}_comb_R1.RG.MD.realign.sort.bam \
	-o ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller.vcf

# Extract SNPs from call set.
${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ${REF} \
	-V ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller.vcf \
	-selectType SNP \
	-o ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller_SNPS.vcf

# Apply filters to call sets.
${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R ${REF} \
	-V ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller_SNPS.vcf \
	--filterExpression \
	"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "GATK_Best_Practices_default" \
	-o ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller_SNPS_Tagged.vcf

# Select only variants that have passed the filters.
${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ${REF} \
	-V ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller_SNPS_Tagged.vcf \
	-select 'vc.isNotFiltered()' \
	-o ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller_SNPS_Tagged_Filtered.vcf

# Variants to table
${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
        -R $REF \
        -T VariantsToTable \
        -V ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller_SNPS_Tagged_Filtered.vcf \
        -F CHROM -F POS -GF AD \
        -o ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller_SNPS_Tagged_Filtered.table


# Remove intermediates
rm ${SAMPLE}_R1.bam.intervals
rm ${SAMPLE}_R1.MD.bam
rm ${SAMPLE}_R1.RG.MD.bam
rm ${SAMPLE}_R1.RG.MD.realign.bam
rm ${SAMPLE}_R1.RG.MD.realign.bai
rm ${SAMPLE}_R1_sort.bam
rm ${SAMPLE}_R1_sort.bam.bai
