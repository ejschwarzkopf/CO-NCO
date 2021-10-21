#!/bin/bash
#SBATCH -J AlignCall_P
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/AlignCall_P_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/AlignCall_P_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file

## Alignment and variant calling for S. uvarum population sequencing 
## To run on NCSU BRC with SLURM
## Chris Large and Caiti S. Heil
## Uses the recommended SNP calling pipeline from Samtools
## Modified version for the CO-NCO project - Enrique J. Schwarzkopf
## This is the both-parents version

# Load programs:

module load bwa/0.7.17

module load samtools/1.12

module load picard/2.25.6

module load gatk/4.2.0.0

N=1
SAMPLEFILE1=/home4/eschwar3/CO_NCO/1.scripts/2.aux/uva_SE_parent1.txt
SAMPLEFILE2=/home4/eschwar3/CO_NCO/1.scripts/2.aux/uva_PE_parent2.txt
SAMPLE1=$(awk 'NR ==  '$N' {print}' $SAMPLEFILE1) # Passed sample prefix (ex: Sample-01)
SAMPLE2=$(awk 'NR ==  '$N' {print}' $SAMPLEFILE2) # Passed sample prefix (ex: Sample-01)
SAMPLE=$(echo ${SAMPLE1}X${SAMPLE2})
DIR=/home4/eschwar3/CO_NCO
WORKDIR=/home4/eschwar3/CO_NCO/3.output/1.align_varcall/4.parents # Where files will be created
SEQDIR=/home4/eschwar3/CO_NCO/2.data/2.sequence_data/1.parents_fastq/s_uvarum_seq/ # Location of Fastqs
MODDIR=/usr/local/bin
JAVA=/home2/cheil/programs/jre1.8.0_231/bin #GATK only works with Java v1.8
SEQID=Bread_insectyeast
REF=/home4/eschwar3/CO_NCO/2.data/1.reference_genome/Sbay.ultrascaf.fasta # Reference genome
#REF=/home2/cheil/reference_seq/sacCer3.fasta
SCRIPTS=/home4/eschwar3/CO_NCO/1.scripts # Location of custom scripts
VCFDIR=${WORKDIR}/vcf
BCFTOOLS=/home4/eschwar3/bin/bin #Location of bcftools

cd ${WORKDIR}

# Parent 1
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${SAMPLE1}.1_1.fastq.gz > ${SAMPLE1}_R1.sam

(>&2 echo ***Samtools - View***)
samtools view -b ${SAMPLE1}_R1.sam -o ${SAMPLE1}_R1.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE1}_R1.bam -o ${SAMPLE1}_R1_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE1}_R1_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${SAMPLE1}_R1_sort.bam

# Remove intermediate files
#rm ${SAMPLE}_R1.sam
#rm ${SAMPLE}_R1.bam



mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${SAMPLE1}_R1_sort.bam \
        OUTPUT=${SAMPLE1}_comb_R1.MD.bam \
        METRICS_FILE=dup_metrics/${SAMPLE1}_comb_R1.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${SAMPLE1}_comb_R1.MD.bam \
        O=${SAMPLE1}_comb_R1.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${SAMPLE1} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${SAMPLE1}_comb_R1.RG.MD.bam \
        -o ${SAMPLE1}_comb_R1.RG.MD.sort.bam
samtools index ${SAMPLE1}_comb_R1.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REF} \
        -I ${WORKDIR}/${SAMPLE1}_comb_R1.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${SAMPLE1}_GATK_HaplotypeCaller.vcf


# Parent 2
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${SAMPLE2}.1_1.fastq.gz ${SEQDIR}/${SAMPLE2}.1_2.fastq.gz > ${SAMPLE2}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${SAMPLE2}_R1R2.sam -o ${SAMPLE2}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE2}_R1R2.bam -o ${SAMPLE2}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE2}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${SAMPLE2}_R1R2_sort.bam

# Remove intermediate files
rm ${SAMPLE2}_R1R2.sam
rm ${SAMPLE2}_R1R2.bam



mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${SAMPLE2}_R1R2_sort.bam \
        OUTPUT=${SAMPLE2}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${SAMPLE2}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${SAMPLE2}_comb_R1R2.MD.bam \
        O=${SAMPLE2}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${SAMPLE2} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${SAMPLE2}_comb_R1R2.RG.MD.bam \
        -o ${SAMPLE2}_comb_R1R2.RG.MD.sort.bam
samtools index ${SAMPLE2}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REF} \
        -I ${WORKDIR}/${SAMPLE2}_comb_R1R2.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${SAMPLE2}_GATK_HaplotypeCaller.vcf


# All together now

(>&2 echo ***GATK - CombineGVCFs***)
gatk CombineGVCFs \
   -R ${REF} \
   --variant ${VCFDIR}/${SAMPLE1}_GATK_HaplotypeCaller.vcf \
   --variant ${VCFDIR}/${SAMPLE2}_GATK_HaplotypeCaller.vcf \
   -O ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller.vcf

(>&2 echo ***GATK - GenotypeGVCFs***)
gatk GenotypeGVCFs \
        -R ${REF} \
        -V ${VCFDIR}/${SAMPLE}_GATK_HaplotypeCaller.vcf \
        -O ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs.vcf

(>&2 echo ***bcftools - filter - indel distance***)
${BCFTOOLS}/bcftools filter \
        --SnpGap 10 \
        -o ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_IndelFiltered.vcf \
        ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs.vcf

(>&2 echo ***GATK - SelectVariants***)
# Extract SNPs from call set.
gatk SelectVariants \
        -R ${REF} \
        -V ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_IndelFiltered.vcf \
        --select-type-to-include SNP \
        -O ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS.vcf

(>&2 echo ***GATK - VariantFiltration***)
# Apply filters to call sets.
gatk VariantFiltration \
        -R ${REF} \
        -V ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS.vcf \
        -filter "QUAL < 100" --filter-name "QUAL100" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS_Tagged.vcf

(>&2 echo ***GATK - SelectVariants***)
# Select only variants that have passed the filters.
gatk SelectVariants \
        -R ${REF} \
        -V ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS_Tagged.vcf \
        --selectExpressions 'vc.isNotFiltered()' \
        -O ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf

(>&2 echo ***GATK - VariantsToTable***)
# Variants to table
gatk VariantsToTable \
        -R $REF \
        -V ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf \
        -F CHROM -F POS -GF AD \
        -O ${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.table


# Remove intermediates
rm ${SAMPLE}_R1.bam.intervals
rm ${SAMPLE}_R1.MD.bam
rm ${SAMPLE}_R1.RG.MD.bam
rm ${SAMPLE}_R1.RG.MD.realign.bam
rm ${SAMPLE}_R1.RG.MD.realign.bai
rm ${SAMPLE}_R1_sort.bam
rm ${SAMPLE}_R1_sort.bam.bai
