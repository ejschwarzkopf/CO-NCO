#!/bin/bash

#SBATCH -J AlignCall_O
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/AlignCall_O_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/AlignCall_O_%J.out
##SBATCH --array=1-4 # Num of jobs = Num of rows in sample list file

## Alignment and variant calling for S. uvarum population sequencing 
## To run on NCSU BRC with SLURM
## Chris Large and Caiti S. Heil
## Uses the recommended SNP calling pipeline from Samtools
## Modified version for the CO-NCO project - Enrique J. Schwarzkopf


# Load programs:

module load bwa/0.7.17

module load samtools/1.12

module load picard/2.25.6

module load gatk/4.2.0.0

N=${SLURM_ARRAY_TASK_ID}
SAMPLEFILE=/home4/eschwar3/CO_NCO/1.scripts/2.aux/offspring.txt
SAMPLE1=$(awk 'NR ==  1 {print}' $SAMPLEFILE) # Passed sample prefix (ex: Sample-01)
SAMPLE2=$(awk 'NR ==  2 {print}' $SAMPLEFILE) # Passed sample prefix (ex: Sample-01)
SAMPLE3=$(awk 'NR ==  3 {print}' $SAMPLEFILE) # Passed sample prefix (ex: Sample-01)
SAMPLE4=$(awk 'NR ==  4 {print}' $SAMPLEFILE) # Passed sample prefix (ex: Sample-01)
SAMPLE=CSH345X347-All
DIR=/home4/eschwar3/CO_NCO
WORKDIR=/home4/eschwar3/CO_NCO/3.output/1.align_varcall/3.offspring # Where files will be created
SEQDIR=/home4/eschwar3/CO_NCO/2.data/2.sequence_data/2.offspring_seq/BaseCalls/ # Location of Fastqs
MODDIR=/usr/local/bin
JAVA=/home2/cheil/programs/jre1.8.0_231/bin #GATK only works with Java v1.8
SEQID=Bread_insectyeast
REF=/home4/eschwar3/CO_NCO/2.data/1.reference_genome/Sbay.ultrascaf.fasta # Reference genome
#REF=/home2/cheil/reference_seq/sacCer3.fasta
SCRIPTS=/home4/eschwar3/CO_NCO/1.scripts # Location of custom scripts
VCFDIR=${WORKDIR}/vcf
BCFTOOLS=/home4/eschwar3/bin/bin #Location of bcftools

cd ${WORKDIR}

# Tetrad 1
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${SAMPLE1}_1.fastq.gz ${SEQDIR}/${SAMPLE1}_2.fastq.gz > ${SAMPLE1}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${SAMPLE1}_R1R2.sam -o ${SAMPLE1}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE1}_R1R2.bam -o ${SAMPLE1}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE1}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${SAMPLE1}_R1R2_sort.bam

# Remove intermediate files
rm ${SAMPLE1}_R1R2.sam
rm ${SAMPLE1}_R1R2.bam


#make directory for duplicate calls 
mkdir -p ${WORKDIR}/dup_metrics #no error if existing

(>&2 echo ***Picard - MarkDuplicates***)
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${SAMPLE1}_R1R2_sort.bam \
        OUTPUT=${SAMPLE1}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${SAMPLE1}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${SAMPLE1}_comb_R1R2.MD.bam \
        O=${SAMPLE1}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${SAMPLE1} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${SAMPLE1}_comb_R1R2.RG.MD.bam \
        -o ${SAMPLE1}_comb_R1R2.RG.MD.sort.bam
samtools index ${SAMPLE1}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REF} \
        -I ${WORKDIR}/${SAMPLE1}_comb_R1R2.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${SAMPLE1}_GATK_HaplotypeCaller.vcf

# Tetrad 2
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${SAMPLE2}_1.fastq.gz ${SEQDIR}/${SAMPLE2}_2.fastq.gz > ${SAMPLE2}_R1R2.sam

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


#make directory for duplicate calls 
mkdir -p ${WORKDIR}/dup_metrics #no error if existing

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

# Tetrad 3
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${SAMPLE3}_1.fastq.gz ${SEQDIR}/${SAMPLE3}_2.fastq.gz > ${SAMPLE3}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${SAMPLE3}_R1R2.sam -o ${SAMPLE3}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE3}_R1R2.bam -o ${SAMPLE3}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE3}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${SAMPLE3}_R1R2_sort.bam

# Remove intermediate files
rm ${SAMPLE3}_R1R2.sam
rm ${SAMPLE3}_R1R2.bam


#make directory for duplicate calls 
mkdir -p ${WORKDIR}/dup_metrics #no error if existing

(>&2 echo ***Picard - MarkDuplicates***)
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${SAMPLE3}_R1R2_sort.bam \
        OUTPUT=${SAMPLE3}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${SAMPLE3}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${SAMPLE3}_comb_R1R2.MD.bam \
        O=${SAMPLE3}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${SAMPLE3} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${SAMPLE3}_comb_R1R2.RG.MD.bam \
        -o ${SAMPLE3}_comb_R1R2.RG.MD.sort.bam
samtools index ${SAMPLE3}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REF} \
        -I ${WORKDIR}/${SAMPLE3}_comb_R1R2.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${SAMPLE3}_GATK_HaplotypeCaller.vcf

# Tetrad 4
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${SAMPLE4}_1.fastq.gz ${SEQDIR}/${SAMPLE4}_2.fastq.gz > ${SAMPLE4}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${SAMPLE4}_R1R2.sam -o ${SAMPLE4}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE4}_R1R2.bam -o ${SAMPLE4}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE4}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${SAMPLE4}_R1R2_sort.bam

# Remove intermediate files
rm ${SAMPLE4}_R1R2.sam
rm ${SAMPLE4}_R1R2.bam


#make directory for duplicate calls 
mkdir -p ${WORKDIR}/dup_metrics #no error if existing

(>&2 echo ***Picard - MarkDuplicates***)
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${SAMPLE4}_R1R2_sort.bam \
        OUTPUT=${SAMPLE4}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${SAMPLE4}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${SAMPLE4}_comb_R1R2.MD.bam \
        O=${SAMPLE4}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${SAMPLE4} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${SAMPLE4}_comb_R1R2.RG.MD.bam \
        -o ${SAMPLE4}_comb_R1R2.RG.MD.sort.bam
samtools index ${SAMPLE4}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REF} \
        -I ${WORKDIR}/${SAMPLE4}_comb_R1R2.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${SAMPLE4}_GATK_HaplotypeCaller.vcf

# All together now

gatk CombineGVCFs \
	-R ${REF} \
	--variant ${VCFDIR}/${SAMPLE1}_GATK_HaplotypeCaller.vcf \
	--variant ${VCFDIR}/${SAMPLE2}_GATK_HaplotypeCaller.vcf \
	--variant ${VCFDIR}/${SAMPLE3}_GATK_HaplotypeCaller.vcf \
	--variant ${VCFDIR}/${SAMPLE4}_GATK_HaplotypeCaller.vcf \
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
        -filter "QUAL < 30" --filter-name "QUAL30" \
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
rm ${SAMPLE}_R1R2.bam.intervals
rm ${SAMPLE}_R1R2.MD.bam
rm ${SAMPLE}_R1R2.RG.MD.bam
rm ${SAMPLE}_R1R2.RG.MD.realign.bam
rm ${SAMPLE}_R1R2.RG.MD.realign.bai
rm ${SAMPLE}_R1R2_sort.bam
rm ${SAMPLE}_R1R2_sort.bam.bai

