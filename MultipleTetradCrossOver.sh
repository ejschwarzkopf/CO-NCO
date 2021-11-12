#!/bin/bash
#SBATCH -J AlignCall_P
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/AlignCall_P_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/AlignCall_P_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file

while getopts parent1:parent2:se1:se2: flag
do
	case "${flag}" in
		parent1) PARENT1=${OPRARG};;
		parent2) PARENT2=${OPTARG};;
		se1) 

# Load programs:

module load bwa/0.7.17

module load samtools/1.12

module load picard/2.25.6

module load gatk/4.2.0.0

#Think about inputs for a second

#Parent names, that lead to fastq files

#Tetrad names, that also lead to fastq files

#Idea: start a naming convention
#	Parent1XParent2_SampleDetails.fastq
#	This lets me cut -d '_' 1 | cut -d 'X' 1,2 to get the parents and use grep to ls the offspring

# Some file preparation
# For the American parents:
# for filename in `ls YMD32343X3245_*`; do newname=$(echo 'SRR1119189XSRR1119180_'$(echo $filename | cut -d '_' -f 2-)); mv $filename $newname; done
# For the European parents:
# for filename in `ls CSH35XPYCC6876_*`; do newname=$(echo 'SRR1119200XSRR1119199_'$(echo $filename | cut -d '_' -f 2-)); mv $filename $newname; done

#### First, we genotype the parents ####

N=1
PARENT1=SRR1119180 # Passed sample prefix (ex: Sample-01)
PARENT2=SRR1119189 # Passed sample prefix (ex: Sample-01)
SAMPLE=$(echo ${PARENT1}X${PARENT2})
DIR=/home4/eschwar3/CO_NCO
WORKDIR=/home4/eschwar3/CO_NCO/3.output/1.align_varcall/4.parents/$SAMPLE # Where files will be created
mkdir $WORKDIR
SEQDIR=/home4/eschwar3/CO_NCO/2.data/2.sequence_data/1.parents_fastq/s_uvarum_seq/ # Location of Fastqs
MODDIR=/usr/local/bin
JAVA=/home2/cheil/programs/jre1.8.0_231/bin #GATK only works with Java v1.8
#SEQID=Bread_insectyeast
REF=/home4/eschwar3/CO_NCO/2.data/1.reference_genome/Sbay.ultrascaf.fasta # Reference genome
#REF=/home2/cheil/reference_seq/sacCer3.fasta
SCRIPTS=/home4/eschwar3/CO_NCO/1.scripts # Location of custom scripts
VCFDIR=${WORKDIR}/vcf
mkdir $VCFDIR
BCFTOOLS=/home4/eschwar3/bin/bin #Location of bcftools

cd ${WORKDIR}

# Parent 1
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${PARENT1}.1_1.fastq.gz > ${PARENT1}_R1.sam

(>&2 echo ***Samtools - View***)
samtools view -b ${PARENT1}_R1.sam -o ${PARENT1}_R1.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${PARENT1}_R1.bam -o ${PARENT1}_R1_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${PARENT1}_R1_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${PARENT1}_R1_sort.bam

# Remove intermediate files
#rm ${SAMPLE}_R1.sam
#rm ${SAMPLE}_R1.bam



mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${PARENT1}_R1_sort.bam \
        OUTPUT=${PARENT1}_comb_R1.MD.bam \
        METRICS_FILE=dup_metrics/${PARENT1}_comb_R1.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${PARENT1}_comb_R1.MD.bam \
        O=${PARENT1}_comb_R1.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${PARENT1} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${PARENT1}_comb_R1.RG.MD.bam \
        -o ${PARENT1}_comb_R1.RG.MD.sort.bam
samtools index ${PARENT1}_comb_R1.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REF} \
        -I ${WORKDIR}/${PARENT1}_comb_R1.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${PARENT1}_GATK_HaplotypeCaller.vcf

# Parent 2
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem ${REF} ${SEQDIR}/${PARENT2}.1_1.fastq.gz ${SEQDIR}/${PARENT2}.1_2.fastq.gz > ${PARENT2}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${PARENT2}_R1R2.sam -o ${PARENT2}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${PARENT2}_R1R2.bam -o ${PARENT2}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${PARENT2}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${PARENT2}_R1R2_sort.bam

# Remove intermediate files
rm ${PARENT2}_R1R2.sam
rm ${PARENT2}_R1R2.bam



mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} MarkDuplicates \
        INPUT=${PARENT2}_R1R2_sort.bam \
        OUTPUT=${PARENT2}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${PARENT2}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
${JAVA}/java -Xmx2g -jar ${PICARD_JAR} AddOrReplaceReadGroups \
        I=${PARENT2}_comb_R1R2.MD.bam \
        O=${PARENT2}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${PARENT2} \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${PARENT2}_comb_R1R2.RG.MD.bam \
        -o ${PARENT2}_comb_R1R2.RG.MD.sort.bam
samtools index ${PARENT2}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - HaplotypeCaller***)
gatk HaplotypeCaller \
        -R ${REF} \
        -I ${WORKDIR}/${PARENT2}_comb_R1R2.RG.MD.sort.bam \
        -ERC GVCF \
        -O ${VCFDIR}/${PARENT2}_GATK_HaplotypeCaller.vcf

### The ${VCFDIR}/${PARENT1}_GATK_HaplotypeCaller.vcf and ${VCFDIR}/${PARENT2}_GATK_HaplotypeCaller.vcf files are the ones I need for the joint genotype calling for the tetrads

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




Now I want to separate the tetrads

A line like: 

# ls SRR1119189XSRR1119180_* | cut -d '_' -f 2,3 | sed 's/.$//' | uniq 

will let me see all the tetrad groups. Some of these will be one tetrad (4 files) and some will be two (8 files)

I can separate them by making an if statement where if SR...$lineG* exists, then this set has two tetrads and I treat it one way, else it has one tetrad and I treat it differently