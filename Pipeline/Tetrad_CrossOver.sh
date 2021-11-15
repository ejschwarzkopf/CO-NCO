#!/bin/bash

# This script takes all the tetrads from a cross, aligns, calls variants, and makes the seg file for each set, then runs CrossOver on the all the tetrads.

PARENT1=
PARENT2=
PARENTDIR=
TETRADDIR=
REFERENCE=
OUTPUTDIR=
CROSSOVERDIR=

while :; do
	case $1 in
		--parent1) PARENT1=$2 ;;
		--parent2) PARENT2=$2 ;;
		--parent-directory) PARENTDIR=$2 ;;
		--tetrad-directory) TETRADDIR=$2 ;;
		--reference-genome) REFERENCE=$2 ;;
		--output-directory) OUTPUTDIR=$2 ;;
		--crossover-directory) CROSSOVERDIR=$2 ;;
	esac

	shift
done

cd ${OUTPUTDIR}

PARENTS=${PARENT1}X${PARENT2}


mkdir ${OUTPUTDIR}/segfiles
mkdir ${OUTPUTDIR}/out

TETRAD_COUNT=$(ls ${TETRADDIR}/${PARENTS}_*A_R1.fastq.gz | wc -l)

touch tmp_tetradnames_${PARENTS}.txt

for ID in `seq 1 ${TETRAD_COUNT};
do
	SAMPLE1=${PARENTS}_${ID}A_R1.fastq.gz # Passed sample prefix (ex: Sample-01)
	SAMPLE2=${PARENTS}_${ID}B_R1.fastq.gz # Passed sample prefix (ex: Sample-01)
	SAMPLE3=${PARENTS}_${ID}C_R1.fastq.gz # Passed sample prefix (ex: Sample-01)
	SAMPLE4=${PARENTS}_${ID}D_R1.fastq.gz # Passed sample prefix (ex: Sample-01)
	SAMPLE=${PARENTS}_${ID}
	DIR=${OUTPUTDIR}
	WORKDIR=${OUTPUTDIR} # Where files will be created
	SEQDIR=${TETRADDIR} # Location of Fastqs
	MODDIR=/usr/local/bin
	JAVA=/home2/cheil/programs/jre1.8.0_231/bin #GATK only works with Java v1.8
	SEQID=Yeasttetrad
	REF=${REFERENCE} # Reference genome
	#REF=/home2/cheil/reference_seq/sacCer3.fasta
	SCRIPTS=/home4/eschwar3/CO_NCO/1.scripts # Location of custom scripts
	VCFDIR=${WORKDIR}/vcf
	BCFTOOLS=/home4/eschwar3/bin/bin #Location of bcftools
	
	cd ${WORKDIR}

	mkdir ${VCFDIR}

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
		--variant ${VCFDIR}/${PARENT1}_GATK_HaplotypeCaller.vcf \
		--variant ${VCFDIR}/${PARENT2}_GATK_HaplotypeCaller.vcf \
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

# Make seg file

rm tmp.txt
touch tmp.txt

vcf=${VCFDIR}/${SAMPLE}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf
loc=${OUTPUTDIR}/PARENTAL_FILTERED_LOCS_FILE_GENOTYPE.txt

while read line 
do
	POS=$(echo $line | sed 's/ /	/g' | cut -f 1,2)
	P1=$(echo $line | sed 's/ /	/g' | cut -f 3)
	TETRAD=$(grep "$POS	" $vcf)
	REF=$(echo $TETRAD | sed 's/ /	/g' | cut -f 4)
	ALT=$(echo $TETRAD | sed 's/ /	/g' | cut -f 5)
	T1=$(echo $TETRAD | sed 's/ /	/g' | cut -f 10 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')
	T2=$(echo $TETRAD | sed 's/ /	/g' | cut -f 11 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')
	T3=$(echo $TETRAD | sed 's/ /	/g' | cut -f 12 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')
	T4=$(echo $TETRAD | sed 's/ /	/g' | cut -f 13 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')

	>&2 echo $POS

	if [ $T1 != 10 ] && [ $T1 != 01 ] && [ $T2 != 10 ] && [ $T2 != 01 ] && [ $T3 != 10 ] && [ $T3 != 01 ] && [ $T4 != 10 ] && [ $T4 != 01 ]
	then
		if [ $T1 == 00 ]
		then
			T1=$REF
		else
			T1=$ALT
		fi
		if [ $T1 == $P1 ]
		then
			T1=0
		else
			T1=1
		fi
		if [ $T2 == 00 ]
		then
			T2=$REF
		else
			T2=$ALT
		fi
		if [ $T2 == $P1 ]
		then
			T2=0
		else
			T2=1
		fi
		if [ $T3 == 00 ]
		then
			T3=$REF
		else
			T3=$ALT
		fi
		if [ $T3 == $P1 ]
		then
			T3=0
		else
			T3=1
		fi
		if [ $T4 == 00 ]
		then
			T4=$REF
		else
			T4=$ALT
		fi
		if [ $T4 == $P1 ]
		then
			T4=0
		else
			T4=1
		fi
		echo "$POS	.	$T1	$T2	$T3	$T4" >> tmp.txt
	fi
done < $loc


sed 's/Suva_//g' tmp.txt > ${OUTPUTDIR}/segfiles/${PARENTS}_${ID}.txt

echo '"${PARENTS}_${ID}",' >> tmp_tetradnames_${PARENTS}.txt

done

${TETRAD_NAMES}=$(cat tmp_tetradnames_${PARENTS}.txt)

sed -i "s/tetradname=['test1']/tetradname=['"${PARENTS}"']/g" ${CROSSOVERDIR}/crossOver.py

sed -i 's/files = [["test1_sample1"]]/files = [['${TETRAD_NAMES%,}']]/g' ${CROSSOVERDIR}/crossOver.py

CROSSOVERNAME=$(echo ${CROSSOVERDIR} | rev | cut -d '/' -f 1 | rev)

python2 ${OUTPUTDIR}/${CROSSOVERNAME}/crossOver.py

mv -R ${OUTPUTDIR}/${CROSSOVERNAME}/out/${PARENTS} ${OUTPUTDIR}