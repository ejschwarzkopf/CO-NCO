#!/bin/bash

PARENT1=
PARENT2=
PARENTDIR=
TETRADDIR=
REFERENCE=
OUTPUTDIR=

while :; do
    case $1 in
        --parent1 ?*)
            PARENT1=${1#* }
            ;;
        --parent2 ?*)
            PARENT2=${1#* }
            ;;
        --directory ?*)
            DIR=${1#* }
            ;;
        --reference-genome ?*)
            REFERENCE=${1#* }
            ;;
        --output-directory ?*)
            OUTPUTDIR=${1#* }
            ;;
    esac

    shift
done

VCFDIR=${OUTPUTDIR}/vcf

cd ${OUTPUTDIR}

mkdir ${VCFDIR}

PARENTS=$(echo ${PARENT1}X${PARENT2})

# All together now

(>&2 echo ***GATK - CombineGVCFs***)
gatk CombineGVCFs \
   -R ${REFERENCE} \
   --variant ${VCFDIR}/${PARENT1}_GATK_HaplotypeCaller.vcf \
   --variant ${VCFDIR}/${PARENT2}_GATK_HaplotypeCaller.vcf \
   -O ${VCFDIR}/${PARENTS}_GATK_HaplotypeCaller.vcf

(>&2 echo ***GATK - GenotypeGVCFs***)
gatk GenotypeGVCFs \
        -R ${REFERENCE} \
        -V ${VCFDIR}/${PARENTS}_GATK_HaplotypeCaller.vcf \
        -O ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs.vcf

(>&2 echo ***bcftools - filter - indel distance***)
${BCFTOOLS}/bcftools filter \
        --SnpGap 10 \
        -o ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_IndelFiltered.vcf \
        ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs.vcf

(>&2 echo ***GATK - SelectVariants***)
# Extract SNPs from call set.
gatk SelectVariants \
        -R ${REFERENCE} \
        -V ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_IndelFiltered.vcf \
        --select-type-to-include SNP \
        -O ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS.vcf

(>&2 echo ***GATK - VariantFiltration***)
# Apply filters to call sets.
gatk VariantFiltration \
        -R ${REFERENCE} \
        -V ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS.vcf \
        -filter "QUAL < 100" --filter-name "QUAL100" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged.vcf

(>&2 echo ***GATK - SelectVariants***)
# Select only variants that have passed the filters.
gatk SelectVariants \
        -R ${REFERENCE} \
        -V ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged.vcf \
        --selectExpressions 'vc.isNotFiltered()' \
        -O ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf

(>&2 echo ***GATK - VariantsToTable***)
# Variants to table
gatk VariantsToTable \
        -R $REFERENCE \
        -V ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf \
        -F CHROM -F POS -GF AD \
        -O ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.table

rm ${VCFDIR}/${PARENTS}_GATK_HaplotypeCaller.vcf
rm ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs.vcf
rm ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_IndelFiltered.vcf
rm ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs.vcf
rm ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS.vcf
rm ${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged.vcf

# Additional filtering to obtain unambiguous sites

vcf=${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf
rm tmp.txt
touch tmp.txt

tail +65 $vcf | while read line
do
    #>&2 echo $line
    POS=$(echo $line | sed 's/ /    /g'| cut -f 1,2)
    GEN1=$(echo $line | sed 's/ /   /g'| cut -f 10 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')
    GEN2=$(echo $line | sed 's/ /   /g'| cut -f 11 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')
    REFBASE=$(echo $line | sed 's/ /    /g' | cut -f 4)
    ALTBASE=$(echo $line | sed 's/ /    /g' | cut -f 5)

    >&2 echo "$GEN1 $GEN2"
    if [ $GEN1 == '11' ] && [ $GEN2 == '00' ]
    then
        echo "$POS  $ALTBASE    $REFBASE" >> tmp.txt
    fi

    if [ $GEN1 == '00' ] && [ $GEN2 == '11' ]
    then
        echo "$POS  $REFBASE    $ALTBASE" >> tmp.txt
    fi
done

sort -V -k1,1 -k2,2 tmp.txt | grep -v "," | grep -v "#" > ${OUTPUTDIR}/PARENTAL_FILTERED_LOCS_FILE_GENOTYPE.txt
rm tmp.txt