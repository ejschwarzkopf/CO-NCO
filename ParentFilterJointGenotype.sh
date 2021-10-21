#!/bin/bash
#SBATCH -J ParentFilterJointGenotype
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/ParentFilterJointGenotype_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/ParentFilterJointGenotype_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file

cd ~/CO_NCO

vcf=3.output/1.align_varcall/4.parents/vcf/SRR1119180XSRR1119189_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf


rm tmp.txt
touch tmp.txt


tail +65 $vcf | while read line
do
	#>&2 echo $line
	POS=$(echo $line | sed 's/ /	/g'| cut -f 1,2)
	GEN1=$(echo $line | sed 's/ /	/g'| cut -f 10 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')
	GEN2=$(echo $line | sed 's/ /	/g'| cut -f 11 | cut -d ':' -f 1 | sed 's*/**g' | sed 's/|//g')
	REFBASE=$(echo $line | sed 's/ /	/g' | cut -f 4)
	ALTBASE=$(echo $line | sed 's/ /	/g' | cut -f 5)

	>&2 echo "$GEN1 $GEN2"
	if [ $GEN1 == '11' ] && [ $GEN2 == '00' ]
	then
		echo "$POS	$ALTBASE	$REFBASE" >> tmp.txt
	fi

	if [ $GEN1 == '00' ] && [ $GEN2 == '11' ]
	then
		echo "$POS	$REFBASE	$ALTBASE" >> tmp.txt
	fi
done

sort -V -k1,1 -k2,2 tmp.txt | grep -v "," | grep -v "#" > PARENTAL_FILTERED_LOCS_FILE_GENOTYPE.txt
