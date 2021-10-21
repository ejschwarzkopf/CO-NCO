#!/bin/bash
#SBATCH -J ParentFilter
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/ParentFilter_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/ParentFilter_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file

cd ~/CO_NCO

table1=3.output/1.align_varcall/1.parent1/vcf/SRR1119180_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.table
table2=3.output/1.align_varcall/2.parent2/vcf/SRR1119189_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.table
vcf1=3.output/1.align_varcall/1.parent1/vcf/SRR1119180_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf
vcf2=3.output/1.align_varcall/2.parent2/vcf/SRR1119189_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf


rm tmp.txt
touch tmp.txt


tail +2 $table1 | while read line
do
	#>&2 echo $line
	POS1=$(echo $line | sed 's/ /	/g'| cut -f 1,2)
	REF1=$(echo $line | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 1)
	ALT1=$(echo $line | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 2)

	ISIN=$(grep "$POS1	" $table2)
	#>&2 echo $ISIN

	if [ $REF1 -eq 0 ]
	then
		if [ ! -z "$ISIN" ]
		then
			REF2=$(echo $ISIN | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 1)
			ALT2=$(echo $ISIN | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 2)
			if [ $REF2 -eq 0 ]
			then
				ALTBASE1=$(grep "$POS1	" $vcf1 | cut -f 5)
				ALTBASE2=$(grep "$POS1	" $vcf2 | cut -f 5)
				if [ "$ALTBASE1" != "$ALTBASE2" ]
				then
					echo "$POS1	$ALTBASE1	$ALTBASE2" >> tmp.txt
				fi
			fi
		else
			ALTBASE1=$(grep "$POS1	" $vcf1 | cut -f 5)
			REFBASE=$(grep "$POS1	" $vcf1 | cut -f 4)
			echo "$POS1	$ALTBASE1	$REFBASE" >> tmp.txt
		fi
	fi
done

tail +2 $table2 | while read line
do
	#>&2 echo $line
	POS2=$(echo $line | sed 's/ /	/g'| cut -f 1,2)
	REF2=$(echo $line | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 1)
	ALT2=$(echo $line | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 2)

	ISIN=$(grep "$POS2	" $table1)
	#>&2 echo $ISIN

	if [ $REF2 -eq 0 ]
	then
		if [ ! -z "$ISIN" ]
		then
			REF1=$(echo $ISIN | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 1)
			ALT1=$(echo $ISIN | sed 's/ /	/g'| cut -f 3 | cut -d ',' -f 2)
			if [ $REF1 -eq 0 ]
			then
				ALTBASE1=$(grep "$POS2	" $vcf1 | cut -f 5)
				ALTBASE2=$(grep "$POS2	" $vcf2 | cut -f 5)
				if [ "$ALTBASE1" != "$ALTBASE2" ]
				then
					echo "$POS2	$ALTBASE1	$ALTBASE2" >> tmp.txt
				fi
			fi
		else
			ALTBASE2=$(grep "$POS2	" $vcf2 | cut -f 5)
			REFBASE=$(grep "$POS2	" $vcf2 | cut -f 4)
			echo "$POS2	$REFBASE	$ALTBASE2" >> tmp.txt
		fi
	fi
done

sort -k1,1V tmp.txt > PARENTAL_FILTERED_LOCS_FILE.txt
