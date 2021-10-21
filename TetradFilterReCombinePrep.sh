#!/bin/bash
#SBATCH -J TetradFilterJointGenotype
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/TetradFilterJointGenotype_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/TetradFilterJointGenotype_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file

cd ~/CO_NCO

rm tmp.txt
touch tmp.txt

vcf=3.output/1.align_varcall/3.offspring/vcf/CSH345X347-All_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf
loc=PARENTAL_FILTERED_LOCS_FILE_GENOTYPE.txt

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

sed 's/Suva_//g' tmp.txt > ReCombineInput.txt

