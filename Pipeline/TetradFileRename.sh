#!/bin/bash

while [ $# -gt 0 ] ; do
	case $1 in
		--dir) DIR=$2 ;;
	esac

	shift
done

echo $DIR
cd $DIR

cross1=SRR1119189XSRR1119180
cross2=SRR1119200XSRR1119199

# Process for the fist cross

rm tmp.txt

ls ${cross1}_[A-C]_*A_S* | cut -d '_' -f 2,3 | rev | cut -c2- | rev | uniq > tmp.txt

cross1_countA=$(wc -l tmp.txt | cut -d ' ' -f 1)

for ID in `seq 1 ${cross1_countA}`
do
	line=$(awk 'NR=='$ID tmp.txt)
	mv ${cross1}_${line}A_*_R1_001.fastq.gz ${cross1}_${ID}A_R1.fastq.gz
	mv ${cross1}_${line}B_*_R1_001.fastq.gz ${cross1}_${ID}B_R1.fastq.gz
	mv ${cross1}_${line}C_*_R1_001.fastq.gz ${cross1}_${ID}C_R1.fastq.gz
	mv ${cross1}_${line}D_*_R1_001.fastq.gz ${cross1}_${ID}D_R1.fastq.gz
	mv ${cross1}_${line}A_*_R2_001.fastq.gz ${cross1}_${ID}A_R2.fastq.gz
	mv ${cross1}_${line}B_*_R2_001.fastq.gz ${cross1}_${ID}B_R2.fastq.gz
	mv ${cross1}_${line}C_*_R2_001.fastq.gz ${cross1}_${ID}C_R2.fastq.gz
	mv ${cross1}_${line}D_*_R2_001.fastq.gz ${cross1}_${ID}D_R2.fastq.gz
done

ls ${cross1}_[A-C]_*E_S* | cut -d '_' -f 2,3 | rev | cut -c2- | rev | uniq > tmp.txt

cross1_countE=$(wc -l tmp.txt)

for ID in `seq 1 ${cross1_countE}`;
do
	IDE=$((${ID}+${cross1_countA}+1))
	line=$(awk 'NR=='$ID tmp.txt)
	mv ${cross1}_${line}E_*_R1_001.fastq.gz ${cross1}_${IDE}A_R1.fastq.gz
	mv ${cross1}_${line}F_*_R1_001.fastq.gz ${cross1}_${IDE}B_R1.fastq.gz
	mv ${cross1}_${line}G_*_R1_001.fastq.gz ${cross1}_${IDE}C_R1.fastq.gz
	mv ${cross1}_${line}H_*_R1_001.fastq.gz ${cross1}_${IDE}D_R1.fastq.gz
	mv ${cross1}_${line}E_*_R2_001.fastq.gz ${cross1}_${IDE}A_R2.fastq.gz
	mv ${cross1}_${line}F_*_R2_001.fastq.gz ${cross1}_${IDE}B_R2.fastq.gz
	mv ${cross1}_${line}G_*_R2_001.fastq.gz ${cross1}_${IDE}C_R2.fastq.gz
	mv ${cross1}_${line}H_*_R2_001.fastq.gz ${cross1}_${IDE}D_R2.fastq.gz
done

# Process for the second cross

rm tmp.txt

ls ${cross2}_[A-C]_*A_S* | cut -d '_' -f 2,3 | rev | cut -c2- | rev | uniq > tmp.txt

cross2_countA=$(wc -l tmp.txt)

for ID in `seq 1 ${cross2_countA}`
do
	line=$(awk 'NR=='$ID tmp.txt)
	mv ${cross2}_${line}A_*_R1_001.fastq.gz ${cross2}_${ID}A_R1.fastq.gz
	mv ${cross2}_${line}B_*_R1_001.fastq.gz ${cross2}_${ID}B_R1.fastq.gz
	mv ${cross2}_${line}C_*_R1_001.fastq.gz ${cross2}_${ID}C_R1.fastq.gz
	mv ${cross2}_${line}D_*_R1_001.fastq.gz ${cross2}_${ID}D_R1.fastq.gz
	mv ${cross2}_${line}A_*_R2_001.fastq.gz ${cross2}_${ID}A_R2.fastq.gz
	mv ${cross2}_${line}B_*_R2_001.fastq.gz ${cross2}_${ID}B_R2.fastq.gz
	mv ${cross2}_${line}C_*_R2_001.fastq.gz ${cross2}_${ID}C_R2.fastq.gz
	mv ${cross2}_${line}D_*_R2_001.fastq.gz ${cross2}_${ID}D_R2.fastq.gz
done

ls ${cross2}_[A-C]_*E_S* | cut -d '_' -f 2,3 | rev | cut -c2- | rev | uniq > tmp.txt

cross2_countE=$(wc -l tmp.txt)

for ID in `seq 1 ${cross2_countE}`;
do
	IDE=$((${ID}+${cross2_countA}+1))
	line=$(awk 'NR=='$ID tmp.txt)
	mv ${cross2}_${line}E_*_R1_001.fastq.gz ${cross2}_${IDE}A_R1.fastq.gz
	mv ${cross2}_${line}F_*_R1_001.fastq.gz ${cross2}_${IDE}B_R1.fastq.gz
	mv ${cross2}_${line}G_*_R1_001.fastq.gz ${cross2}_${IDE}C_R1.fastq.gz
	mv ${cross2}_${line}H_*_R1_001.fastq.gz ${cross2}_${IDE}D_R1.fastq.gz
	mv ${cross2}_${line}E_*_R2_001.fastq.gz ${cross2}_${IDE}A_R2.fastq.gz
	mv ${cross2}_${line}F_*_R2_001.fastq.gz ${cross2}_${IDE}B_R2.fastq.gz
	mv ${cross2}_${line}G_*_R2_001.fastq.gz ${cross2}_${IDE}C_R2.fastq.gz
	mv ${cross2}_${line}H_*_R2_001.fastq.gz ${cross2}_${IDE}D_R2.fastq.gz
done

rm tmp.txt