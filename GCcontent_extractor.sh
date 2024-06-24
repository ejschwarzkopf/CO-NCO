#!/bin/bash
#SBATCH -J GCcontent
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/GCcontent_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/GCcontent_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file

cd ~/CO_NCO

file=~/CO_NCO/3.output/CO_NCO_Model_Table_noGC.txt

ref=~/CO_NCO/2.data/1.reference_genome/Sbay.ultrascaf.fasta

k=$(wc -l $file | cut -d ' ' -f 1)

echo "GC" > tmp.txt

for i in `seq 2 $k`; do
	chr=$(awk 'NR == '$i' { print $1 }' $file)
	start=$(awk 'NR == '$i' { print $2 }' $file)
	string=$(awk 'NR == '$((2 * $chr))' { print }' $ref)
	string=${string:$start:20000}
	Gcount=${string//[^G]}
	Gcount=${#Gcount}
	Ccount=${string//[^C]}
        Ccount=${#Ccount}
	GCcontent=$(awk 'BEGIN { print ('$Gcount'+'$Ccount')/20000 }')
	echo $GCcontent >> tmp.txt
done

rows_file=$(wc -l $file | cut -d ' ' -f 1)

rows_GC=$(wc -l tmp.txt | cut -d ' ' -f 1)

if [ ${rows_file} -eq ${rows_GC} ] 
then
	paste ~/CO_NCO/3.output/CO_NCO_Model_Table_noGC.txt tmp.txt > ~/CO_NCO/3.output/CO_NCO_Model_Table_final.txt
fi
