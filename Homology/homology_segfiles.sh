#!/bin/bash
#SBATCH -J homology_segfiles
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/homology_segfiles_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/homology_segfiles_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file


PARENTS=
VCFDIR=
OUTPUTDIR=

while [ $# -gt 0 ] ; do
	case $1 in
		--parents) PARENTS=$2 ;;
		--vcf-directory) VCFDIR=$2 ;;
		--output-directory) OUTPUTDIR=$2 ;;
	esac

	shift
done

# USAGE

# sbatch homology_segfiles.sh --parents 'SRR1119189XSRR1119180' --vcf-directory '/home4/eschwar3/CO_NCO/3.output/2.CrossOver/SRR1119189XSRR1119180/vcf/' --output-directory '/home4/eschwar3/CO_NCO/3.output/2.CrossOver/SRR1119189XSRR1119180/'

vcf=${VCFDIR}/${PARENTS}_GATK_GenotypeGVCFs_SNPS_Tagged_Filtered.vcf
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

	if [ $GEN1 == '01' ] && [ $GEN2 == '11' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi

	if [ $GEN1 == '01' ] && [ $GEN2 == '00' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '10' ] && [ $GEN2 == '11' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '10' ] && [ $GEN2 == '00' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '11' ] && [ $GEN2 == '01' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '11' ] && [ $GEN2 == '10' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '00' ] && [ $GEN2 == '01' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '00' ] && [ $GEN2 == '10' ]
        then
                echo "$POS      0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '10' ] && [ $GEN2 == '10' ]
        then
                echo "$POS	0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '10' ] && [ $GEN2 == '01' ]
        then
                echo "$POS	0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '01' ] && [ $GEN2 == '01' ]
        then
                echo "$POS	0.5" >> tmp.txt
        fi
	
	if [ $GEN1 == '01' ] && [ $GEN2 == '10' ]
        then
                echo "$POS	0.5" >> tmp.txt
        fi
        
	if [ $GEN1 == '11' ] && [ $GEN2 == '00' ]
        then
                echo "$POS	1" >> tmp.txt
        fi
        
	if [ $GEN1 == '00' ] && [ $GEN2 == '11' ]
        then
                echo "$POS	1" >> tmp.txt
        fi
done

sort -V -k1,1 -k2,2 tmp.txt | grep -v "," | grep -v "#" > ${OUTPUTDIR}/${PARENTS}_HOMOLOGY_SEGTFILE.txt
rm tmp.txt

sed -i 's/Suva_//g' ${OUTPUTDIR}/${PARENTS}_HOMOLOGY_SEGTFILE.txt
