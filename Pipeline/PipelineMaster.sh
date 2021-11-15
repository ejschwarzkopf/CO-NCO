#!/bin/bash
#SBATCH -J CO_NCO
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/CO_NCO_%J.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/CO_NCO_%J.out
##SBATCH --array=1-1 # Num of jobs = Num of rows in sample list file

# Input:
# parent SRR number 
# whether parents are paired-end sequenced
# parent fastqc directory
# tetrad fastqc directory
# reference genome file
# output folder

die() {
	printf '%s\n' "$1" >&2
	exit 1
}

PARENT1=
PARENT2=
SE1=0
SE2=0
PARENTDIR=
TETRADDIR=
REFERENCE=
OUTPUTDIR=
CROSSOVERDIR=

while [ $# -gt 0 ] ; do
	case $1 in
		--dir) DIR=$2 ;;
		--parent1) PARENT1=$2 ;;
		--parent2) PARENT2=$2 ;;
		--se1) SE1=1 ;;
		--se2) SE2=1 ;;
		--parent-directory) PARENTDIR=$2 ;;
		--tetrad-directory) TETRADDIR=$2 ;;
		--reference-genome) REFERENCE=$2 ;;
		--output-directory) OUTPUTDIR=$2 ;;
		--crossover-directory) CROSSOVERDIR=$2 ;;
	esac

	shift
done

if [[ $SE1 -eq 1 ]]
then
	./SE_Parent_Align.sh --parent ${PARENT1} --directory ${PARENTDIR} --reference-genome ${REFERENCE} --output-directory ${OUTPUTDIR}
else
	./PE_Parent_Align.sh --parent ${PARENT1} --directory ${PARENTDIR} --reference-genome ${REFERENCE} --output-directory ${OUTPUTDIR}
fi

if [[ $SE2 -eq 1 ]]
then
	./SE_Parent_Align.sh --parent ${PARENT2} --directory ${PARENTDIR} --reference-genome ${REFERENCE} --output-directory ${OUTPUTDIR}
else
	./PE_Parent_Align.sh --parent ${PARENT2} --directory ${PARENTDIR} --reference-genome ${REFERENCE} --output-directory ${OUTPUTDIR}
fi

cp -R ${CROSSOVERDIR} ${OUTPUTDIR}

CROSSOVERNAME=$(echo ${CROSSOVERDIR} | rev | cut -d '/' -f 1 | rev)

PARENTS=${PARENT1}X${PARENT2}

RUN TETRAD SCRIPT --tetrad-directory ${TETRADDIR} --output-directory ${OUTPUTDIR} --parent1 ${PARENT1} --parent2 ${PARENT2} --crossover-name ${CROSSOVERNAME}

rm -R ${OUTPUTDIR}/${CROSSOVERNAME}













