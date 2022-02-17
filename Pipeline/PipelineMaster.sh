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

# usage example:
sbatch PipelineMaster.sh \
--parent1 SRR1119189 \
--parent2 SRR1119180 \
--se2 \
--parent-directory /home4/eschwar3/CO_NCO/2.data/2.sequence_data/1.parents_fastq/s_uvarum_seq \
--tetrad-directory /home4/eschwar3/CO_NCO/2.data/2.sequence_data/2.offspring_seq/BigRunSequences/fastq_NVS114A_Heilv2/Heil \
--reference-genome /home4/eschwar3/CO_NCO/2.data/1.reference_genome/Sbay.ultrascaf.fasta \
--output-directory /home4/eschwar3/CO_NCO/3.output/2.CrossOver/SRR1119189XSRR1119180 \
--crossover-directory /home4/eschwar3/CrossOver_v6.3/

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

./Parent_UnambiguousSites.sh --parent1 ${PARENT1} --parent2 ${PARENT2} --reference-genome ${REFERENCE} --output-directory ${OUTPUTDIR}

cp -R ${CROSSOVERDIR} ${OUTPUTDIR}

CROSSOVERNAME=$(echo ${CROSSOVERDIR} | rev | cut -d '/' -f 2 | rev)

PARENTS=${PARENT1}X${PARENT2}

./Tetrad_CrossOver.sh --tetrad-directory ${TETRADDIR} --output-directory ${OUTPUTDIR} --parent1 ${PARENT1} --parent2 ${PARENT2} --crossover-directory ${CROSSOVERDIR} --reference-genome ${REFERENCE}

rm -R ${OUTPUTDIR}/${CROSSOVERNAME}













