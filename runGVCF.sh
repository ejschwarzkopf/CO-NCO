#!/bin/bash

#SBATCH -J seqalign
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out

## Variant calling for S. uvarum population sequencing 
## To run on NCSU BRC with SLURM
## Chris Large and Caiti S. Heil
## Uses the recommended GATK pipeline for multiple samples
##Input is an aligned, processed, sorted and indexed bam file

SAMPLE=$1 # Passed sample prefix (ex: Sample-01)
DIR=/home2/cheil
WORKDIR=/home2/cheil/S_uva_Rec/vcf # Where files will be created
BAMDIR=/home2/cheil/S_uva_Rec/bam
MODDIR=/usr/local/bin
JAVA=/home2/cheil/programs/jre1.8.0_231/bin #GATK only works with Java v1.8
REF=/home2/cheil/reference_seq/Sbay.ultrascaf.fasta # Reference genome
SCRIPTS=/home2/cheil/scripts # Location of custom scripts
GVCFLIST=${WORKDIR}/tetrad_gvcf.list

cd ${WORKDIR}

# HaplotypeCaller in multiple sample mode



${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
	-R ${REF} \
	-T HaplotypeCaller \
	-I ${BAMDIR}/${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam \
	--emitRefConfidence GVCF \
	-o ${SAMPLE}_GATK_HaplotypeCaller.g.vcf




#${JAVA}/java -Xmx2g -jar ${MODDIR}/GenomeAnalysisTK.jar \
#        -R ${REF} \
#        -T HaplotypeCaller \
#        -I ${WORKDIR}/$line.1_comb_R1R2.RG.MD.realign.sort.bam \
#        --emitRefConfidence GVCF \
#        -o $line.1_GATK_HaplotypeCaller.g.vcf



