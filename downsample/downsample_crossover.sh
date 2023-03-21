#!/bin/bash
#SBATCH -J Downsample
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_%a.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_%a.out
#SBATCH --array=1-180

module load python/2.7.18
module load R/4.1.0

cd ~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling/

k=$SLURM_ARRAY_TASK_ID

i=$(awk 'NR=='$k' {print $1}' downsample_aux.txt)
j=$(awk 'NR=='$k' {print $2}' downsample_aux.txt)

Rscript -e 'allSNPs<-read.table("subsample_'$i'/replicate_'$j'/downsample_SNPs.txt")
i='$i'
if(i==1){l=194}else if(i==1.58){l=307}else if(i==2){l=38}else if(i==4){l=776}else if(i==8){l=1552}else if(i==11.9){l=2309}else if(i==16){l=3104}else if(i==32){l=6208}else if(i==40){l=7760}else{NULL}
sample<-sort(runif(n=l, min=1, max=10609))
sample_list=allSNPs[sample,]
# This is where I need to load all of the seg files, extract the sampled snps and write the downsampled seg files
segfiles<-list()
for(k in c(1:26, 28:48)){
filename1=paste0("~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling/segfiles/SRR1119200XSRR1119199_", k, ".txt")
segfile<-read.table(filename1, header=F)
segfiles[[k]]<-segfile
segfiles[[k]]=rbind(segfiles[[k]][which(segfiles[[k]][,1]<14),], segfiles[[k]][which(segfiles[[k]][,1] == 14 & segfiles[[k]][,2] %in% sample_list[,2]),], segfiles[[k]][which(segfiles[[k]][,1]>14),])
filename=paste0("subsample_'$i'/replicate_'$j'/CrossOver_v6.3/segfiles/SRR1119200XSRR1119199_", k, ".txt")
write.table(segfiles[[k]], filename, quote=F, col.names=F, row.names=F)
}

write.table(sample_list, "subsample_'$i'/replicate_'$j'/sample_list.txt", quote=F, col.names=F, row.names=F)'

cd subsample_$i/replicate_$j/CrossOver_v6.3

for file in `ls segfiles`; do sed -i 's/ /	/g' segfiles/$file; done

python2 crossOver.py

cd out/SRR119200XSRR1119199

sed -i 's/\[/"\[/g' CoList_SRR119200XSRR1119199.txt
sed -i 's/\]/\]"/g' CoList_SRR119200XSRR1119199.txt
sed -i 's/#//g' CoList_SRR119200XSRR1119199.txt






# Do whatever we're going to do with the output


