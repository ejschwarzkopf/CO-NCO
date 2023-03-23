#!/bin/bash
#SBATCH -J Downsample
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_%a.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_%a.out

module load R/4.1.0


cd ~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling/subsample_$k/tables

Rscript -e 'CO_NCO_list<-list()
CO_NCO_table<-data.frame(CO_mean=c(), CO_se=c(), NCO_mean=c(), CO_mean=c())
std.err <- function(x) sd(x)/sqrt(length(x))
dens<-c(1,1.58,2,4,8,16,32,40)
for(k in 1:8){
filename<-paste0("subsample_", dens[k], "_table.txt")
CO_NCO_list[[k]]<-read.table(filename)
CO_NCO_table[k,]<-data.frame(CO_mean=mean(CO_NCO_list[[k]][,1]), CO_se=ste.err(CO_NCO_list[[k]][,1]), NCO_mean=mean(CO_NCO_list[[k]][,3]), CO_mean=std.err(CO_NCO_list[[k]][,3]))
}
rownames(CO_NCO_table)=dens
write.table(CO_NCO_table, "CO_NCO_table.txt, quote=F")'



