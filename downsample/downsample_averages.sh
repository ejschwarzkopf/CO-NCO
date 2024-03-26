#!/bin/bash
#SBATCH -J Downsample
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_avg.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_avg.out

module load R/4.1.0


cd ~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling

for i in 1 1.58 2 4 8 11.9 16 32 40; do cat subsample_$i/tables/CO_NCO_table_* > subsample_$i/tables/CO_NCO_${i}_table.txt; done
	
for i in 1 1.58 2 4 8 11.9 16 32 40; do for j in `seq 1 20`; do echo $i; done > tmp.txt; paste tmp.txt subsample_$i/tables/CO_NCO_${i}_table.txt; done > CO_NCO_table_all.txt

Rscript -e 'CO_NCO_list<-list()
CO_NCO_table<-data.frame(CO_mean=c(), CO_se=c(), NCO_mean=c(), CO_mean=c())
std.err <- function(x) sd(x)/sqrt(length(x))
dens<-c(1,1.58,2,4,8,11.9,16,32,40)
for(k in 1:9){
filename<-paste0("subsample_", dens[k], "/tables/CO_NCO_", dens[k], "_table.txt")
CO_NCO_list[[k]]<-read.table(filename)
<<<<<<< HEAD
CO_NCO_table<-rbind(CO_NCO_table,data.frame(CO_mean=mean(CO_NCO_list[[k]][,1]), CO_se=std.err(CO_NCO_list[[k]][,1]), NCO_mean=mean(CO_NCO_list[[k]][,3]), NCO_se=std.err(CO_NCO_list[[k]][,3])))
=======
CO_NCO_table[k,]<-data.frame(CO_mean=mean(CO_NCO_list[[k]][,1]), CO_se=ste.err(CO_NCO_list[[k]][,1]), NCO_mean=mean(CO_NCO_list[[k]][,3]), CO_mean=std.err(CO_NCO_list[[k]][,3]))
>>>>>>> 23f734d4b228d29a2242253d2b66713d3296589d
}
rownames(CO_NCO_table)=dens
write.table(CO_NCO_table, "CO_NCO_table.txt", quote=F)'

sed -i '1s/^/Density /' CO_NCO_table.txt


