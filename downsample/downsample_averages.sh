#!/bin/bash
#SBATCH -J Downsample
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_avg.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_avg.out

module load R/4.1.0


cd ~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling

for i in 0.05 0.75 1.95 2.95 4.45 7.87 14 1.35 3.225 1000; do cat subsample_$i/tables/CO_NCO_table_* > subsample_$i/tables/CO_NCO_${i}_table.txt; done
	
for i in 0.05 0.75 1.95 2.95 4.45 7.87 14 1.35 3.225 1000; do for j in `seq 1 20`; do echo $i; done > tmp.txt; paste tmp.txt subsample_$i/tables/CO_NCO_${i}_table.txt; done > CO_NCO_table_all.txt

Rscript -e 'CO_NCO_list<-list()
CO_NCO_table<-data.frame(CO_mean=c(), CO_se=c(), NCO_mean=c(), NCO_se=c())
std.err <- function(x) sd(x, na.rm=T)/sqrt(length(x))
dens<-c(0.05,0.75,1.95,2.95,4.45,7.87,14,1.35,3.225,1000)
for(k in 1:length(dens)){
filename<-paste0("subsample_", dens[k], "/tables/CO_NCO_", dens[k], "_table.txt")
CO_NCO_list[[k]]<-read.table(filename)
CO_NCO_table<-rbind(CO_NCO_table, data.frame(CO_mean=mean(CO_NCO_list[[k]][,1], na.rm=T), CO_se=std.err(CO_NCO_list[[k]][,1]), NCO_mean=mean(CO_NCO_list[[k]][,3], na.rm=T), NCO_se=std.err(CO_NCO_list[[k]][,3]), Tract_mean=mean(CO_NCO_list[[k]][,5], na.rm=T), Tract_se=std.err(CO_NCO_list[[k]][,5]), TractSD_mean=mean(CO_NCO_list[[k]][,7], na.rm=T), TractSD_se=std.err(CO_NCO_list[[k]][,7]), Count_2kb_mean=mean(CO_NCO_list[[k]][,8], na.rm=T),Count_2kb_se=std.err(CO_NCO_list[[k]][,8]),Count_5kb_mean=mean(CO_NCO_list[[k]][,9], na.rm=T),Count_5kb_se=std.err(CO_NCO_list[[k]][,9]),Prop_2kb_mean=mean(CO_NCO_list[[k]][,10], na.rm=T),Prop_2kb_se=std.err(CO_NCO_list[[k]][,10]),Prop_5kb_mean=mean(CO_NCO_list[[k]][,11], na.rm=T),Prop_5kb_se=std.err(CO_NCO_list[[k]][,11])))
}
rownames(CO_NCO_table)=dens
write.table(CO_NCO_table, "CO_NCO_table.txt", quote=F)'

sed -i '1s/^/Density /' CO_NCO_table.txt


