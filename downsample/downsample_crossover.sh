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

write.table(sample_list, "subsample_'$i'/replicate_'$j'/sample_list.txt", quote=F, col.names=F, row.names=F)
q()'

cd subsample_$i/replicate_$j/CrossOver_v6.3

for file in `ls segfiles`; do sed -i 's/ /	/g' segfiles/$file; done

python2 crossOver.py

cd out/SRR1119200XSRR1119199

sed -i 's/\[/"\[/g' CoList_SRR1119200XSRR1119199.txt
sed -i 's/\]/\]"/g' CoList_SRR1119200XSRR1119199.txt
sed -i 's/#//g' CoList_SRR1119200XSRR1119199.txt

sed -i 's/#//g' TractList_SRR1119200XSRR1119199.txt
sed -i 's/YJM SNP/YJM.SNP/g' TractList_SRR1119200XSRR1119199.txt
sed -i 's/_//g' TractList_SRR1119200XSRR1119199.txt
sed -i 's/\[/"\[/g' TractList_SRR1119200XSRR1119199.txt
sed -i 's/\]/\]"/g' TractList_SRR1119200XSRR1119199.txt

Rscript -e 'COs1<-read.table("CoList_SRR1119200XSRR1119199.txt",header=T)
NCOs1_raw<-read.table("TractList_SRR1119200XSRR1119199.txt",header=T)
NCOs1<-NCOs1_raw[which(NCOs1_raw$marker>=3),]
Suva_chrom_length<-read.table("~/CO_NCO/2.data/Suva_chrom_length.txt")
COs1_1<-COs1[which(COs1$chr==1),]
COs1_2<-COs1[which(COs1$chr==2),]
COs1_3<-COs1[which(COs1$chr==3),]
COs1_4<-COs1[which(COs1$chr==4),]
COs1_5<-COs1[which(COs1$chr==5),]
COs1_6<-COs1[which(COs1$chr==6),]
COs1_7<-COs1[which(COs1$chr==7),]
COs1_8<-COs1[which(COs1$chr==8),]
COs1_9<-COs1[which(COs1$chr==9),]
COs1_10<-COs1[which(COs1$chr==10),]
COs1_11<-COs1[which(COs1$chr==11),]
COs1_12<-COs1[which(COs1$chr==12),]
COs1_13<-COs1[which(COs1$chr==13),]
COs1_14<-COs1[which(COs1$chr==14),]
COs1_15<-COs1[which(COs1$chr==15),]
COs1_16<-COs1[which(COs1$chr==16),]
NCOs1_1<-NCOs1[which(NCOs1$chr==1),]
NCOs1_2<-NCOs1[which(NCOs1$chr==2),]
NCOs1_3<-NCOs1[which(NCOs1$chr==3),]
NCOs1_4<-NCOs1[which(NCOs1$chr==4),]
NCOs1_5<-NCOs1[which(NCOs1$chr==5),]
NCOs1_6<-NCOs1[which(NCOs1$chr==6),]
NCOs1_7<-NCOs1[which(NCOs1$chr==7),]
NCOs1_8<-NCOs1[which(NCOs1$chr==8),]
NCOs1_9<-NCOs1[which(NCOs1$chr==9),]
NCOs1_10<-NCOs1[which(NCOs1$chr==10),]
NCOs1_11<-NCOs1[which(NCOs1$chr==11),]
NCOs1_12<-NCOs1[which(NCOs1$chr==12),]
NCOs1_13<-NCOs1[which(NCOs1$chr==13),]
NCOs1_14<-NCOs1[which(NCOs1$chr==14),]
NCOs1_15<-NCOs1[which(NCOs1$chr==15),]
NCOs1_16<-NCOs1[which(NCOs1$chr==16),]
COs1_count_1=data.frame(chr=1,start=10000*((1:(Suva_chrom_length[1,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[1,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[1,1]/10000+1))-1)+5000.5)
COs1_count_2=data.frame(chr=2,start=10000*((1:(Suva_chrom_length[2,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[2,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[2,1]/10000+1))-1)+5000.5)
COs1_count_3=data.frame(chr=3,start=10000*((1:(Suva_chrom_length[3,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[3,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[3,1]/10000+1))-1)+5000.5)
COs1_count_4=data.frame(chr=4,start=10000*((1:(Suva_chrom_length[4,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[4,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[4,1]/10000+1))-1)+5000.5)
COs1_count_5=data.frame(chr=5,start=10000*((1:(Suva_chrom_length[5,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[5,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[5,1]/10000+1))-1)+5000.5)
COs1_count_6=data.frame(chr=6,start=10000*((1:(Suva_chrom_length[6,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[6,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[6,1]/10000+1))-1)+5000.5)
COs1_count_7=data.frame(chr=7,start=10000*((1:(Suva_chrom_length[7,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[7,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[7,1]/10000+1))-1)+5000.5)
COs1_count_8=data.frame(chr=8,start=10000*((1:(Suva_chrom_length[8,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[8,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[8,1]/10000+1))-1)+5000.5)
COs1_count_9=data.frame(chr=9,start=10000*((1:(Suva_chrom_length[9,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[9,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[9,1]/10000+1))-1)+5000.5)
COs1_count_10=data.frame(chr=10,start=10000*((1:(Suva_chrom_length[10,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[10,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[10,1]/10000+1))-1)+5000.5)
COs1_count_11=data.frame(chr=11,start=10000*((1:(Suva_chrom_length[11,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[11,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[11,1]/10000+1))-1)+5000.5)
COs1_count_12=data.frame(chr=12,start=10000*((1:(Suva_chrom_length[12,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[12,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[12,1]/10000+1))-1)+5000.5)
COs1_count_13=data.frame(chr=13,start=10000*((1:(Suva_chrom_length[13,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[13,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[13,1]/10000+1))-1)+5000.5)
COs1_count_14=data.frame(chr=14,start=10000*((1:(Suva_chrom_length[14,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[14,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[14,1]/10000+1))-1)+5000.5)
COs1_count_15=data.frame(chr=15,start=10000*((1:(Suva_chrom_length[15,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[15,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[15,1]/10000+1))-1)+5000.5)
COs1_count_16=data.frame(chr=16,start=10000*((1:(Suva_chrom_length[16,1]/10000+1))-1)+1,end=10000*((1:(Suva_chrom_length[16,1]/10000+1))),mid=10000*((1:(Suva_chrom_length[16,1]/10000+1))-1)+5000.5)
COs1_count_1$COcount=apply(COs1_count_1,1,function(x){length(which(COs1_1$pos.bp.>=x[2]&COs1_1$pos.bp.<=x[3]))})
COs1_count_2$COcount=apply(COs1_count_2,1,function(x){length(which(COs1_2$pos.bp.>=x[2]&COs1_2$pos.bp.<=x[3]))})
COs1_count_3$COcount=apply(COs1_count_3,1,function(x){length(which(COs1_3$pos.bp.>=x[2]&COs1_3$pos.bp.<=x[3]))})
COs1_count_4$COcount=apply(COs1_count_4,1,function(x){length(which(COs1_4$pos.bp.>=x[2]&COs1_4$pos.bp.<=x[3]))})
COs1_count_5$COcount=apply(COs1_count_5,1,function(x){length(which(COs1_5$pos.bp.>=x[2]&COs1_5$pos.bp.<=x[3]))})
COs1_count_6$COcount=apply(COs1_count_6,1,function(x){length(which(COs1_6$pos.bp.>=x[2]&COs1_6$pos.bp.<=x[3]))})
COs1_count_7$COcount=apply(COs1_count_7,1,function(x){length(which(COs1_7$pos.bp.>=x[2]&COs1_7$pos.bp.<=x[3]))})
COs1_count_8$COcount=apply(COs1_count_8,1,function(x){length(which(COs1_8$pos.bp.>=x[2]&COs1_8$pos.bp.<=x[3]))})
COs1_count_9$COcount=apply(COs1_count_9,1,function(x){length(which(COs1_9$pos.bp.>=x[2]&COs1_9$pos.bp.<=x[3]))})
COs1_count_10$COcount=apply(COs1_count_10,1,function(x){length(which(COs1_10$pos.bp.>=x[2]&COs1_10$pos.bp.<=x[3]))})
COs1_count_11$COcount=apply(COs1_count_11,1,function(x){length(which(COs1_11$pos.bp.>=x[2]&COs1_11$pos.bp.<=x[3]))})
COs1_count_12$COcount=apply(COs1_count_12,1,function(x){length(which(COs1_12$pos.bp.>=x[2]&COs1_12$pos.bp.<=x[3]))})
COs1_count_13$COcount=apply(COs1_count_13,1,function(x){length(which(COs1_13$pos.bp.>=x[2]&COs1_13$pos.bp.<=x[3]))})
COs1_count_14$COcount=apply(COs1_count_14,1,function(x){length(which(COs1_14$pos.bp.>=x[2]&COs1_14$pos.bp.<=x[3]))})
COs1_count_15$COcount=apply(COs1_count_15,1,function(x){length(which(COs1_15$pos.bp.>=x[2]&COs1_15$pos.bp.<=x[3]))})
COs1_count_16$COcount=apply(COs1_count_16,1,function(x){length(which(COs1_16$pos.bp.>=x[2]&COs1_16$pos.bp.<=x[3]))})
COs1_count_1$NCOcount=apply(COs1_count_1,1,function(x){length(which(NCOs1_1$pos.bp.>=x[2]&NCOs1_1$pos.bp.<=x[3]&NCOs1_1$tracttype!=1&NCOs1_1$tracttype!=10))})
COs1_count_2$NCOcount=apply(COs1_count_2,1,function(x){length(which(NCOs1_2$pos.bp.>=x[2]&NCOs1_2$pos.bp.<=x[3]&NCOs1_2$tracttype!=1&NCOs1_2$tracttype!=10))})
COs1_count_3$NCOcount=apply(COs1_count_3,1,function(x){length(which(NCOs1_3$pos.bp.>=x[2]&NCOs1_3$pos.bp.<=x[3]&NCOs1_3$tracttype!=1&NCOs1_3$tracttype!=10))})
COs1_count_4$NCOcount=apply(COs1_count_4,1,function(x){length(which(NCOs1_4$pos.bp.>=x[2]&NCOs1_4$pos.bp.<=x[3]&NCOs1_4$tracttype!=1&NCOs1_4$tracttype!=10))})
COs1_count_5$NCOcount=apply(COs1_count_5,1,function(x){length(which(NCOs1_5$pos.bp.>=x[2]&NCOs1_5$pos.bp.<=x[3]&NCOs1_5$tracttype!=1&NCOs1_5$tracttype!=10))})
COs1_count_6$NCOcount=apply(COs1_count_6,1,function(x){length(which(NCOs1_6$pos.bp.>=x[2]&NCOs1_6$pos.bp.<=x[3]&NCOs1_6$tracttype!=1&NCOs1_6$tracttype!=10))})
COs1_count_7$NCOcount=apply(COs1_count_7,1,function(x){length(which(NCOs1_7$pos.bp.>=x[2]&NCOs1_7$pos.bp.<=x[3]&NCOs1_7$tracttype!=1&NCOs1_7$tracttype!=10))})
COs1_count_8$NCOcount=apply(COs1_count_8,1,function(x){length(which(NCOs1_8$pos.bp.>=x[2]&NCOs1_8$pos.bp.<=x[3]&NCOs1_8$tracttype!=1&NCOs1_8$tracttype!=10))})
COs1_count_9$NCOcount=apply(COs1_count_9,1,function(x){length(which(NCOs1_9$pos.bp.>=x[2]&NCOs1_9$pos.bp.<=x[3]&NCOs1_9$tracttype!=1&NCOs1_9$tracttype!=10))})
COs1_count_10$NCOcount=apply(COs1_count_10,1,function(x){length(which(NCOs1_10$pos.bp.>=x[2]&NCOs1_10$pos.bp.<=x[3]&NCOs1_10$tracttype!=1&NCOs1_10$tracttype!=10))})
COs1_count_11$NCOcount=apply(COs1_count_11,1,function(x){length(which(NCOs1_11$pos.bp.>=x[2]&NCOs1_11$pos.bp.<=x[3]&NCOs1_11$tracttype!=1&NCOs1_11$tracttype!=10))})
COs1_count_12$NCOcount=apply(COs1_count_12,1,function(x){length(which(NCOs1_12$pos.bp.>=x[2]&NCOs1_12$pos.bp.<=x[3]&NCOs1_12$tracttype!=1&NCOs1_12$tracttype!=10))})
COs1_count_13$NCOcount=apply(COs1_count_13,1,function(x){length(which(NCOs1_13$pos.bp.>=x[2]&NCOs1_13$pos.bp.<=x[3]&NCOs1_13$tracttype!=1&NCOs1_13$tracttype!=10))})
COs1_count_14$NCOcount=apply(COs1_count_14,1,function(x){length(which(NCOs1_14$pos.bp.>=x[2]&NCOs1_14$pos.bp.<=x[3]&NCOs1_14$tracttype!=1&NCOs1_14$tracttype!=10))})
COs1_count_15$NCOcount=apply(COs1_count_15,1,function(x){length(which(NCOs1_15$pos.bp.>=x[2]&NCOs1_15$pos.bp.<=x[3]&NCOs1_15$tracttype!=1&NCOs1_15$tracttype!=10))})
COs1_count_16$NCOcount=apply(COs1_count_16,1,function(x){length(which(NCOs1_16$pos.bp.>=x[2]&NCOs1_16$pos.bp.<=x[3]&NCOs1_16$tracttype!=1&NCOs1_16$tracttype!=10))})
COs1_count_all=rbind(COs1_count_1,COs1_count_2,COs1_count_3,COs1_count_4,COs1_count_5,COs1_count_6,COs1_count_7,COs1_count_8,COs1_count_9,COs1_count_10,COs1_count_11,COs1_count_12,COs1_count_13,COs1_count_14,COs1_count_15,COs1_count_16)
windows=which(COs1_count_all$chr==4&COs1_count_all$start>=236000&COs1_count_all$end<=430000)
d=COs1_count_all[windows,]
std.err <- function(x) sd(x)/sqrt(length(x))
Count_output=data.frame(CO_mean=mean(d$COcount),CO_se=std.err(d$COcount),NCO_mean=c(d$NCOcount),NCO_se=c(d$NCOcount))
write.table(Count_output,"CO_NCO_table.txt",quote=F,row.names=F,col.names=F)'

mkdir ~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling/subsample_$i/tables
mv CO_NCO_table.txt ~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling/subsample_$i/tables/CO_NCO_table_$j.txt




