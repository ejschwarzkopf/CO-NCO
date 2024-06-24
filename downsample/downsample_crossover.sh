#!/bin/bash
#SBATCH -J Downsample
#SBATCH --error=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_%a.err
#SBATCH --output=/home4/eschwar3/CO_NCO/1.scripts/1.logs/Downsample_%a.out
#SBATCH --array=1-200

# External variables

die() {
        printf '%s\n' "$1" >&2
        exit 1
}

help() {
  echo "Usage: $(basename $0) [-h] [-c CHROMOSOME] [-s START_POS] [-e END_POS] [-d DIRECTORY] [-a AUX_FILE_DIRECTORY]"
  echo "Options:"
  echo "  -h, --help     Display this help message"
  echo "  -c, --chrom    Set chromosome for the introgression"
  echo "  -s, --start    Set starting position for the introgression"
  echo "  -e, --end      Set end position for the introgression"
  echo "  -d, --dir      Set output directory"
  echo "  -a, --aux      Set auxiliary file for the script"
  echo "  -p, --pos      Set possible SNPs file"
  echo "  -x, --cross    Set cross name"
  exit 0
}

chrom=14
start=18500
end=586500
dir=~/CO_NCO/3.output/2.CrossOver/SRR1119200XSRR1119199_Unmasked/subsampling/
#dir=~/CO_NCO/1.scripts
aux=~/CO_NCO/2.aux/downsample_aux.txt
pos=possibleSNPs.txt
cross=SRR1119200XSRR1119199

while getopts ":c:s:e:d:a:" opt ; do
	case $opt in
		h) help ;;
		c) chrom=$OPTARG ;;
                s) start=$OPTARG ;;
                e) end=$OPTARG ;;
                d) dir=$OPTARG ;;
                a) aux=$OPTARG ;;
		p) pos=$OPTARG ;;
		x) cross=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        esac
done

echo "chrom = "$chrom
echo "start = "$start
echo "end = "$end
shift "$((OPTIND-1))"


# Load necessary modules

module load python/2.7.18
module load R/4.1.0

# Move to location of downsample output files

cd $dir

echo "dir = "$dir

# Sets the scenario for the downsampling. Equivalent to the row number of the "downsample_aux.txt" file

k=$SLURM_ARRAY_TASK_ID
echo "k = "$k

# From the current line number (k), we extract the SNP density (i) and replicate number (j)

i=$(awk 'NR=='$k' {print $1}' $aux)
j=$(awk 'NR=='$k' {print $2}' $aux)

if [ $k -eq 1 ]
then
	awk '$1=='$chrom' && $2>'$start' && $2<'$end' {print}' $pos > downsample_SNPs.txt
	
	for n in 0.05 0.75 1.95 2.95 4.45 7.87 14 1.35 3.225 1000; do
		for m in `seq 1 20`; do
			cp downsample_SNPs.txt subsample_$n/replicate_$m
		done
	done
else
	sleep 1.5m
fi

echo $i
echo $j

#seed=${chrom}${i}${j}${k}
seed=$(awk 'NR=='$k' {print $3}' $aux)

Rscript -e 'allSNPs<-read.table("subsample_'$i'/replicate_'$j'/downsample_SNPs.txt") # Read in the positions of SNPs to be downsampled (chr pos)
i='$i' # Set the SNP density
set.seed('$seed') # Sets a seed for the specific run
print("'$seed'")
length=(allSNPs[nrow(allSNPs),2] - allSNPs[1,2])
l=round((length * i)/1000)
#if(i==1){l=194}else if(i==1.58){l=307}else if(i==2){l=388}else if(i==4){l=776}else if(i==8){l=1552}else if(i==11.9){l=2309}else if(i==16){l=3104}else if(i==32){l=6208}else if(i==40){l=7760}else{NULL} # Establishes the number of SNPs to be sampled (# of SNPs in the region * (desired SNP density / observed SNP density))
#sample<-sort(floor(runif(n=l, min=1, max=nrow(allSNPs)))) # Samples without replacement from the number of SNPs
if(l>nrow(allSNPs)){l=nrow(allSNPs)}
print(l)
sample<-sort(sample(x=1:nrow(allSNPs), size = l, replace = FALSE)) # Similar to the commented out line above, but this time sample without replacement
sample_list=allSNPs[sample,] # Extracts the sampled SNPs into a table
# Load all of the seg files, extract the sampled snps and write the downsampled seg files
segfiles<-list()
segfiles_down<-list()
for(k in c(1:26, 28:48)){
filename1=paste0("'$dir'segfiles/'$cross'_", k, ".txt")
segfile<-read.table(filename1, header=F)
segfiles[[k]]<-segfile
segfiles_down[[k]]=rbind(segfiles[[k]][which(segfiles[[k]][,1]<'$chrom'),], segfiles[[k]][which(segfiles[[k]][,1]=='$chrom' & segfiles[[k]][,2]<'$start'),], segfiles[[k]][which(segfiles[[k]][,1] == '$chrom' & segfiles[[k]][,2] %in% sample_list[,2]),], segfiles[[k]][which(segfiles[[k]][,1]=='$chrom' & segfiles[[k]][,2]>'$end'),], segfiles[[k]][which(segfiles[[k]][,1]>'$chrom'),])
filename=paste0("subsample_'$i'/replicate_'$j'/CrossOver_v6.3/segfiles/'$cross'_", k, ".txt")
write.table(segfiles_down[[k]], filename, quote=F, col.names=F, row.names=F)
}
# This last part saves the sampled SNPs to a file
write.table(sample_list, "subsample_'$i'/replicate_'$j'/sample_list.txt", quote=F, col.names=F, row.names=F)
q()'

# Move to the CrossOver folder and run CrossOver on the subsampled segfiles

cd subsample_$i/replicate_$j/CrossOver_v6.3

for file in `seq 1 26` `seq 28 48`; do sed -i 's/ /	/g' segfiles/${cross}_$file.txt; done

python2 crossOver.py

# Prep CrossOver output for loading into R

cd out/$cross

sed -i 's/\[/"\[/g' CoList_$cross.txt
sed -i 's/\]/\]"/g' CoList_$cross.txt
sed -i 's/#//g' CoList_$cross.txt

sed -i 's/#//g' TractList_$cross.txt
sed -i 's/YJM SNP/YJM.SNP/g' TractList_$cross.txt
sed -i 's/_//g' TractList_$cross.txt
sed -i 's/\[/"\[/g' TractList_$cross.txt
sed -i 's/\]/\]"/g' TractList_$cross.txt

# Loads output into R and extract the average number of CO and NCO for subsampled and not subsampled segfiles

Rscript -e 'COs1<-read.table("CoList_'$cross'.txt",header=T)
NCOs1_raw<-read.table("TractList_'$cross'.txt",header=T)
NCOs1<-NCOs1_raw[which(NCOs1_raw$marker>=3&NCOs1_raw$tractminlen!=0&NCOs1_raw$tractlen<=10000),]
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
COs1_count_1=data.frame(chr=1,start=20000*((1:(Suva_chrom_length[1,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[1,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[1,1]/20000+1))-1)+10000.5)
COs1_count_2=data.frame(chr=2,start=20000*((1:(Suva_chrom_length[2,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[2,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[2,1]/20000+1))-1)+10000.5)
COs1_count_3=data.frame(chr=3,start=20000*((1:(Suva_chrom_length[3,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[3,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[3,1]/20000+1))-1)+10000.5)
COs1_count_4=data.frame(chr=4,start=20000*((1:(Suva_chrom_length[4,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[4,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[4,1]/20000+1))-1)+10000.5)
COs1_count_5=data.frame(chr=5,start=20000*((1:(Suva_chrom_length[5,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[5,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[5,1]/20000+1))-1)+10000.5)
COs1_count_6=data.frame(chr=6,start=20000*((1:(Suva_chrom_length[6,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[6,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[6,1]/20000+1))-1)+10000.5)
COs1_count_7=data.frame(chr=7,start=20000*((1:(Suva_chrom_length[7,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[7,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[7,1]/20000+1))-1)+10000.5)
COs1_count_8=data.frame(chr=8,start=20000*((1:(Suva_chrom_length[8,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[8,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[8,1]/20000+1))-1)+10000.5)
COs1_count_9=data.frame(chr=9,start=20000*((1:(Suva_chrom_length[9,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[9,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[9,1]/20000+1))-1)+10000.5)
COs1_count_10=data.frame(chr=10,start=20000*((1:(Suva_chrom_length[10,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[10,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[10,1]/20000+1))-1)+10000.5)
COs1_count_11=data.frame(chr=11,start=20000*((1:(Suva_chrom_length[11,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[11,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[11,1]/20000+1))-1)+10000.5)
COs1_count_12=data.frame(chr=12,start=20000*((1:(Suva_chrom_length[12,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[12,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[12,1]/20000+1))-1)+10000.5)
COs1_count_13=data.frame(chr=13,start=20000*((1:(Suva_chrom_length[13,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[13,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[13,1]/20000+1))-1)+10000.5)
COs1_count_14=data.frame(chr=14,start=20000*((1:(Suva_chrom_length[14,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[14,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[14,1]/20000+1))-1)+10000.5)
COs1_count_15=data.frame(chr=15,start=20000*((1:(Suva_chrom_length[15,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[15,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[15,1]/20000+1))-1)+10000.5)
COs1_count_16=data.frame(chr=16,start=20000*((1:(Suva_chrom_length[16,1]/20000+1))-1)+1,end=20000*((1:(Suva_chrom_length[16,1]/20000+1))),mid=20000*((1:(Suva_chrom_length[16,1]/20000+1))-1)+10000.5)
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
write.table(COs1_count_all, "COs_count_all.txt", quote=F,row.names=F,col.names=T)'

Rscript	-e 'COs1_count_all<-read.table("COs_count_all.txt", header=T)
NCOs1_raw<-read.table("TractList_'$cross'.txt",header=T)
NCOs1<-NCOs1_raw[which(NCOs1_raw$marker>=3&NCOs1_raw$tractminlen!=0&NCOs1_raw$tractlen<=10000),]
d=COs1_count_all[which(COs1_count_all$chr=='$chrom'&COs1_count_all$end>='$start'&COs1_count_all$start<='$end'),]
std.err<-function(x) sd(x)/sqrt(length(x))
Count_output=data.frame(CO_mean=mean(d$COcount),CO_se=std.err(d$COcount),NCO_mean=mean(d$NCOcount),NCO_se=std.err(d$NCOcount))
t_d=NCOs1[which(NCOs1$chr=='$chrom'&NCOs1$pos.bp>='$start'&NCOs1$pos.bp<='$end'),]
Tract_output=data.frame(Tract_mean=mean(t_d$tractlen),Tract_se=std.err(t_d$tractlen),Tract_sd=sd(t_d$tractlen),Count_2kb=length(which(t_d$tractlen>=2000)),Count_5kb=length(which(t_d$tractlen>=5000)),Prop_2kb=length(which(t_d$tractlen>=2000))/nrow(t_d),Prop_5kb=length(which(t_d$tractlen>=5000))/nrow(t_d))
Full_output=cbind(Count_output,Tract_output)' -e 'write.table(Full_output,"CO_NCO_table.txt",quote=F,row.names=F,col.names=F)
print(Full_output)'

#cut -d ' ' -f 2-8 CO_NCO_table.txt | tail -n 1 > tmp.txt
#mv tmp.txt CO_NCO_table.txt
mkdir ${dir}subsample_$i/tables
mv CO_NCO_table.txt ${dir}subsample_$i/tables/CO_NCO_table_$j.txt



