gff<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
gff<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
gff<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")

#GWAS_Peak<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_310kin_LD02/Peak_SNPs.txt",header=T,sep=",")
#GWAS_Peak<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_310kin_LD02/Peak_SNPs.txt",header=T,sep=",")

#GWAS_Peak<-read.table(paste("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_2/peak_snps_v3_ce_manual.csv",sep=""),header=T,sep=",") # modified based on "ce_SD18_Peak_SNPs.txt" results
GWAS_Peak<-read.table(paste("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_3_005/peak_snps_v3_ce_manual.csv",sep=""),header=T,sep=",") # modified based on "ce_SD18_Peak_SNPs.txt" results
GWAS_Peak<-read.table(paste("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_3_01/peak_snps_v3_ce_manual.txt",sep=""),header=T,sep="\t") # modified based on "ce_SD18_Peak_SNPs.txt" results
#GWAS_Peak<-GWAS_Peak[which(GWAS_Peak$PeakSNP=="Y"),]

## exclude helicopters
GWAS_Peak<-GWAS_Peak[which(GWAS_Peak$comments==""),]

GWAS_Peak$Position<-as.numeric(as.character(GWAS_Peak$Position))
GWAS_Peak<-GWAS_Peak[order(GWAS_Peak$Chromosome,GWAS_Peak$Position),]
GWAS_Peak$SNP<-as.character(GWAS_Peak$SNP)  # v4 results
GWAS_Peak$left250<-GWAS_Peak$Position-200000
GWAS_Peak$right250<-GWAS_Peak$Position+200000

#exp<-"ce_SD18_both" # "ce_SD18_both" or "ft_SD18_both"

####### gene position v4.42 #####
v4_gene_gtf<-read.table("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/fromPengfei/plot_orth_candidates/Zea_mays.AGPv4.42.gtf.short.txt",header=F,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)
v4_gene_gtf<-read.table("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/fromPengfei/plot_orth_candidates/Zea_mays.AGPv4.42.gtf.short.txt",header=F,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)

v4_gene_gtf$V9<-as.character(v4_gene_gtf$V9)
#test<-strsplit(v4_gene_gtf[1,9]," ")
gene_name<-strsplit(v4_gene_gtf[,9]," ") #the v4 gene names are in the last anotation column and need be pick out
gene_name<-sapply(gene_name, "[", 2)  #take the 2nd element of each component in a list
gene_name<- substr(gene_name,1,nchar(gene_name)-1)
v4_gene_gtf<-cbind(v4_gene_gtf[,-c(6:9)],gene_name)
colnames(v4_gene_gtf)<-c("chr","source","type","start","end","gene_name")

########################
gwas_hit<-GWAS_Peak
all_genes<-matrix(ncol=15,nrow=0)
for (i in 1:nrow(gwas_hit)){
  one_hit<-gwas_hit[i,]
  sub_genes<-v4_gene_gtf[which(v4_gene_gtf$chr==one_hit$Chromosome & v4_gene_gtf$start<one_hit$right250 & v4_gene_gtf$end>one_hit$left250),c(1,4:6)]
  nm_subgenes<-nrow(sub_genes)
  hit_info<-one_hit[rep(1, each=nm_subgenes),]
  sub_genes<-cbind(hit_info,sub_genes)
  all_genes<-rbind(all_genes,sub_genes)
}

## remove duplicated genes
length(unique(all_genes$gene_name))
keeplist2<-unique(all_genes$gene_name)
all_genes<-all_genes[match(keeplist2,all_genes$gene_name),]

setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02")
#write.table(all_genes,"ce_Top0001SNPs_Peak59_200K_cand_genes.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(all_genes,"ce_Top0005SNPs_Peak32_200K_cand_genes.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(all_genes,"ce_Top001SNPs_Peak22_200K_cand_genes.txt",col.names=T,row.names=F,quote=F,sep="\t")

################################################################################
#cand_genes<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_310kin_LD02/Peak8_LD05_cand_genes.txt",header=T,sep="\t")
cand_gwas<-read.table(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/ce_Top0005SNPs_Peak32_200K_cand_genes.txt",sep=""),header=T,sep="\t")
cand_gwas<-read.table(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/ce_Top001SNPs_Peak22_200K_cand_genes.txt",sep=""),header=T,sep="\t")
#cand_gwas<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/ce_Top0005SNPs_Peak32_200K_cand_genes.txt",sep=""),header=T,sep="\t")

#annotation<-read.delim("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene"
## need import from csv

cand_gwas$gene<-as.character(cand_gwas$gene)

all_genes_v4ANNO<-merge(cand_gwas,annotation,by="gene",all.y=F,all.x=T)

#setwd(mainDir)
#write.table(all_genes_v4ANNO,"GAPIT_cand_8Peaks_LD05_v4anno.txt",row.names=F,col.names=T,sep="\t",quote=F)

setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection")
write.table(all_genes_v4ANNO,"GWAS_Top0005SNPs_Peak32_200K_cand_genes_wAnno.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(all_genes_v4ANNO,"GWAS_Top001SNPs_Peak22_200K_cand_genes_wAnno.txt",row.names=F,col.names=T,sep="\t",quote=F)


####### Distance ############
#mainDir<-"/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_310kin_LD02" # adjust based on SNP set used for GWAS
mainDir<-"/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection" # adjust based on SNP set used for GWAS
setwd(mainDir)
#all_genes_v4ANNO<-read.delim("GAPIT_cand_8Peaks_LD05_v4anno.txt",header=T,sep="\t")
#all_genes_v4ANNO<-read.delim("GWAS_Top0005SNPs_Peak32_200K_cand_genes_wAnno.txt",header=T,sep="\t")
all_genes_v4ANNO<-read.delim("GWAS_Top001SNPs_Peak22_200K_cand_genes_wAnno.txt",header=T,sep="\t")

all_genes_v4ANNO$Dist<-NA
for (i in 1:nrow(all_genes_v4ANNO)){
  if (all_genes_v4ANNO$Position[i]>all_genes_v4ANNO$start[i]&all_genes_v4ANNO$Position[i]<all_genes_v4ANNO$end[i]){
    all_genes_v4ANNO$Dist[i]<-0
  }else {
    all_genes_v4ANNO$Dist[i]<-min(c(abs(all_genes_v4ANNO$Position[i]-all_genes_v4ANNO$start[i]),abs(all_genes_v4ANNO$Position[i]-all_genes_v4ANNO$end[i])))
  }

}

###### r2 ################
install.packages("lubridate")
library(lubridate)
library(data.table)

numericGeno<-function(geno){
  geno<-t(geno)
  colnames(geno)<-geno[1,];geno<-as.matrix(geno[-1,])

  alleles<-c("A","T","G","C")

  for (j in 1:ncol(geno)){

    freq<-table(geno[,j])
    freq<-freq[which(names(freq) %in% alleles)]

    freq<-freq[order(freq,decreasing=T)]
    major<-names(freq[1])
    minor<-names(freq[2])
    geno[,j]<-as.character(geno[,j])
    geno[which(geno[,j]==major),j]<-0
    geno[which(geno[,j]==minor),j]<-2
    geno[which(geno[,j]!=0&geno[,j]!=2),j]<-NA
  }
  return(geno)
}

both_end<-c("start","end")
all_genes_v4ANNO$LD_r2_max<-NA
all_genes_v4ANNO$LD_r2_median<-NA
all_genes_v4ANNO$LD_r2_mean<-NA
all_genes_v4ANNO$SNPs_in_gene<-NA

all_genes_v4ANNO<-all_genes_v4ANNO[order(all_genes_v4ANNO$Chromosome),]

snp_chr<-all_genes_v4ANNO$Chromosome[1]
# geno.all<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",snp_chr,"_prun09.hmp.txt",sep=""),
#                      header=T,sep="\t",comment.char="")
#geno.all<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",snp_chr,"_prun09.hmp.txt",sep=""),
#                     header=T,sep="\t",comment.char="")
geno.all<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",snp_chr,"_filter.hmp.txt",sep=""),
                     header=T,sep="\t",comment.char="")
geno.all<-fread(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",snp_chr,"_filter.hmp.txt",sep=""),
                     sep="\t",data.table=F)

for (i in 1:nrow(all_genes_v4ANNO)){
  SNP<-as.character(all_genes_v4ANNO$SNP[i])
  snp_chr<-all_genes_v4ANNO$Chromosome[i]
  snp_pos<-all_genes_v4ANNO$Position[i]

  ## need to change genotype file when change chromosomes (all_genes_v4ANNO has been ordered by chromosome)
  if (i>1){
    if(snp_chr!=all_genes_v4ANNO$Chromosome[i-1]){
      # geno.all<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",snp_chr,"_prun09.hmp.txt",sep=""),
      #                      header=T,sep="\t",comment.char="")
      #geno.all<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",snp_chr,"_prun09.hmp.txt",sep=""),
      #                     header=T,sep="\t",comment.char="")
      geno.all<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",snp_chr,"_filter.hmp.txt",sep=""),
                           header=T,sep="\t",comment.char="")


    }
  }

  gene_start<-all_genes_v4ANNO$start[i]
  gene_end<-all_genes_v4ANNO$end[i]

  Interval<-seq(from=min(c(gene_start,gene_end)),to=max(c(gene_start,gene_end)),by=1)
  snps_in_gene<-geno.all[which(geno.all$chrom==snp_chr&geno.all$pos %in% Interval),-(2:11)]

  if (all_genes_v4ANNO$Dist[i]==0){
    all_genes_v4ANNO$LD_r2_max[i]<-1
    all_genes_v4ANNO$LD_r2_median[i]<-1
    all_genes_v4ANNO$LD_r2_mean[i]<-1
  }else{
    if (nrow(snps_in_gene)>0){
      geno<-snps_in_gene
      all_genes_v4ANNO$SNPs_in_gene[i]<-"Y"
    } else {
      pick_end<-both_end[which.min(c(abs(snp_pos-gene_start),abs(snp_pos-gene_end)))]
      Interval2<-seq(from=min(c(all_genes_v4ANNO[i,pick_end],snp_pos)),to=max(c(all_genes_v4ANNO[i,pick_end],snp_pos)),by=1)
      geno.pre<-geno.all[which(geno.all$chrom==snp_chr,geno.all$pos %in% Interval2),-(2:11)]
      if (snp_pos>gene_start){
        geno<-geno.pre[1:10,]
      }else {
        geno<-tail(geno.pre,10)
      }
    }

    geno<-numericGeno(geno)
    ## geno.1 is the genotype of top SNP
    geno.1<-geno.all[c(which(geno.all$rs==SNP)),-(2:11)]
    geno.1<-numericGeno(geno.1)


    LD<-matrix(nrow=ncol(geno),ncol=2)
    for (j in 1:ncol(geno)){

      geno.2<-geno[,j]
      geno.11<-geno.1[which(!is.na(geno.1)&!is.na(geno.2))]
      geno.22<-geno.2[which(!is.na(geno.1)&!is.na(geno.2))]
      #temp<-as.data.frame(cbind(as.numeric(geno.11),as.numeric(geno.22)))
      #temp$comb<-paste(temp[,1],temp[,2],sep="_")

      len<-length(geno.11)
      p1<-length(which(geno.11==0))/len;p2<-1-p1
      q1<-length(which(geno.22==0))/len;q2<-1-q1
      p1q1<-length(which(geno.11==0&geno.22==0))/len # 0 is major allele, so there has to be individuals have 00 genotype
      r2<-round((p1q1-p1*q1)^2/(p1*q1*p2*q2),3)

      #LD[j,1]<-colnames(geno)[j]
      LD[j,2]<-r2
    }
    all_genes_v4ANNO$LD_r2_max[i]<-round(max(LD[,2]),3)
    all_genes_v4ANNO$LD_r2_median[i]<-round(median(LD[,2]),3)
    all_genes_v4ANNO$LD_r2_mean[i]<-round(mean(LD[,2]),3)
  }

}

setwd(mainDir)
#write.table(all_genes_v4ANNO,"GAPIT_cand_8Peaks_LD05_v4anno_LD_Dist.txt",row.names=F,col.names=T,sep="\t",quote=F)
#write.table(all_genes_v4ANNO,"GWAS_Top0005SNPs_Peak32_200K_cand_genes_wAnno_LD_Dist.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(all_genes_v4ANNO,"GWAS_Top001SNPs_Peak22_200K_cand_genes_wAnno_LD_Dist.txt",row.names=F,col.names=T,sep="\t",quote=F)
