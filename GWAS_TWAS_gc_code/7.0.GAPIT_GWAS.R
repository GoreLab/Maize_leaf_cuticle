### generate FT as covairate for AZ and AllEnv
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT")
pheno.all<-read.table("CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
Taxa<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/filter_prune_468K/MLC_GBSSNP310_468K_v4_chr1_prun09.012.indv",stringsAsFactors=F)
FT<-pheno.all[pheno.all$MLC_STANDARD %in% Taxa[,1],c(1,7,11)]
FT[is.na(FT$ft_AZ_untr),2]<-mean(FT$ft_AZ_untr,na.rm = T)
write.table(FT,"FT_as_covariate_1617.txt",col.names=T,row.names=F,sep="\t",quote=F)

## for Table S1
pheno<-pheno.all[which(pheno.all$MLC_STANDARD %in% Taxa[,1]),c(1,23,26)]
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT")
write.table(pheno,"gc_FT_310_TableS1.txt",col.names=T,row.names=F,sep="\t",quote=F)

###############################################################
###############################################################
# HapMap3 GWAS
###############################################################
################################################################
chr=2
trait<-c(1:2) #ce_SD18_both_tr
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT")
#pheno.all<-read.table("CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
pheno.all<-read.table("Transformed_BLUP_ce_SD18_both.txt",header=T,sep="\t")


## if need FT as covariate, use the following line
#FT<-read.table("FT_as_covariate.txt",header=T,sep="\t")

## FT as covariate for AZ and AllEnv
#FT<-read.table("FT_as_covariate_1617.txt",header=T,sep="\t") # for 310 lines
# FT<-pheno.all[,c(1,7,11)]
# FT[is.na(FT[,2]),2]<-mean(FT[,2],na.rm=T)
# FT[is.na(FT[,3]),3]<-mean(FT[,3],na.rm=T)
#####################################################

## no pruning, 450 imputation, 310 filter
geno.all<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",chr,"_filter.hmp.txt",sep=""),
                    header=F,sep="\t",comment.char="")


myG<-geno.all

#trait<-c(2,12:20)# 09242018 CE rate, just AZ, SD and all4
#trait<-c(2,12:19) # 12:15 for single env, 18:19 for 16 combined and 17 combined

myY<-pheno.all[,trait]
##### if log CE rate ###
#trait<-c(21:23)
#myY<-cbind(pheno.all[,1],log(pheno.all[,trait]))
#colnames(myY)[1]<-"Taxa"
#######################

## 450 imputation, 310 filtering
myKI<-read.table("centeredIBS_HMP3_AGPv4_LD02_450to310.txt")

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("multtest")
# install.packages("gplots")
# install.packages("LDheatmap")
# install.packages("genetics")
# install.packages("EMMREML")
# install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
#source("http://zzlab.net/GAPIT/gapit_functions.txt")
#source("/Users/Meng/Google Drive/MLC_AZ_2017/GAPIT/gapit_functions.R")
source("gapit_functions.R") # 2017 GAPIT
#source('http://www.zzlab.net/GAPIT/previous/gapit_functions20190714.txt') ## 2018 GAPIT, from Di
source("http://zzlab.net/GAPIT/emma.txt")

################################
# re-direct to folders of chr
###############################
mainDir<-"/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02"

subDir<-paste("chr",chr,sep="")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))
#########################################
#FT[is.na(FT[,2]),2]<-mean(FT[,2],na.rm=T)
#FT[is.na(FT[,3]),3]<-mean(FT[,3],na.rm=T)

#### without FT ######################################
nm_ind<-nrow(myKI)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  #CV=myCV, ## to be ditermined
  KI=myKI,
  #kinship.algorithm="Zhang",
  group.from=nm_ind, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=nm_ind,
  group.by=1,
  Major.allele.zero=T,
  Model.selection=FALSE
) #automatically include best number of K groups

####################################
# stacking result and Manhattan plot
####################################
library(qqman)
library(qvalue)

#Files<-c("SD18_1","SD18_2","SD18_both")
Files<-c("ce_AZ","ce_SD","ce_All","ce_SD18_both","ft_SD18_both")
#Files<-c("AZ","SD","All")
Files<-c("ce_SD18_both_tr")

for (f in Files){
  results<-matrix(nrow=0,ncol=10)
  mainDir="/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02"

  for (i in 1:10){
    subDir<-paste("chr",i,sep="")
    setwd(file.path(mainDir, subDir))

    file<-paste("GAPIT.MLM.",f,"_untr.GWAS.Results.csv",sep="")
    #file<-paste("GAPIT.CMLM.ce_",f,"_untr.GWAS.Results.csv",sep="")
    #file<-paste("GAPIT.MLM.",f,".GWAS.Results.csv",sep="")
    gemma<-read.csv(file,header=T)
    gemma$Chromosome<-i
    results<-rbind(results,gemma)
  }
  results<-results[order(results$Chromosome,results$Position),]

  setwd(mainDir)

  qobj <- qvalue(p = results$P.value)
  qvalues <- qobj$qvalues
  results$qvalues<-qvalues
  #results<-results[order(results$qvalues),]

   ## pruning
  GI.MP=results
  GI.MP=GI.MP[order(GI.MP$P.value),]
  topSNP<-GI.MP[1:5000,] # top 5000, ~0.1%
  topSNP.2<-GI.MP[1:ceiling(nrow(GI.MP)*0.1),] # top 10%
  topSNP.3<-GI.MP[1:ceiling(nrow(GI.MP)*0.0001),] # top 0.01%
  restSNP<-GI.MP[-(1:5000),]
  set.seed(89898)
  keep<-sample(1:nrow(restSNP),round(nrow(restSNP)*0.1,0),replace=F)
  restSNP<-restSNP[keep,]
  GI.MP.pruning<-rbind(topSNP,restSNP)

  # output for top SNPs
  #write.table(topSNP,paste("CE_",f,"_HapMap3_topSNPs.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables
  #write.table(topSNP.2,paste("CE_",f,"_HapMap3_topSNPs01.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for finsher's combined test

  write.table(topSNP,paste(f,"_HapMap3_topSNPs.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables
  write.table(topSNP.2,paste(f,"_HapMap3_topSNPs01.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for finsher's combined test
  write.table(topSNP.3,paste(f,"_HapMap3_topSNPs_forCand.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables

  #}
  ### plot
  if (length(results$P.value[which(results$qvalues<0.1)])==0){
    #pdf(paste("CE_",f,"_HapMap3_man_01.pdf",sep=""),width=8,height=4)  # "_1": 1% SNPs are highlighted; "_01": 5000 SNPs are highlighted
    pdf(paste(f,"_man_SNPforCand.pdf",sep=""),width=8,height=4)
    manhattan(GI.MP.pruning, chr = "Chromosome", bp = "Position", p = "P.value", snp = "SNP",
            col = c("navy", "darkorange1"), chrlabs = NULL,
            suggestiveline = FALSE,genomewideline =FALSE,
            highlight = as.character(topSNP.3$SNP),
            main=paste("GWAS ",f," HMP3",sep=""))
    dev.off()
  }else {
    threshold1<-(max(GI.MP$P.value[which(GI.MP$qvalues<0.1)])+min(GI.MP$P.value[which(GI.MP$qvalues>0.1)]))/2
    threshold2<-(max(GI.MP$P.value[which(GI.MP$qvalues<0.05)])+min(GI.MP$P.value[which(GI.MP$qvalues>0.05)]))/2

    #pdf(paste("CE_",f,"_HapMap3_man_01.pdf",sep=""),width=8,height=4) # "_1": 1% SNPs are highlighted; "_01": 5000 SNPs are highlighted
    pdf(paste(f,"_man_SNPforCand.pdf",sep=""),width=8,height=4)
    manhattan(GI.MP.pruning, chr = "Chromosome", bp = "Position", p = "P.value", snp = "SNP",
            col = c("navy", "darkorange1"), chrlabs = NULL,
            suggestiveline = -log10(threshold2), genomewideline = -log10(threshold1),
            highlight = as.character(topSNP.3$SNP),
            main=paste("GWAS ",f," HMP3",sep=""))
    dev.off()
  }
  #jpeg(paste("CE_",f,"_HapMap3_qq.jpg",sep=""))
  jpeg(paste(f,"_HapMap3_qq.jpg",sep=""))
  qq(results$P.value)
  dev.off()

}
