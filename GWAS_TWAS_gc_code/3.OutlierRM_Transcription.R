###############################################
# for PEER residuals of BLUPs
###############################################
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel")

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel")
#BLUP<-read.table("log_hisat_PEER25_PEERres_repBLUP.txt",header=T,row.names=1,sep="\t")
BLUP<-read.table("rlog_hisat_PEER20_PEERres_repBLUP.txt",header=T,row.names=1,sep="\t")
BLUP<-cbind.data.frame(rownames(BLUP),BLUP)
colnames(BLUP)[1]<-"MLC_STANDARD"

################################################
# outlier removal
################################################
source("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/script/outlier_removal_function_Meng.R")

nbvariables=1
nbtraits=ncol(BLUP)-1
random=c()
fixed=c()
BLUx="BLUE" # because the model used for outlier removal only contains grand mean, and there is no random factors in the model

clean.dataset <- BLUP[,c(1:nbvariables)]
alltraitnames <- vector()
#for (i in 1:10){
for (i in 1:nbtraits){
  curr.trait <- colnames(BLUP[ (nbvariables + i) ])
  alltraitnames <- c(alltraitnames, curr.trait)

  transfpheno=BLUP[,c(1:nbvariables,i+nbvariables)]
  initial.outliers(transfpheno, curr.trait, random, fixed, nbvariables, BLUx) -> out
  as.matrix(out) -> out
  colnames(out) <- curr.trait

  if(i == 1) {
    cleanedpheno <- cbind(clean.dataset, out)
  } else {
    cleanedpheno <- cbind(cleanedpheno, out)
  }

}
cleanedpheno<-as.data.frame(cleanedpheno)
cleanedpheno[,1]<-BLUP[,1]
colnames(cleanedpheno)[1]<-"MLC_STANDARD"
write.table(cleanedpheno,"rlog_hisat_PEER20_PEERres_repBLUP_OLRM.txt",col.names=T,row.names=F,sep="\t",quote=F)

########################################
#simple correlation between gene expression vs gc for plausible candidate genes
########################################
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel")
GE_clean<-read.table("rlog_hisat_PEER20_PEERres_repBLUP_OLRM.txt",header=T,row.names=1,sep="\t")
GE_clean2<-read.table("rlog_hisat_PEER20_PEERres_repBLUP_OLRM_final.txt",header=T,row.names=1,sep="\t")

candID<-c("Zm00001d005087","Zm00001d010426","Zm00001d015477","Zm00001d022364","Zm00001d032788","Zm00001d036765",
          "Zm00001d039411","Zm00001d043509","Zm00001d044162","Zm00001d049479","Zm00001d051599")
GE_clean_cand<-GE_clean[,which(colnames(GE_clean) %in% candID)]
GE_clean_cand<-cbind(rownames(GE_clean_cand),GE_clean_cand)
colnames(GE_clean_cand)[1]<-c("MLC_STANDARD")

GE_clean_cand2<-GE_clean2[,which(colnames(GE_clean2)=="Zm00001d049479")]
GE_clean_cand2<-cbind.data.frame(rownames(GE_clean2),GE_clean_cand2)
colnames(GE_clean_cand2)<-c("MLC_STANDARD","Zm00001d049479_2")

pheno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
#pheno.all<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
#pheno.all<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")

pheno.all<-pheno.all[,c(1,23)] # only test SD18_both for CE and FT

info<-merge(pheno.all,GE_clean_cand,by="MLC_STANDARD",all.y=T)
info<-merge(info,GE_clean_cand2,by="MLC_STANDARD",all.y=T) ## last col contain complete Zm00001d049479

## need unit here for gc
reg_plot<-function(vect1,vect2,Title){
  plot(vect1,vect2,main=Title,
       #ylab=expression(paste(italic("g")["c"],sep="")),
       ylab=expression(italic("g")["c"]~(g~h^{-1}~g^{-1})),
       xlab="PEER residuals of BLUPs for gene expression levels")
  fit=lm(vect2~vect1)
  abline(fit,col="red")
  r<-format(cor(vect1,vect2,use="complete.obs"), digits=2)
  legend("topleft",bty="n", legend=bquote(italic(r)~"="~.(r)))

  #legend("topleft",legend=expression(paste(italic(r)," = ",format(sqrt(summary(fit)$r.squared),digits=2),sep="")))
}

#setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")
pdf("Pearsons_corr_gc_vs_GE_candidates.pdf")
for ( i in 3:13){
  vect1<-info[,i]
  vect2<-info[,2] # gc
  Title=colnames(info)[i]

  if(i==13){
    plot(vect1,vect2,main=Title,
         #ylab=expression(paste(italic("g")["c"],sep="")),
         ylab=expression(italic("g")["c"]~(g~h^{-1}~g^{-1})),
         xlab="PEER residuals of BLUPs for gene expression levels")
    fit=lm(vect2~vect1)
    abline(fit,col="red")
    r<-format(cor(info[,14],vect2,use="complete.obs"), digits=2)
    legend("topleft",bty="n", legend=bquote(italic(r)~"="~.(r)))

  }else{
    reg_plot(vect1,vect2,Title)

  }

}
dev.off()


##########################################

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel")
BLUP<-read.table("rlog_hisat_PEER20_PEERres_repBLUP.txt",header=T,row.names=1,sep="\t")

cleanedpheno<-read.table("rlog_hisat_PEER20_PEERres_repBLUP_OLRM.txt",header=T,sep="\t")
nm_rm=0
for (i in 2:ncol(cleanedpheno)){
  rm<-length(cleanedpheno[is.na(cleanedpheno[,i]),i])
  nm_rm<-nm_rm+rm
}
nm_rm/((ncol(cleanedpheno)-1)*nrow(cleanedpheno))
# nm_rm=14560, 0.002346275 of the total data points

nm_NA<-c()
#cleanedpheno<-cleanedpheno[,-1]
for (i in 2:ncol(cleanedpheno)){
  nm<-length(cleanedpheno[is.na(cleanedpheno[,i]),i])
  nm_NA<-c(nm_NA,nm)
}
max(nm_NA)
range(nm_NA)

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel")
pdf("nm_NA_distribution.pdf")
hist(nm_NA,breaks=20)
hist(nm_NA[which(nm_NA>2)],breaks=20)
hist(nm_NA[which(nm_NA>4)],breaks=20)
dev.off()

length(nm_NA[which(nm_NA>31)])
flagged_gene<- colnames(cleanedpheno)[which(nm_NA>31)]
cleanedpheno<-cleanedpheno[,-which(colnames(cleanedpheno) %in% flagged_gene)] ## remove number of NA > 10% population (31)

write.table(cleanedpheno,"rlog_hisat_PEER20_PEERres_repBLUP_OLRM_final.txt",col.names=T,row.names=F,sep="\t",quote=F)

########### ENDS HERE ##############
####################################
