library(qvalue)
library(rrBLUP)
library(matrixcalc)
#install.packages("Matrix")
library(Matrix)
############################
align<-"hisat"
rep<-"BLUP" #
transf<-"rlog"
K=20

## kinship
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/")
prun<-"LD02" # "LD01", "LD02"

geno_LD02<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt",header=F,sep="\t")
#geno_LD02<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt",header=F,sep="\t")#
#geno_LD02<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt",header=F,sep="\t")
rownames(geno_LD02)<-geno_LD02[,1]
geno_LD02<-geno_LD02[,-1]
colnames(geno_LD02)<-rownames(geno_LD02)

### gene position
v4_gene_gff<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)
v4_gene_gff<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)
v4_gene_gff<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)

v4_gene_gff<-v4_gene_gff[which(v4_gene_gff$chr %in% c(1:10)),]
v4_gene_gff$chr<-as.numeric(as.character(v4_gene_gff$chr))
v4_gene_gff$start<-as.numeric(as.character(v4_gene_gff$start))
v4_gene_gff<-v4_gene_gff[,c(1,4,7)]

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/")
#taxa<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
taxa<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
#taxa<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
taxa<-as.character(taxa[,1])

kin<-geno_LD02[which(rownames(geno_LD02) %in% taxa),which(colnames(geno_LD02) %in% taxa)]
kin<-kin[order(rownames(kin)),order(colnames(kin))]
kin<-as.matrix(kin)

posdefmat <- function(mat) {
  if (is.positive.definite(round(mat, 18))) {
    g = mat
  }
  else {
    g <-nearPD(mat)$mat
    warning("The matrix was adjusted for the nearest positive definite matrix")
  }
  return(g)
}

kin.adj<-posdefmat(kin)
kin.test<-as.matrix(kin.adj)

### phenotype
pheno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
pheno.all<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
pheno.all<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")

pheno.all<-pheno.all[which(pheno.all$MLC_STANDARD %in% taxa),c(1,23,26)] # only test SD18_both for CE and FT
pheno.all<-pheno.all[order(pheno.all$MLC_STANDARD),]
pheno.all[is.na(pheno.all[,3]),3]<-mean(pheno.all[,3],na.rm=T)

## counts.log should be 310x20018 for the rlog data

## for comparison purpose
# K=10
# transf="log"
# counts.log<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel/",
#                              transf,"_",align,"_PEER",K,"_PEERres_rep",rep,".txt",sep=""),header=T,sep="\t",row.names=1)


counts.log<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel/",
                             transf,"_",align,"_PEER",K,"_PEERres_rep",rep,"_OLRM_final.txt",sep=""),header=T,sep="\t",row.names=1)
counts.log<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel/",
                             transf,"_",align,"_PEER",K,"_PEERres_rep",rep,"_OLRM_final.txt",sep=""),header=T,sep="\t",row.names=1)


counts.log<-t(counts.log) # row: marker; col: observation
counts.log<-counts.log[,order(colnames(counts.log))]
which(colnames(counts.log)!=pheno.all$MLC_STANDARD)
which(colnames(counts.log)!=colnames(kin))
counts.log<-cbind.data.frame(rownames(counts.log),counts.log)
colnames(counts.log)[1]<-"gene_name"
#counts.log<-merge(v4_gene_gff,counts.log,by="gene_name",all=F) # test 19857 genes that have position info
counts.log<-merge(v4_gene_gff,counts.log,by="gene_name",all.y=T) # test 19857 genes that have position info
nm_noPos<-nrow(counts.log[is.na(counts.log$chr),])
counts.log$chr[is.na(counts.log$chr)]<-11
counts.log$start[is.na(counts.log$start)]<-c(1:nm_noPos)


#######################################
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")
scores<-GWAS(pheno.all, counts.log, fixed=NULL, K=kin.test, n.PC=0,min.MAF=-Inf, n.core=10, P3D=FALSE, plot=TRUE)
#write.table(scores,"PEER_K10_P3D_kin_ce_ft_SD18_TWAS.txt",col.names=T,row.names=F,sep="\t",quote=F)
#write.table(scores,"PEER_K10_noP3D_kin_ce_ft_SD18_TWAS.txt",col.names=T,row.names=F,sep="\t",quote=F)
#write.table(scores,"rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(scores,"rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS_MtCp.txt",col.names=T,row.names=F,sep="\t",quote=F) ## Test all genes, include Mt,Cp and contigs

######################################
# Manhattan plots
######################################
library(qqman)
library(qvalue)

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")

#twas0<-read.table("rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS.txt",header=T,sep="\t")
twas0<-read.table("rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS_MtCp.txt",header=T,sep="\t")

unique(twas0$chr)

for(i in 4:ncol(twas0)){
  twas0[,i]<-10^((-1)*twas0[,i])
}

twas00<-twas0[which(twas0$chr %in% c(1:10)),] #19857 genes

TWAS_Cand<-vector()

#pdf(paste("ce_ft_TWAS_nonP3D_PEER20_man.pdf",sep=""),width=8,height=4)
pdf(paste("ce_ft_TWAS_nonP3D_PEER20_MtCp_man.pdf",sep=""),width=8,height=4)
for (i in 4:ncol(twas00)){
  twas_1<-twas00[,c(1:3,i)]
  twas_1<-twas_1[order(twas_1[,4]),]

  qobj <- qvalue(p = twas_1[,4])
  qvalues <- qobj$qvalues
  twas_1$qvalues<-qvalues

  #twas_1[which(twas_1$gene_name=="Zm00001d025180"),]
  twas.cand<-twas_1[1:ceiling(nrow(twas0)*0.01),]

  trait<-colnames(twas_1)[4]

  twas.cand<-cbind(twas.cand,rep(trait,nrow(twas.cand)))
  colnames(twas.cand)[c(4,6)]<-c("pvalue","trait")

  TWAS_Cand<-rbind(TWAS_Cand,twas.cand)


  twas<-twas_1[which(twas_1$chr %in% c(1:10)),]
  colnames(twas)[4]<-"pvalue"


  ## start of if else
  if (length(twas[which(twas$qvalues<0.05),4])==0){

    manhattan(twas, chr = "chr", bp = "start", p = "pvalue", snp = "gene_name",
              col = c("navy", "darkorange1"), chrlabs = NULL,
              suggestiveline = FALSE,genomewideline =FALSE,
              main=paste("TWAS kinship model, rlog, PEER K=20: ",trait,sep=""))

  }else {
    threshold<-(max(twas$pvalue[which(twas$qvalues<0.05)])+min(twas$pvalue[which(twas$qvalues>0.05)]))/2
    manhattan(twas, chr = "chr", bp = "start", p = "pvalue", snp = "gene_name",
              col = c("navy", "darkorange1"), chrlabs = NULL,
              suggestiveline = FALSE, genomewideline = -log10(threshold),
              main=paste("TWAS kinship model, rlog, PEER K=20: ",trait,sep=""))
  }
  ## end of if else
}
dev.off()
#write.table(TWAS_Cand,paste("TWAS_cand_001_ce_ft_noP3D_rlog_K=20.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(TWAS_Cand,paste("TWAS_cand_001_ce_ft_noP3D_rlog_K=20_MtCp.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

pdf(paste("ce_ft_TWAS_noP3D_PEER20_qq.pdf",sep=""),width=4,height=4)
for (i in 4:ncol(twas0)){
  twas<-twas0[,c(1:3,i)]
  wax<-colnames(twas)[4]
  colnames(twas)[4]<-"pvalue"

  qq(as.numeric(twas$pvalue),main=wax)
}
dev.off()
##############################
# add gene annotation
##############################
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D")

#file<-"TWAS_cand_001_ce_ft_noP3D_rlog_K=20"
file<-"TWAS_cand_001_ce_ft_noP3D_rlog_K=20_MtCp"

TWAS<-read.table(paste(file,".txt",sep=""),header=T,sep="\t")

#annotation<-read.delim("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
#annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/Users/zhenghao/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene_name"
## need import from csv
TWAS$gene_name<-as.character(TWAS$gene_name)

all_genes_v4ANNO<-merge(TWAS,annotation,by="gene_name",all.y=F,all.x=T)

setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection")

#write.table(all_genes_v4ANNO,paste(file,"_wAnno_final.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
write.table(all_genes_v4ANNO,paste(file,"_wAnno_MtCp.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
