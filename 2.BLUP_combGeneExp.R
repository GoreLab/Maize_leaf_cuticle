# caluclate BLUPs on the rlog_clean data, 
# may need further filtering on n% zero values, but can subset them from BLUP results

########################
align<-"hisat"
transf<-"rlog" 
#########################
# reads_cutoff<-2 # 1 or 2 for 1M (2M) cutoff: reads per library
counts_1<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_cleanedQC_counts_v4_rep1_rlog_clean.txt",sep=""),row.names=1)
counts_2<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_cleanedQC_counts_v4_rep2_rlog_clean.txt",sep=""),row.names=1)

## not correct to use _1M. Need to check !!!!
# counts_1<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_QC_cmp1_cor74_rep1_1M.txt",sep=""),header=T,sep="\t",row.names=1)
# counts_2<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_QC_cmp1_cor74_rep2_1M.txt",sep=""),header=T,sep="\t",row.names=1)

# local
counts_1<-read.table(paste("/Users/Meng/Desktop/RNAseq_temp/Hisat_cleanedQC_counts_v4_rep1_rlog_clean.txt",sep=""),row.names=1)
counts_2<-read.table(paste("/Users/Meng/Desktop/RNAseq_temp/Hisat_cleanedQC_counts_v4_rep2_rlog_clean.txt",sep=""),row.names=1)

cm_gene<-intersect(rownames(counts_1),rownames(counts_2))
nm_gene<-length(cm_gene) # 21435 common genes, may need further filtering
#cm_line<-intersect(colnames(counts_1),colnames(counts_2))

#taxa<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt",header=F,stringsAsFactors=F)
taxa<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt",header=F,stringsAsFactors=F)

taxa<-taxa[,1]

counts_1_e<-counts_1[cm_gene,which(colnames(counts_1) %in% taxa)]
counts_1_c<-counts_1[cm_gene,grep("MO17",colnames(counts_1),fixed=T)]
cm_counts_1<-cbind(counts_1_e,counts_1_c)

counts_2_e<-counts_2[cm_gene,which(colnames(counts_2) %in% taxa)]
counts_2_c<-counts_2[cm_gene,grep("MO17",colnames(counts_2),fixed=T)]
cm_counts_2<-cbind(counts_2_e,counts_2_c)

#### for uploading RNA-seq data ####
setwd ("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/sequence_upload")
write.table(c(colnames(counts_1_e),colnames(counts_1_c)),"RNA-seq_samples_to_upload_rep1.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(c(colnames(counts_2_e),colnames(counts_2_c)),"RNA-seq_samples_to_upload_rep2.txt",col.names=F,row.names=F,sep="\t",quote=F)
length(c(colnames(counts_1_e),colnames(counts_1_c)))
#####################################

which(rownames(cm_counts_1)!=rownames(cm_counts_2))
#length(union(colnames(cm_counts_1),colnames(cm_counts_2)))
taxa[which(taxa %in% union(colnames(cm_counts_1),colnames(cm_counts_2)))]
#### run rep 1 and rep 2 sequencially, and stack the two results ##################
for (rep in 1:2){
  if (rep==1){
    counts=cm_counts_1
  } else {
    counts=cm_counts_2
  }
  ############################
  counts.log<-counts
  
  cand_counts <- as.data.frame(t(counts.log))
  gene_name<-colnames(cand_counts)
  
 #######################################
  
  ## experimental design ################
  #design<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/SD18_Design_Chk_Barcode_forBLUP.txt",sep=""),header=T,sep="\t")
  #design<-read.table(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/SD18_Design_Chk_Barcode_forBLUP.txt",sep=""),header=T,sep="\t")
  design<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/SD18_Design_Chk_Barcode_forBLUP.txt",sep=""),header=T,sep="\t")
  
  design$MLC_mf<-as.character(design$MLC_mf)
  design<-design[which(design$book.replication==rep),]
  design$MLC_mf[grep("^[0-9]",design$MLC_mf)]<-paste("X",design$MLC_mf[grep("^[0-9]",design$MLC_mf)],sep="")
  
  design$CHECK<-99
  design$CHECK[which(design$MLC_STANDARD=="MO17")]<-rep
  
  design$IS_EXPERIMENTAL<-1
  design$IS_EXPERIMENTAL[which(design$MLC_STANDARD=="MO17")]<-0
  
  design$COL1<-design$book.cols
  ## SD18 only
  design$COL1[which(design$book.cols>9)]<-19-design$book.cols[which(design$book.cols>9)]
  design<-design[,c(6,9:12,3,5)]
  colnames(design)<-c("MLC_STANDARD", "MLC_mf", "CHECK", "IS_EXPERIMENTAL", "COL1","BLOCK","rep1")
  
  # add plate information
  #plate_info<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/RNAseq_PlateInfo_v4_rep",rep,".txt",sep=""),header=T,sep="\t")
  #plate_info<-read.table(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/RNAseq_PlateInfo_v4_rep",rep,".txt",sep=""),header=T,sep="\t")
  plate_info<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/RNAseq_PlateInfo_v4_rep",rep,".txt",sep=""),header=T,sep="\t")
  
  ### for checking
  # taxa_t<-rownames(cand_counts)
  # taxa_pl<-as.character(plate_info[,1])
  # taxa_d<-as.character(design$MLC_mf)
  # taxa_t_pl<-intersect(taxa_t,taxa_pl)
  # taxa_t_d<-intersect(taxa_t,taxa_d)
  # taxa_pl_d<-intersect(taxa_pl,taxa_d)
  # taxa_t.1<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Individuals_after_filter_rep_1.txt",header=T,sep="\t")
  # taxa_t.2<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Individuals_after_filter_rep_2.txt",header=T,sep="\t")
  # taxa_t.both<-union(taxa_t.1[,1],taxa_t.2[,1])
  # taxa_final<-intersect(taxa_t.both,taxa_d)
  # taxa_pheno<-as.character(pheno.all$MLC_STANDARD[!is.na(pheno.all$ft_SD18_both_untr)])
  # taxa_final.2<-intersect(taxa_t.both,taxa_pheno)
  ## end of checking
  
  
  colnames(plate_info)[1]<-"MLC_mf"
  design<-merge(design,plate_info,by="MLC_mf",all=F)
  
  cand_counts<-cbind(rownames(cand_counts),cand_counts)
  colnames(cand_counts)[1]<-"MLC_mf"
  rownames(cand_counts)<-NULL
  #cand_counts_test<-merge(cand_counts[,1:5],design,by="MLC_mf",all=T)
  
  cand_counts<-merge(cand_counts,design,by="MLC_mf",all=F)
  #colnames(cand_counts)[ncol(cand_counts)-1]<-"CEadj"
  cand_counts$rep<-rep
  
  if(rep==1){
    cand_counts_1<-cand_counts
  }else{
    cand_counts_2<-cand_counts
  }
}
which(colnames(cand_counts_1)!=colnames(cand_counts_2))
############################################################
# Plot gene-based correlation between rep1 & rep2, before PEER
##############################################################

# cm_taxa<-intersect(as.character(cand_counts_1[,1]),as.character(cand_counts_2[,1]))
# cm_gene_taxa_1<-cand_counts_1[which(cand_counts_1[,1] %in% cm_taxa),]
# cm_gene_taxa_2<-cand_counts_2[which(cand_counts_2[,1] %in% cm_taxa),]
# which(as.character(cm_gene_taxa_1[,1])!=as.character(cm_gene_taxa_2[,1]))
# Correlations<-vector()
# for (i in 2:(1+length(cm_gene))){
#   correlation<-round(cor(cm_gene_taxa_1[,i],cm_gene_taxa_2[,i]),3)
#   Correlations<-c(Correlations,correlation)
# }
# setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/H2_est")
# pdf("gene-based_cor_rep1vsrep2_bfPEER.pdf",height=6,width=8)
# hist(Correlations,xlab="Correlation",main="gene-based correlation between rep1 & rep2, before PEER")
# dev.off()

###########################################################

############# after running rep 1 & 2, row bind cand_counts1 & 2

cand_counts_both<-rbind.data.frame(cand_counts_1,cand_counts_2) # 313 unique lines, 245 lines in both environments; 282 rep1 and 302 rep 2
length(unique(as.character(cand_counts_both$MLC_STANDARD))) # should be 311 including MO17, 310 experimental lines
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/comb_GeneExp")
write.table(cand_counts_both,"Hisat_cleanedQC_counts_v4_repBoth_rlog_clean_wDesign.txt",row.names=F,sep="\t",quote=F)

######## BLUPs calculation (only for combined exp) #############

#### Path of the license:
lic_path = "/workdir/ml2498/ASReml/License/"

# Load asreml:
setwd(lic_path)
library(asreml)
asreml.lic(license = "asreml.lic", install = TRUE)
#################################

cand_counts=cand_counts_both
for (i in 2:(nm_gene+1)){ # for subset H2
  #for (i in 2:(length(gene_name)+1)){ # for whole set BLUPs
  cand_counts[,i]<-as.numeric(as.character(cand_counts[,i]))
}
cand_counts$COL1<-as.factor(cand_counts$COL1)
cand_counts$BLOCK<-as.factor(cand_counts$BLOCK)
cand_counts$CHECK<-as.factor(cand_counts$CHECK)
cand_counts$Plate<-as.factor(cand_counts$Plate)
cand_counts$rep<-as.factor(cand_counts$rep)
cand_counts$MLC_STANDARD<-as.factor(as.character(cand_counts$MLC_STANDARD))
cand_counts$IS_EXPERIMENTAL<-as.numeric(cand_counts$IS_EXPERIMENTAL)

Res<-vector()
Gvar<-vector()
GXEvar<-vector()
plate<-vector()
block<-vector()
colm<-vector()
rep<-vector()

######################################
# If only want to have BLUPs
######################################
# only run this line when BLUPs are needed
count_Blup<-vector()
#Residual<-vector()

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/comb_GeneExp")
#setwd("/home/ml2498/Desktop/GeneExpression/comb_GeneExp")

#pdf("Hist_BLUPs_res_GeneExpression.pdf",height = 4,width = 5)
#for (i in 101:200){
for (i in 2:(nm_gene+1)){
  
  GeneID<-colnames(cand_counts)[i]
  
  # ## if run it at local computer
  # tryCatch({
  #   fit.asr <- eval(parse(text=paste("asreml(fixed = ",GeneID," ~ CHECK, random = ~ MLC_STANDARD:IS_EXPERIMENTAL+rep+rep/BLOCK+rep/COL1+rep/Plate+MLC_STANDARD:IS_EXPERIMENTAL:rep,na.action=na.method(x='omit'),data = cand_counts)",sep="")))
  # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # ##################################

    ### if run it on server
  fit.asr <- eval(parse(text=paste("asreml(fixed = ",GeneID," ~ CHECK, random = ~ MLC_STANDARD:IS_EXPERIMENTAL+rep+rep/BLOCK+rep/COL1+rep/Plate+MLC_STANDARD:IS_EXPERIMENTAL:rep,na.method.X='omit',data = cand_counts)",sep="")))
  
  #print(summary(fit.asr)$varcomp)
  var_res<-summary(fit.asr)$varcomp[7,1]
  Res<-c(Res,var_res)
  
  var_geno<-summary(fit.asr)$varcomp[5,1]
  Gvar<-c(Gvar,var_geno)
  
  var_gxe<-summary(fit.asr)$varcomp[6,1]
  GXEvar<-c(GXEvar,var_gxe)
  
  var_r<-summary(fit.asr)$varcomp[1,1]
  rep<-c(rep,var_r)
  
  var_p<-summary(fit.asr)$varcomp[2,1]
  plate<-c(plate,var_p)
  
  var_b<-summary(fit.asr)$varcomp[4,1]
  block<-c(block,var_b)
  
  var_c<-summary(fit.asr)$varcomp[3,1]
  colm<-c(colm,var_c)
  
  blup<-fit.asr$coefficients$random
  
  ## laptop (asreml4)
  #blup<-blup[-c(1:112,grep("MO17",names(blup))),]
  
  ## server (asreml3)
  blup<-blup[-c(1:112,grep("MO17",names(blup)))]
  
  intercept<-fit.asr$coefficients$fixed[4]
  check99<-fit.asr$coefficients$fixed[3]
  blup<-round(blup+intercept+check99,3)
  count_Blup<-cbind(count_Blup,blup)
  
  #residual<-fit.asr$residuals
  #residual<-residual[-grep("MO17",cand_counts$MLC_mf)]
  #residual<-round(residual,3)
  #Residual<-cbind(Residual,residual)
  
}

colnames(count_Blup)<-colnames(cand_counts)[2:(nm_gene+1)]
#colnames(Residual)<-colnames(Residual)[2:(nm_gene+1)]

#rownames(Residual)<-cand_counts$MLC_mf[-grep("MO17",cand_counts$MLC_mf)]

count_Blup1<-cbind(as.character(names(blup)),count_Blup)
rownames(count_Blup1)<-NULL
count_Blup1<-count_Blup1[grep("MLC_STANDARD_",count_Blup1[,1]),]
count_Blup1[,1]<-sub("MLC_STANDARD_","",count_Blup1[,1])
count_Blup1[,1]<-sub(":IS_EXPERIMENTAL","",count_Blup1[,1])

colnames(count_Blup1)[1]<-"MLC_STANDARD"
write.table(count_Blup1,"hisat_BLUPs_rlog_21Kgenes.txt",col.names=T,row.names=F,sep="\t",quote=F)

## store variance components
VC<-cbind(colnames(cand_counts)[2:(nm_gene+1)],rep, block,colm,plate,Gvar,GXEvar,Res)
colnames(VC)<-c("GeneID","var_rep","var_block","var_col","var_plate","var_geno","var_gxe","var_residual")
write.table(VC,"variance_components_proportional_trpt_abundance_rlog.txt",col.names=T,row.names=F,sep="\t",quote=F)

################################################################################
## Determine genes with too many zero's using "filter_out_genes_w_many_zeros.R"
## filter out those 1435 genes from BLUPs
################################################################################

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/comb_GeneExp")
rm<-read.table("remove_genes_with_many0s_fromBLUPs.txt",header=T,sep="\t")
count_Blup1<-read.table("hisat_BLUPs_rlog_21Kgenes.txt",header=T,sep="\t")

count_Blup1<-count_Blup1[,-which(colnames(count_Blup1) %in% as.character(rm$geneID))]

#############################################################
###### PEER
#############################################################
reads_cutoff=2

count_Blup<-count_Blup1

count_Blup<-count_Blup[-grep("rep",count_Blup$MLC_STANDARD),] # 312x17281, MO17 not included
rownames(count_Blup)<-count_Blup[,1]
write.table(count_Blup,"hisat_BLUPs_rlog_20019_genes.txt",row.names=F,sep="\t",quote=F)

count_Blup<-count_Blup[,-1]


K_test<-min(round(dim(count_Blup)[1]/4),100)
align<-"hisat"
rep<-"BLUP"
transf<-"rlog" 


library(peer)
set.seed(2010)
# build the model
model = PEER()
# run model
PEER_setNk(model,K_test)
# PEER_setNk(model,4)

PEER_setPhenoMean(model,as.matrix(count_Blup))
PEER_update(model) # Converged (var(residuals)) after 206 iterations, on server about 2.5-3 hours

# get precision
precision = PEER_getAlpha(model)
precision = cbind(paste("PEER",1:K_test,sep=""),precision)
colnames(precision)<-c("factor","precision")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel")
if (reads_cutoff==2){
  write.table(precision,paste(transf,"_",align,"_PEERfactor",K_test,"_precision_rep",rep,".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
} else if (reads_cutoff==1){
  write.table(precision,paste(transf,"_",align,"_PEERfactor",K_test,"_precision_rep",rep,"_1M.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}
  
###### plot 1/precision ############################
pdf(paste(transf,"_PEER_precision_rep",rep,".pdf",sep=""))
precision<-read.table(paste(transf,"_",align,"_PEERfactor",K_test,"_precision_rep",rep,".txt",sep=""),header=T,sep="\t")
precision$var<-1/as.numeric(as.character(precision$precision))
precision$factor<-as.numeric(precision$factor)
precision<-precision[order(precision$factor),]
prec_sub<-precision[2:40,]
prec_sub1<-precision[-1,]
plot(precision$factor,precision$var,main=paste("Precision of PEER factors: rep ",rep,sep=""),pch=1,xlab=NULL,ylab="1/Precision")
plot(prec_sub$factor,prec_sub$var,main=paste("Precision of PEER factors: rep ",rep,sep=""),
     ylim=c(0,0.05),
     pch=19,xlab=NULL,ylab="1/Precision",cex=0.5,lty=1)
plot(prec_sub1$factor,prec_sub1$var,main=paste("Precision of PEER factors: rep ",rep,sep=""),
     ylim=c(0,0.04),
     pch=19,xlab="Number of PEER factors",ylab="1/Precision",cex=0.5,lty=1)

abline(v=c(20),col="red")
dev.off()
#################################################
K<-20
align<-"hisat"
rep<-"BLUP"
transf<-"rlog"


library(peer)
set.seed(1987)
# build the model
model = PEER()
# run model
PEER_setNk(model,K)
# PEER_setNk(model,4)

PEER_setPhenoMean(model,as.matrix(count_Blup))
PEER_update(model) # Converged (var(residuals)) after 206 iterations, on server about 2.5-3 hours

PEERres = as.data.frame(PEER_getResiduals(model))
colnames(PEERres)<-colnames(count_Blup)
rownames(PEERres)<-rownames(count_Blup)
write.table(PEERres,paste(transf,"_",align,"_PEER",K,"_PEERres_rep",rep,".txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)

#############################
#############################





















