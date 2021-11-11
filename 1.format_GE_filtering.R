library(edgeR)

#align="star" # hisat or star
align="hisat"
##################################
######################################################
######################################################
# raw counts
######################################################
###################################################
raw_counts<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/meng-all-plates-raw-counts.txt",header=T,sep="\t",stringsAsFactors=F)
#raw_counts<-read.table("~/Desktop/GeneExpression/v4_counts/meng-all-plates-raw-counts.txt",header=T,sep="\t",,stringsAsFactors=F)

raw_counts<-raw_counts[,-c(grep("BLANK",colnames(raw_counts)),grep("Check",colnames(raw_counts)))]# remove BLANK and Check
#seperate raw_counts into rep1 and rep2
raw1<-raw_counts[,1]
raw2<-raw_counts[,1]
nm1<-"GeneID"
nm2<-"GeneID"
for (i in 2:ncol(raw_counts)){
  plot<-as.numeric(as.character(substr(colnames(raw_counts)[i],6,8)))
  if (plot<343){
    raw1<-cbind(raw1,raw_counts[,i])
    nm1<-c(nm1,colnames(raw_counts)[i])
  }else {
    raw2<-cbind(raw2,raw_counts[,i])
    nm2<-c(nm2,colnames(raw_counts)[i])
  }
}
colnames(raw1)<-nm1
colnames(raw2)<-nm2
##### remove the lines with alignment rate < 0.6
#setwd(paste("~/Desktop/GeneExpression/AlignmentRate_v4/",sep=""))
setwd(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/AlignmentRate_v4/",sep=""))

# rep 1
align_rate1<-read.csv("meng_P1_alignment_rate.csv",stringsAsFactors=F)
align_rate1$X[which(align_rate1$X=="P1-18SD156")]<-"P1-18SD156du"
align_rate1$X[which(align_rate1$X=="P1-18SD190")]<-"P1-18SD190du"

align_rate2<-read.csv("meng_P2_alignment_rate.csv",stringsAsFactors=F)
align_rate2$X[which(align_rate2$X=="P2-18SD223")]<-"P2-18SD223du"

align_rate3<-read.csv("meng_P3_alignment_rate.csv",stringsAsFactors=F)
align_rate4<-read.csv("meng_P4_alignment_rate.csv",stringsAsFactors=F)

align_rate<-rbind(align_rate1,align_rate2,align_rate3,align_rate4)
align_rate$Uniquely.alinged=as.numeric(sub("%", "", align_rate$Uniquely.alinged))
align_rate$Mutiple.aligned=as.numeric(sub("%", "", align_rate$Mutiple.aligned))
align_rate$total<-align_rate$Uniquely.alinged+align_rate$Mutiple.aligned
align_rate$X<-substr(align_rate$X,4,nchar(align_rate$X))
colnames(align_rate)[1]<-"Library"

QC1<-read.table("QC_rep1.txt",header=T,sep="\t")
QC_align_1<-merge(QC1,align_rate,by="Library",all=T)
write.table(QC_align_1,"QC_align_rep1.txt",col.names=T,row.names=F,sep="\t",quote=F)

# rep 2
align_rate5<-read.csv("meng_P5_alignment_rate.csv",stringsAsFactors=F)
align_rate6<-read.csv("meng_P6_alignment_rate.csv",stringsAsFactors=F)
align_rate7<-read.csv("meng_P7_alignment_rate.csv",stringsAsFactors=F)
align_rate8<-read.csv("meng_P8_alignment_rate.csv",stringsAsFactors=F)

align_rate<-rbind(align_rate5,align_rate6,align_rate7,align_rate8)
align_rate$Uniquely.alinged=as.numeric(sub("%", "", align_rate$Uniquely.alinged))
align_rate$Mutiple.aligned=as.numeric(sub("%", "", align_rate$Mutiple.aligned))
align_rate$total<-align_rate$Uniquely.alinged+align_rate$Mutiple.aligned
align_rate$X<-substr(align_rate$X,9,nchar(align_rate$X))
colnames(align_rate)[1]<-"Library"

QC2<-read.table("QC_rep2.txt",header=T,sep="\t")
QC_align_2<-merge(QC2,align_rate,by="Library",all=T)
write.table(QC_align_2,"QC_align_rep2.txt",col.names=T,row.names=F,sep="\t",quote=F)
########################################################

reads_cutoff<-2 # 1 or 2 for 1M (2M) cutoff: reads per library

for (rep in 1:2){
  # server
  setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/AlignmentRate_v4/")
  QC_align<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/AlignmentRate_v4/QC_align_rep",rep,".txt",sep=""),header=T,sep="\t")

  # local
  # setwd("~/Desktop/GeneExpression/AlignmentRate_v4/")
  # QC_align<-read.table(paste("~/Desktop/GeneExpression/AlignmentRate_v4/QC_align_rep",rep,".txt",sep=""),header=T,sep="\t")

  QC_align<-QC_align[-c(grep("BLANK",QC_align$Library),grep("Check",QC_align$Library)),]

  # pdf(paste("QC_plot_v4_rep",rep,".pdf",sep=""))
  # hist(QC_align$PropPolyA12)
  # hist(QC_align$PropReadThrough)
  # length(QC_align$Library[which(QC_align$TotalReads<2e6)])
  # plot(QC_align$TotalReads,QC_align$total,xlim=c(0,15000000))
  # abline(v=c(1e6,2e6),h=60)
  # hist(QC_align$TotalReads,xlim=c(0,15000000),breaks=50)
  # hist(QC_align$total)
  # plot(QC_align$PropPolyA12,QC_align$PropReadThrough)
  # plot(QC_align$PropPolyA12,QC_align$PropReadThrough,ylim=c(0,0.3))
  # dev.off()

  if (reads_cutoff==2){
  ## 2M cutoff: reads per library
  keep_QC<-QC_align$Library[which(QC_align$TotalReads>2e6 & QC_align$TotalReads<2e7 & QC_align$PropReadThrough<0.6 &QC_align$PropPolyA12<0.6)]
  } else if (reads_cutoff==1){
  ## 1M cutoff: reads per library
  keep_QC<-QC_align$Library[which(QC_align$TotalReads>1e6 & QC_align$TotalReads<2e7 & QC_align$PropReadThrough<0.6 &QC_align$PropPolyA12<0.6)]
  }

  keep_align<-QC_align$Library[which(QC_align$total>0.6)]

  if (rep==1){
    raw=raw1
    raw<-raw[,c(1,which((colnames(raw) %in% paste("X",keep_align,sep="")) & (colnames(raw) %in% paste("X",keep_QC,sep=""))))]
    ################################
    # change name to MLC standard, adjust names for swapped samples
    ################################
    #barcode_map<-read.table("/workdir/ml2498/OfficeCmp/GoogleDrive/TWAS_2018/TWAS/SD18_Barcode_MLC_map_rep1.txt",header=T,sep="\t",stringsAsFactors=F)
    barcode_map<-read.table("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/TWAS_2018/TWAS/SD18_Barcode_MLC_map_rep1.txt",header=T,sep="\t",stringsAsFactors=F)
    #barcode_map<-read.table("~/Desktop/GeneExpression/SD18_Barcode_MLC_map_rep1.txt",header=T,sep="\t",stringsAsFactors=F)
    barcode_map<-read.table("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/TWAS_2018/TWAS/SD18_Barcode_MLC_map_rep1.txt",header=T,sep="\t",stringsAsFactors=F)

    for (i in 2:ncol(raw)){
      colnames(raw)[i]<-barcode_map$MLC_mf[which(barcode_map$Barcode==colnames(raw)[i])]
    }

    plate_info<-as.data.frame(QC_align[,c("Library","Plate")])
    plate_info$Library<-as.character(plate_info$Library)
    for (i in 1:nrow(plate_info)){
      if(paste("X",plate_info$Library[i],sep="") %in% barcode_map$Barcode){
        plate_info$Library[i]<-barcode_map$MLC_mf[which(barcode_map$Barcode==paste("X",plate_info$Library[i],sep=""))]
      }
    }


    badnames<-c("LH123HT","TZU_CHIAO_HSI_WU_105","I29","YONG_28","CO256","CSJ3",
                "NC232","KO679Y","PHK76","N542","N501","A659","IA5125B","R229",
                "B8","MO17_244")

    adjnames<-c("TZU_CHIAO_HSI_WU_105.adj","LH123HT.adj","YONG_28.adj","I29.adj","CSJ3.adj","CO256.adj",
                "KO679Y.adj","NC232.adj","N542.adj","PHK76.adj","A659.adj","N501.adj","R229.adj","IA5125B.adj",
                "MO17_244.adj","B8.adj")

    for (i in 1:length(badnames)){
      colnames(raw)[which(colnames(raw)==badnames[i])]<-adjnames[i]
      plate_info$Library[which(plate_info$Library==badnames[i])]<-adjnames[i]
    }

    ### remove ".adj" from the adjusted names ####
    colnames(raw)<-gsub(".adj", "", colnames(raw))
    plate_info$Library<-gsub(".adj", "", plate_info$Library)
    NotMatchGBS<-c("MO17_110","SD102","A322","A305","K4","A321","KY21","W32","SG30A",
                   "H71","B79","CM99","B_18_INBR_FR_SUPERGOLD","G22_T122")

    raw<-raw[,-which(colnames(raw) %in% NotMatchGBS)]
  }else {
    raw=raw2
    raw<-raw[,c(1,which((colnames(raw) %in% paste("X",keep_align,sep="")) & (colnames(raw) %in% paste("X",keep_QC,sep=""))))]
    ################################
    # change name to MLC standard, adjust names for swapped samples
    ################################
    barcode_map<-read.table("/workdir/ml2498/OfficeCmp/GoogleDrive/TWAS_2018/TWAS/SD18_Barcode_MLC_map_rep2.txt",header=T,sep="\t",stringsAsFactors=F)
    #barcode_map<-read.table("~/Desktop/GeneExpression/SD18_Barcode_MLC_map_rep2.txt",header=T,sep="\t",stringsAsFactors=F)

    for (i in 2:ncol(raw)){
      colnames(raw)[i]<-barcode_map$MLC_mf[which(barcode_map$Barcode==colnames(raw)[i])]
    }

    plate_info<-as.data.frame(QC_align[,c("Library","Plate")])
    plate_info$Library<-as.character(plate_info$Library)
    for (i in 1:nrow(plate_info)){
      if(paste("X",plate_info$Library[i],sep="") %in% barcode_map$Barcode){
        plate_info$Library[i]<-barcode_map$MLC_mf[which(barcode_map$Barcode==paste("X",plate_info$Library[i],sep=""))]
      }
    }
    NotMatchGBS<-c("MO17_614","SG30A","B_18_INBR_FR_SUPERGOLD","W32","B114","MO7","F44")
    raw<-raw[,-which(colnames(raw) %in% NotMatchGBS)]
    }
  setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts")
  #setwd("~/Desktop/GeneExpression/v4_counts/")

  if (reads_cutoff==2){
    write.table(raw,paste("Hisat_cleanedQC_rawcounts_v4_rep",rep,".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
    write.table(plate_info,paste("RNAseq_PlateInfo_v4_rep",rep,".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
  } else if (reads_cutoff==1){
    write.table(raw,paste("Hisat_cleanedQC_rawcounts_v4_rep",rep,"_1M.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
    write.table(plate_info,paste("RNAseq_PlateInfo_v4_rep",rep,"_1M.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
  }
}

######################################################

#######################################################################
# normalize, filter, transform
#######################################################################
## Get scaled (sample-wise scaling, gene-wise scaling does not matter) cpm counts
library(data.table)
library(DESeq2)
reads_cutoff=2

for (rep in 1:2){
  if (reads_cutoff==2){
    #raw<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_cleanedQC_rawcounts_v4_rep",rep,".txt",sep=""),header=T,sep="\t")
    #raw<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_cleanedQC_rawcounts_v4_rep",rep,".txt",sep=""),header=T,sep="\t")
    raw<-fread(paste("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_cleanedQC_rawcounts_v4_rep",rep,".txt",sep=""),data.table=F)
    raw<-read.table(paste("/Users/Meng/Desktop/RNAseq_temp/Hisat_cleanedQC_rawcounts_v4_rep",rep,".txt",sep=""),header=T,sep="\t")

  }else if(reads_cutoff==1){
    #raw<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_cleanedQC_rawcounts_v4_rep",rep,"_1M.txt",sep=""),header=T,sep="\t")
    #raw<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/Hisat_cleanedQC_rawcounts_v4_rep",rep,"_1M.txt",sep=""),header=T,sep="\t")
  }

  raw<-raw[,-grep(".",colnames(raw),fixed=T)]  # 297 include MO17 for rep 1
  #colnames(raw)[grep(".",colnames(raw),fixed=T)]
  #colnames(raw)[which(colnames(raw)=="R229")]
  raw<-raw[-grep("_",raw[,1],fixed=T),]
  #raw[grep("_",raw[,1],fixed=T),1]

  rownames(raw)<-raw[,1]
  raw<-raw[,-1]
  raw<-as.matrix(raw)



  #test<-rlog(raw, blind=FALSE,betaPriorVar,fitType = "parametric")
  counts_rlog<-rlog(raw, blind=FALSE,fitType = "parametric")
  rownames(counts_rlog)<-rownames(raw)
  setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts")
  write.table(counts_rlog,paste("Hisat_cleanedQC_counts_v4_rep",rep,"_rlog.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
}

###### filtering and check for variability #####
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts")
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts")
setwd("/Users/Meng/Desktop/RNAseq_temp/")

## check for distribution of zeros
pdf("dist_num_zero_in_row_col.pdf")

for (rep in 1:2){
  ## row: gene; col: sample
  counts_rlog<-read.table(paste("Hisat_cleanedQC_counts_v4_rep",rep,"_rlog.txt",sep=""),head=T,row.names=1,sep="\t")
  #col_0<-apply(counts_rlog,2,function(c)sum(c=0))
  col_0<-colSums(counts_rlog == 0)/nrow(counts_rlog) # all accession has 13.384% zeros
  hist(col_0,main=paste("Distribution of number of zero values in columns (rep",rep,")",sep=""))
  #row_0<-rowSums(counts_rlog == 0)/ncol(counts_rlog)
  row_0<-rowSums(counts_rlog == 0) # *either the row is all zero's or no zero
  hist(row_0,main=paste("Distribution of number of zero values in rows (rep",rep,")",sep=""))

}
dev.off()

for (rep in 1:2){
  ## row: gene; col: sample
  counts_rlog<-read.table(paste("Hisat_cleanedQC_counts_v4_rep",rep,"_rlog.txt",sep=""),head=T,row.names=1,sep="\t")
  counts_rlog[counts_rlog<0]<-0

  nm_sample<-ncol(counts_rlog)
  row_0<-rowSums(counts_rlog == 0)
  counts_rlog.1<-counts_rlog[-which(row_0==nm_sample),]# remove row's with all zero's
  row_0.1<-rowSums(counts_rlog.1 == 0)
  hist(row_0.1,breaks=50)
  print(length(row_0.1[row_0.1==0]))
  print(length(row_0.1[row_0.1>250]))
  #table(row_0.1)

  ge_var<-apply(counts_rlog.1,1,var)
  print(range(ge_var)) # 1.172146e-10 2.595646e+01 rep1, 2.054193e-09 4.859525e+01 rep2
  print(median(ge_var)) # 0.1738482 rep1, 0.1487046 rep2

  flagged_genes<-rownames(counts_rlog.1)[which(row_0.1>250)]# flagged genes to check after TWAS and Fisher's test
  write.table(flagged_genes,paste("genes_w_0values_in250samples_rep",rep,".txt",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(counts_rlog.1,paste("Hisat_cleanedQC_counts_v4_rep",rep,"_rlog_clean.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
}

##############################################

# calculate pairwise correlations between genotypes
# R2<0.5 should be removed
# do not need scale based on library size for pearson correlation
###################################################
#setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts")
setwd("/Users/Meng/Desktop/RNAseq_temp/")

library(gplots)

rep=2
reads_cutoff=2
for(rep in 1:2){
  y.filtered<-read.table(paste("Hisat_cleanedQC_counts_v4_rep",rep,"_rlog_clean.txt",sep=""),header=T,sep="\t")
  geno_cor<-cor(y.filtered, method = "pearson", use = "complete.obs")
  is.matrix(geno_cor)

  Mo17_cor<-geno_cor[grep("MO17",rownames(geno_cor)),grep("MO17",colnames(geno_cor))]
  pdf(paste("Hisat_MO17_pw_cor_rep",rep,"_rlog.pdf",sep=""))
  heatmap.2(Mo17_cor, trace="none",margins=c(7,8))
  dev.off()

  pdf(paste("Hisat_sample_pw_cor_beforeClean_rep",rep,"_rlog.pdf",sep=""))
  heatmap.2(geno_cor, trace="none",labRow = FALSE, labCol = FALSE)

  heatmap.2(geno_cor, trace="none",labRow = FALSE, labCol = FALSE,
            Rowv=FALSE,Colv=FALSE)
  dev.off()
}

sqrt(0.5)
ColMeans<-colMeans(geno_cor, na.rm = FALSE, dims = 1)
ColMeans<-cbind(colnames(geno_cor),ColMeans)

## Nothing to remove for both reps
######### The END ###############
