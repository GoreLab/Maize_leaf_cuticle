cd /workdir/ml2498/hap_gwas #robbins' server
scp cbsugore02:/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window/robbins_server/* .

cd /workdir/ml2498/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window/robbins_server
mkdir results
cd results
scp cbsurobbins:/workdir/ml2498/hap_gwas/* .

### format genotype, phenotype and kinship matrix for ASReml GWAS

library(data.table)

#setwd("/workdir/ml2498/hap_gwas")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/200K_window")
system("ls")
## no filtering on haplotype frequency
#hap<-fread("/workdir/ml2498/hap_gwas/HMP3_hap_chr4_istl1_region_600K.hmp.txt",data.table=F)
hap<-fread("../200K_window/HMP3_hap_chr4_istl1_region_400K.hmp.txt",data.table=F)
#geno<-fread("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/HMP3_geno_chr4_istl1_region_3e+05.hmp.txt",data.table=F)

hap.in<-t(hap)

rownames(hap.in)[1]<-"rs"

hap.info<-t(hap.in[c(1:4,12,14),])
colnames(hap.in)<-as.character(hap.info[,1])
hap.in<-hap.in[-(1:14),]

hap.in<-hap.in[order(rownames(hap.in)),]

# for(i in 1:ncol(hap.in)){
#   hap.in[,i]<-as.numeric(as.character(hap.in[,i]))
# }
#
# for (i in 1:ncol(hap.in)){
#   if(any(hap.in[,i]>10,na.rm=T)){
#     print(i)
#   }
# }
#
# hap.in<-hap.in*0.1

write.table(hap.in,"HMP3_hap_chr4_istl1_region_400K.txt",quote=F,sep="\t")
write.table(hap.info,"HMP3_hap_chr4_istl1_region_400K.info.txt",row.names=F,quote=F,sep="\t")

#############
# phenotype
#############
#pheno.all<-read.table("CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
pheno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
#Taxa<-read.table("MLC_GBSSNP310_468K_v4_chr1_prun09.012.indv",stringsAsFactors=F)
Taxa<-read.table("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/filter_prune_468K/MLC_GBSSNP310_468K_v4_chr1_prun09.012.indv",stringsAsFactors=F)

pheno<-pheno.all[which(pheno.all$MLC_STANDARD %in% Taxa[,1]),c(1,23)]
pheno<-pheno[order(as.character(pheno$MLC_STANDARD)),]
which(pheno$MLC_STANDARD!=rownames(hap.in))

write.table(pheno,"Pheno310_ce_sd18_both_untr_sorted.txt",row.names=F,quote=F,sep="\t")

################
# kinship
################
#myKI<-read.table("centeredIBS_HMP3_AGPv4_LD02_450to310.txt")
myKI<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt")

rownames(myKI)<-myKI[,1]
myKI<-myKI[,-1]
colnames(myKI)<-rownames(myKI)
myKI<-myKI[order(rownames(myKI)),order(colnames(myKI))]
which(rownames(myKI)!=pheno[,1])
myKI<-cbind.data.frame(rownames(myKI),myKI)

write.table(myKI,"centeredIBS_HMP3_AGPv4_LD02_450to310_sort.txt",col.names=F,row.names=F,quote=F,sep="\t")


###################
# gwas using asreml
###################
library(asreml)
library(MASS)

#setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window/robbins_server/results")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/200K_window")
pheno<-read.table("Pheno310_ce_sd18_both_untr_sorted.txt",header=T,sep="\t")

kin<-read.table("centeredIBS_HMP3_AGPv4_LD02_450to310_sort.txt",header=F,sep="\t",row.names=1)
colnames(kin)<-rownames(kin)

geno<-read.table("HMP3_hap_chr4_istl1_region_400K.txt",header=T,sep="\t",comment.char="")
#geno$interval<-geno$EndPos-geno$pos
#hist(interval)

rownames(pheno)<-pheno[,1]
kin<-as.matrix(kin)
kin.inv<-ginv(kin)

#Row<-rep(1:432,432)
#Column<-rep(1:432,each=432)
#Row<-rep(rownames(kin),432)
#Column<-rep(rownames(kin),each=432)
#Ainverse<-stack(as.data.frame(kin.inv))
#kin.ainv<-cbind(Row,Column,Ainverse[,1])
#rownames(kin.ainv)<-rep(rownames(kin),432)
#colnames(kin.ainv)[3]<-"Ainverse"

toSparse <- function(m) {
  comb <- data.frame(row = rep(1:nrow(m), each = nrow(m)),
                     column = rep.int(1:nrow(m), nrow(m)))
  x <- comb[comb$row >= comb$column, ]
  x$value <- m[cbind(x$row, x$column)]
  attr(x, 'rowNames') <- rownames(m)
  return(x)
}
kin.ainv<-toSparse(kin.inv)
attr(kin.ainv,'rowNames')<-rownames(kin)

#rownames(kin.inv)<-rownames(kin)
#colnames(kin.inv)<-colnames(kin)

hap.geno<-geno

hap.info<-read.table("HMP3_hap_chr4_istl1_region_400K.info.txt",header=T,sep="\t")

#rownames(hap.geno)<-hap.geno[,1]
which(rownames(hap.geno)!=pheno[,1])
hap.geno<-cbind(pheno,hap.geno)
#attr(hap.geno,'rowNames')<-rownames(kin)

p<-vector()
allele<-vector()
Effect<-list()

hap.geno$MLC_STANDARD<-as.factor(hap.geno$MLC_STANDARD)
for (j in 3:ncol(hap.geno)){

  # ## if run with asreml v4
  # tryCatch({
  #   fit.asr <- asreml(fixed = ce_SD18_both_untr ~ 1+as.factor(hap.geno[,j]),random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),na.action=na.method(x='omit'),data = hap.geno)
  #   fit.asr <- asreml(fixed = ce_SD18_both_untr ~ 1+as.factor(hap.geno[,j]),random= ~vm(MLC_STANDARD,kin.ainv),na.action=na.method(x='omit'),data = hap.geno)
  #
  # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # ##################################


  #for (j in 3:10){
  hap.fit<-asreml(ce_SD18_both_untr~1+as.factor(hap.geno[,j]),random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=hap.geno,na.method.X='omit')
  #hap.fit<-asreml(ce_ALL4_REML_tr~1+hap.geno[,j],random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=hap.geno,na.method.X='include')
  #hap.fit<-asreml(ce_ALL4_REML_tr~as.factor(hap.geno[,j]),data=hap.geno,na.method.X='include')
  #print(hap.fit$coefficients$fixed)
  #print(anova(hap.fit))
  effect=hap.fit$coefficients$fixed
  Effect[[colnames(hap.geno)[j]]]<-effect
  fav_allele=names(effect)[which(effect==min(effect[-length(effect)]))]

  #Wald_test<-wald(hap.fit, denDF = "default", ssType = "conditional")
  #pvalue<-Wald_test$Wald[2,"Pr"]

  pvalue=anova(hap.fit)$`Pr(Chisq)`[2]
  allele<-c(allele,fav_allele)
  p<-c(p,pvalue)
}
results<-cbind(hap.info,allele,p)
write.table(results,paste("ISTL1_600K_hap_gwas_results_noFilter.txt",sep=""),col.names = T,row.names = F,sep="\t",quote=FALSE)
lapply(Effect, write, "Allelic_effect_ISTL1_600K_hap_gwas_results_noFilter.txt", append=TRUE, ncolumns=1000)

for (i in 3:ncol(hap.geno)){
  hap.geno[,i]<-as.factor(hap.geno[,i])
}

imp_hap<-c(95,82,76,90,94,42,33,79,74,85) # important haploblock numbers
sink(file="pairwise_haplo_effect_gc_GWAS_310_asreml.txt")
for (hap in imp_hap){
  j=hap+2
  print(paste("**** Testing genetic effect of Block ", hap," ****",sep=""))
  hap.fit<-asreml(ce_SD18_both_untr~1+hap.geno[,j],random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=hap.geno,na.method.X='omit')
  hap.fit_t <- as.asrtests(hap.fit, NULL, NULL)
  diffs <- predictPlus.asreml(classify = "hap.geno[, j]", pairwise=TRUE,
                              asreml.obj = hap.fit,
                              wald.tab = hap.fit_t$wald.tab)

}
closeAllConnections()


#########################################
# Manhattan plot
#########################################
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window/robbins_server/results")
setwd("/Users/zhenghao/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window/robbins_server/results")

library(data.table)
NAM_hap<-fread("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window/NAM_founders_istl1_600K_haps_9999Prun_imp.hmp.txt",
        header=T,sep="\t",data.table = F)
NAM_hap<-fread("/Users/zhenghao/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window/NAM_founders_istl1_600K_haps_9999Prun_imp.hmp.txt",
               header=T,sep="\t",data.table = F)
nam_hap<-NAM_hap[,-(24:31)]

poly_in_nam<-c()
for (i in 1:nrow(nam_hap)){

  single_hap<-t(nam_hap[i,-(1:14)])
  stat_haps<-table(single_hap)
  nm_haps<-length(stat_haps)
  if(nm_haps>1){
    poly_in_nam<-c(poly_in_nam,as.character(nam_hap[i,1]))
  }
}

file<-"ISTL1_600K_hap_gwas_results_noFilter"

##

hap.res<-read.table(paste(file,".txt",sep=""),header=T,sep="\t")
## subset the results for 200Kb window
hap.res<-hap.res[which(hap.res$block_number<=121 & hap.res$block_number>30),]
hap.res$block_number<-1:nrow(hap.res)
hap.res$p[1]<-0.0084
#######################################

hap.res$logp<-(-1)*log10(hap.res$p)
hap.res$poly_in_nam<-"N"
hap.res$poly_in_nam[which(hap.res$rs %in% poly_in_nam)]<-"Y"
hap.res$pos_Mb<-hap.res$pos/1000000
hap.res$EndPos_Mb<-hap.res$EndPos/1000000

setwd("/Users/zhenghao/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/200K_window/")
write.table(hap.res,"ISTL1_400K_hap_gwas_results_noFilter_forPlot.txt",row.names=F,quote=F,sep="\t")

a<-ceiling(max(hap.res$logp))
b<-range(hap.res$pos_Mb)[1]
c<-range(hap.res$EndPos_Mb)[2]

istl1<-c(31729260,31732664)/1000000
#LD_window<-c(31568471,31817229)
#LD_200K<-c(31533883,31933883)/1000000

# hap.res_200K<-hap.res[which(hap.res$pos_Mb<=LD_200K[2]&hap.res$EndPos_Mb>=LD_200K[1]),]
# hap.res_200K$fdr<-p.adjust(hap.res_200K$p, method ="fdr") #fdr 0.05: 37 for MAF 0.05, 42 for MAF 0.03, 37 for MAF 0.04
# sig_hap<-hap.res_200K[which(hap.res_200K$fdr<0.05),c(1,6)]
# setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/NAM_NIL/haploview/300K_window")
# write.table(sig_hap,"sig_hap_200K_window_ByFDR.txt",row.names=F,col.names=T,sep="\t",quote=F)

#hap.res$fdr200K<-"N"
#hap.res$fdr200K[which(hap.res$rs %in% sig_hap[,1])]<-"Y"

#pdf(paste(file,".pdf",sep=""),width=7, height=4.5)
pdf(paste("ISTL1_400K_hap_gwas_results_noFilter.pdf",sep=""),width=7, height=4.5,family="sans")
#pdf(paste(SNP,"_v3_cand2.pdf",sep=""),width=7, height=4.5)
par(mar = c(5,5,2,5))

plot(NULL, xlim=c(b,c), ylim=c(0,a), ylab=expression(-log[10](italic(p))), xlab="Physical position (Mb) on chromosome 4")
#rect(xleft=LD_200K[1], ybottom=0-0.2, xright=LD_200K[2], ytop=a+0.2, density = NULL,col="lightyellow1",border=NA)
rect(xleft=istl1[1], ybottom=0-0.2, xright=istl1[2], ytop=a+0.2, density = NULL,col="grey85",border=NA)
for (i in 1:nrow(hap.res)){
  segments(x0=hap.res$pos_Mb[i],x1=hap.res$EndPos_Mb[i],y0=hap.res$logp[i],y1=hap.res$logp[i],
           #col=ifelse(logp>=threshold,"blue","grey"),
           col=ifelse(hap.res$poly_in_nam[i]=="Y","blue2","deeppink3"),
           lty=1,ylim=c(0,a),
           #lwd=ifelse(logp>threshold1,1,0.75)
           lwd=2)

}
#points(hap.res$logp~hap.res$pos,col="blue",pch=1,cex=1)
dev.off()

hap.res$fdr<-p.adjust(hap.res$p, method ="fdr") #fdr 0.05: 37 for MAF 0.05, 42 for MAF 0.03, 37 for MAF 0.04
