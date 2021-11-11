library(data.table)

gc<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/NAM_NIL/TWAS/nNIL_gc_BLUE_sas_combEnv.txt",header=T,sep="\t")
colData<-gc[grep("B73",gc$Genotype,fixed=T),]
colData$Genotype<-as.character(colData$Genotype)
##
colData$ISTL1_intro<-"B73"

for(i in 1:nrow(colData)){
  if(colData$Genotype[i]!="B73"){
    temp<-strsplit(colData$Genotype[i],"-",fixed=T)
    colData$ISTL1_intro[i]<-temp[[1]]
  }
}
########################
colData$block85<-"haplotype_1"
colData$block94<-"haplotype_1"
colData$block82<-"haplotype_1"
colData$block42<-"haplotype_1"

colData$block79<-"haplotype_1"
colData$block74<-"haplotype_1"
colData$block90<-"haplotype_1"


## block85
for (i in 1:nrow(colData)){
  if(colData$ISTL1_intro[i] %in% c("CML69","Ki3","NC350","Oh43")){
    colData$block85[i]<-"haplotype_6"
  }else if(colData$ISTL1_intro[i] %in% c("Ki11","NC358","Tx303")){
    colData$block85[i]<-"haplotype_7"
  }else if(colData$ISTL1_intro[i] %in% c("Tzi8")){
    colData$block85[i]<-"haplotype_2"
  }
}
colData$block85_collapse<-colData$block85
colData$block85_collapse[which(colData$ISTL1_intro=="Tzi8")]<-"haplotype_1"
## block94
for (i in 1:nrow(colData)){
  if(colData$ISTL1_intro[i]=="Ki11"){
    colData$block94[i]<-"haplotype_4"
  }else if(colData$ISTL1_intro[i]=="Tx303"){
    colData$block94[i]<-"haplotype_2"
  }
}
colData$block94_snp<-"A"
colData$block94_snp[which(colData$ISTL1_intro=="Tx303"|colData$ISTL1_intro=="Ki11")]<-"C"

## block82
for (i in 1:nrow(colData)){
  if(colData$ISTL1_intro[i] %in% c("Ki11","NC358","Tx303")){
    colData$block82[i]<-"haplotype_7"
  }else if(colData$ISTL1_intro[i] %in% c("Tzi8")){
    colData$block82[i]<-"haplotype_5"
  }
}
colData$block82_collapse<-colData$block82
colData$block82_collapse[which(colData$ISTL1_intro=="Tzi8")]<-"haplotype_1"

## block42
for (i in 1:nrow(colData)){
  if(colData$ISTL1_intro[i] %in% c("Ki11","NC358","Tx303")){
    colData$block42[i]<-"haplotype_2"
  }else if(colData$ISTL1_intro[i] %in% c("Oh43")){
    colData$block42[i]<-"haplotype_7"
  }else if(colData$ISTL1_intro[i] %in% c("Tzi8")){
    colData$block42[i]<-"haplotype_3"
  }
}
##block79
for (i in 1:nrow(colData)){
  if(colData$ISTL1_intro[i] %in% c("Ki11", "NC358", "Tx303")){
    colData$block79[i]<-"haplotype_8"
  }
}
## block74
for (i in 1:nrow(colData)){
  if(colData$ISTL1_intro[i] %in% c("B73","CML69","Ki3","NC350","Oh43","Tzi8")){
    colData$block74[i]<-"haplotype_2"
  }
}

## block90
for (i in 1:nrow(colData)){
  if(colData$ISTL1_intro[i] %in% c("Ki11","NC358","Tx303")){
    colData$block90[i]<-"haplotype_2"
  }
}

setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/NAM_NIL/gc_hap_assotiation")
write.table(colData,"nNIL_gc_combENV_10top_haplotype_20SD_SAS_BLUE.txt",row.names=F,col.names=T,quote=F,sep="\t") #BLUE, imcomplete block design,asreml

######################################
library(asreml)
library(asremlPlus)

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/NAM_NIL/gc_hap_assotiation")
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/NAM_NIL/gc_hap_assotiation")
colData<-read.table("nNIL_gc_combENV_10top_haplotype_20SD_SAS_BLUE.txt",header=T,sep="\t")
## gc ~ haplotype
sink(file="pairwise_top10_haplo_effect_gcBLUE_combENV_SD20_asreml.txt")
for (i in 4:ncol(colData)){
  #for (i in c(10:12)){
  colData[,i]<-as.factor(as.character(colData[,i]))
  if(length(levels(colData[,i]))>1){
    fit<-eval(parse(text=paste("asreml(fixed=gc_BLUE~",colnames(colData)[i],",data=colData,na.method.X='omit')",sep="")))
    print(paste("**** Testing genetic effect of ", colnames(colData)[i]," ****",sep=""))
    fit_t <- as.asrtests(fit, NULL, NULL)
    diffs <- predictPlus.asreml(classify = colnames(colData)[i], pairwise=TRUE,
                                asreml.obj = fit, 
                                wald.tab = fit_t$wald.tab)
    print(fit_t$wald.tab)
    #(BLUPs_n<-diffs$differences)  
    #(BLUPs_n_p<-diffs$p.differences)
    
  }
}
closeAllConnections()
############################

setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/NAM_NIL/gc_hap_assotiation")
setwd("/Users/zhenghao/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/NAM_NIL/gc_hap_assotiation")
allelic_effect<-read.table("allelic_effect_comp_widiv_vs_nNIL.txt",header=T,sep="\t")

allelic_effect$haplotype<-as.character(allelic_effect$haplotype)
allelic_effect$haploblock<-as.character(allelic_effect$haploblock)

blocks<-unique(allelic_effect$haploblock)
pdf("allelic_effect_comp_widiv_vs_nNIL.pdf",width=9, height=3,family="serif")
par(mfrow=c(1,3))
for(b in blocks){
  temp<-allelic_effect[which(allelic_effect$haploblock==b),]
  plot(temp$widiv,temp$nNIL,xlab=bquote("Haplotype effect (g "~h^-1~g^-1~"; WiDiv panel)"),ylab="Haplotype effect (g"~h^-1~g^-1~"; nNILs)",main=b,
       xlim=c(range(temp$widiv)[1]-0.005,range(temp$widiv)[2]+0.005),
       ylim=c(range(temp$nNIL)[1]-0.005,range(temp$nNIL)[2]+0.005)
  )
  text(temp$nNIL~temp$widiv, labels=temp$haplotype, cex=0.5, font=2, pos=3,
       xlim=c(range(temp$widiv)[1]-0.005,range(temp$widiv)[2]+0.005),
       ylim=c(range(temp$nNIL)[1]-0.005,range(temp$nNIL)[2]+0.005)
       )
  fit=lm(temp$nNIL~temp$widiv)
  abline(fit,col="red")
  r<-format(cor(temp$widiv,temp$nNIL,use="complete.obs"), digits=3)
  legend("topleft",bty="n", legend=bquote(italic(r)~"="~.(r)),col="red")
  
}
dev.off()

pdf("allelic_effect_comp_widiv_vs_nNIL_block52.pdf",width=5, height=4.5,family="sans")
b="block82"
  temp<-allelic_effect[which(allelic_effect$haploblock==b),]
  plot(temp$widiv,temp$nNIL,xlab=bquote("Haplotype effect (g "~h^-1~g^-1~"; WiDiv panel)"),ylab="Haplotype effect (g"~h^-1~g^-1~"; nNILs)",main=NULL,
       xlim=c(range(temp$widiv)[1]-0.005,range(temp$widiv)[2]+0.005),
       ylim=c(range(temp$nNIL)[1]-0.005,range(temp$nNIL)[2]+0.005)
  )
  # text(temp$nNIL~temp$widiv, labels=temp$haplotype, cex=0.8, font=2, pos=3,
  #      xlim=c(range(temp$widiv)[1]-0.005,range(temp$widiv)[2]+0.005),
  #      ylim=c(range(temp$nNIL)[1]-0.005,range(temp$nNIL)[2]+0.005)
  #)
  fit=lm(temp$nNIL~temp$widiv)
  abline(fit,col="red")
  r<-format(cor(temp$widiv,temp$nNIL,use="complete.obs"), digits=3)
  legend("topleft",bty="n", legend=bquote(italic(r)~"="~.(r)),col="red")

dev.off()





