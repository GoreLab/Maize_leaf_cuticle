st_clean<-vector()
for (rep in c(1,2)){
  #raw<-read.table(paste("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata/raw_data/stomatal_rep",rep,".txt",sep=""),header=T,sep="\t")
  raw<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/raw_data/stomatal_rep",rep,".txt",sep=""),header=T,sep="\t")
  raw<-raw[!is.na(raw$plot),]
  raw$plot<-as.character(raw$plot)
  raw$Stomate<-as.numeric(as.character(raw$Stomate))
  raw$Rep<-rep

  # create barcode
  for (j in 1:nrow(raw)){
    if (nchar(raw$plot[j])==3){
      raw$barcode[j]<-paste("X18SD",raw$plot[j],sep="")
    }else{
      raw$barcode[j]<-paste("X18SD",rep(0,(3-nchar(raw$plot[j]))),raw$plot[j],sep="")
    }
  }

  # at least 2 plants per plot, at least 2 stomata per plant
  raw.clean<-vector()
  Plots<-unique(raw$plot)

  #p=Plots[1]
  p="664"
  for (p in Plots){
    raw.temp<-raw[which(raw$plot==p),]
    Indv<-unique(raw.temp$plant)

    for (indv in Indv){
      if (nrow(raw.temp[which(raw.temp$plant==indv),])<2){
        raw.temp$plant[which(raw.temp$plant==indv)]<-NA
      }
    }

    if(nrow(raw.temp[is.na(raw.temp$plant),])>0){
      print(paste("Remove ",nrow(raw.temp[is.na(raw.temp$plant),])," rows from plot ",p," Indv ",indv,sep=""))
    }

    raw.temp<-raw.temp[!is.na(raw.temp$plant),]

    Indv.2<-unique(raw.temp$plant)# update the Individual list after filter for the number of stomata per plant

    if(length(Indv.2)<2){
      raw.temp$plot[which(raw.temp$plot==p)]<-NA
    }

    if (nrow(raw.temp[is.na(raw.temp$plot),])>0){
      print(paste("Remove ",nrow(raw.temp[is.na(raw.temp$plot),])," rows from plot ",p,sep=""))
    }
    raw.temp<-raw.temp[!is.na(raw.temp$plot),]

    raw.clean<-rbind(raw.clean,raw.temp)
  }

  st_clean<-rbind(st_clean,raw.clean)
}
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata")
write.table(st_clean,"cleaned_stomatal_measure_bothRep.txt",col.names=T,row.names=F,sep="\t",quote=F) # still by stomate

###############################################################################################
# outlier removal
# clean up by trait, if there is only one observation per Indv or one indiv per plot, set to NA for that trait (not doing this step)
# aggregate by plant, then outlier removal, to be consistent with wax data (plant/leaf level)
#################################################################################################
st_clean<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cleaned_stomatal_measure_bothRep.txt",header=T,sep="\t")
st_clean<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/stomata/cleaned_stomatal_measure_bothRep.txt",header=T,sep="\t")
st_clean<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cleaned_stomatal_measure_bothRep.txt",header=T,sep="\t")

### counting
st_clean$plot<-as.factor(st_clean$plot)
plots<-levels(st_clean$plot)
length(levels(st_clean$plot)) # 106 plots were phenotyped

MAX<-10
MIN<-10
NM_plant<-vector()
for (p in unique(st_clean$plot)){
  sub_st<-st_clean[which(st_clean$plot==p),]
  nm_plant<-length(unique(sub_st$plant)) # number of plants in each plot
  NM_plant<-c(NM_plant,nm_plant)

  nm_st<-table(sub_st$plant)
  print(nm_st)
  if(max(nm_st)>=MAX){
    MAX<-max(nm_st)
  }
  if(min(nm_st)<=MIN){
    MIN<-min(nm_st)
  }
}
#
MAX ## 17, maximum number of stomata phenotyped for one plant
MIN ## 2, minimum number of stomata phenotyped for one plant
range(NM_plant) # 2-5 plants per plot were phenotyped


### aggregate
agg = aggregate(.~ barcode+ plant,data=st_clean,FUN = mean,na.action = na.omit,na.rm = TRUE)

colnames(agg)[1]<-"Barcode"

#barcode<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/wax/CE_outliers_wax.txt",header=T,sep="\t")
barcode<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/wax/CE_outliers_wax.txt",header=T,sep="\t")
agg<-merge(barcode[,1:2],agg,by="Barcode",all.y=T)

##########################################
# check for outliers
##########################################
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/")
setwd("/home/ml2498/Desktop/Labserver2/MaizeLeafCuticle/TWAS_2018/stomata")

agg$Rep<-as.factor(agg$Rep)

library(ggplot2)
#setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata/plot")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/plot")
pdf('Check_Outliers_byPlot.pdf',height=4,width=10)
for (i in c(6:9)){
  print(ggplot(agg, aes(x = MLC_STANDARD, y = agg[,i], group = Rep))+
          geom_point(aes(color=Rep))+
          theme(axis.text.x = element_text(angle = 90, hjust = 1))+
          labs(y = colnames(agg)[i])
  )
}
dev.off()

################################################
# outlier removal
################################################
source("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/script/outlier_removal_function_Meng.R")
#source("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/script/outlier_removal_function_Meng.R")

agg<-cbind(agg[,c(1:4,10,6:9)])

library(asreml)

nbvariables=5
nbtraits=4
random=c("Rep")
fixed=c("MLC_STANDARD")
BLUx="BLUE"

clean.dataset <- agg[,c(1:nbvariables)]
alltraitnames <- vector()

for (i in 1:nbtraits){
  curr.trait <- colnames(agg[ (nbvariables + i) ])
  alltraitnames <- c(alltraitnames, curr.trait)

  transfpheno=agg[,c(1:nbvariables,i+nbvariables)]
  initial.outliers(transfpheno, curr.trait, random, fixed, nbvariables, BLUx) -> out
  as.matrix(out) -> out
  colnames(out) <- curr.trait

  if(i == 1) {
    cleanedpheno <- cbind(clean.dataset, out)
  } else {
    cleanedpheno <- cbind(cleanedpheno, out)
  }

}

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/stomata")
write.table(cleanedpheno,"all_stomatal_byLeaf_outlierRM.txt",col.names=T,row.names=F,sep="\t",quote=F)
###########################
# plot outlier removed data
###########################
library(ggplot2)
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata")

cleanedpheno<-read.table("all_stomatal_byLeaf_outlierRM.txt",header=T,sep="\t")
cleanedpheno$Rep<-as.factor(cleanedpheno$Rep)

setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata/plot")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/plot")
pdf('Check_after_outlierRM_byPlot.pdf',height=4,width=10)
for (i in 6:9){
  print(ggplot(cleanedpheno, aes(x = MLC_STANDARD, y = cleanedpheno[,i], group = Rep))+
          geom_point(aes(color=Rep))+
          theme(axis.text.x = element_text(angle = 90, hjust = 1))+
          labs(y = colnames(cleanedpheno)[i])
  )
}
dev.off()

pdf('Stomata_Distrn_after_outlierRM_leafBased.pdf',height=4,width=5)
for (i in 6:9){
  hist(cleanedpheno[,i],xlab=colnames(cleanedpheno)[i],main="")
}
dev.off()
#################################
# aggregate by plot
# then aggregate by taxa
#################################
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/stomata")

cleanedpheno<-read.table("all_stomatal_byLeaf_outlierRM.txt",header=T,sep="\t")
byPlot<-aggregate(.~ Barcode,data=cleanedpheno,FUN = mean,na.action = na.omit,na.rm = TRUE)

byPlot<-byPlot[,-(2:3)] # remove the MLC_STANDARD column, need to rematch

barcode<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/wax/CE_outliers_wax.txt",header=T,sep="\t")
#barcode<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/wax/CE_outliers_wax.txt",header=T,sep="\t")
byPlot<-merge(barcode[,1:2],byPlot,by="Barcode",all.y=T)

# may need to update this for wax data
library(asreml)

st_Blup<-vector()
for (t in 5:8){

  fit_as<-eval(parse(text=paste("asreml(fixed = ",colnames(byPlot)[t]," ~ 1,random= ~Rep+MLC_STANDARD,na.action=na.method(x='omit'),data = byPlot)"))) #asreml 4
  #fit_as<-eval(parse(text=paste("asreml(fixed = ",colnames(byPlot)[t]," ~ 1,random= ~Rep+MLC_STANDARD,na.method.X='omit',data = byPlot)"))) #asreml 3

  blup<-fit_as$coefficients$random
  blup<-blup[-1,]

  intercept<-fit_as$coefficients$fixed[1]
  blup<-round(blup+intercept,5)
  st_Blup<-cbind(st_Blup,blup)

}
st_Blup<-cbind.data.frame(rownames(st_Blup),st_Blup)
colnames(st_Blup)<-c("MLC_STANDARD",colnames(byPlot)[5:8])
st_Blup$MLC_STANDARD<-as.character(st_Blup$MLC_STANDARD)
st_Blup$MLC_STANDARD<-substr(st_Blup$MLC_STANDARD,14,nchar(st_Blup$MLC_STANDARD))

write.table(st_Blup,"stomata_BLUPs_tobeformatted.txt",col.names = T,row.names = F,quote=F,sep="\t")
