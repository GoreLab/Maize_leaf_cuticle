suppressMessages(library(rpart))
suppressMessages(library(ggplot2))
suppressMessages(library(mlbench))
suppressMessages(library(caret))
suppressMessages(library(party))
suppressMessages(library(pROC))
suppressMessages(library(pROC))

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

# wax_ce<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata/CE_stomata_wax_PC_Hansey_FS.txt",header=T,sep="\t")
#wax_ce<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/CE_stomata_wax_PC_Hansey_FS.txt",header=T,sep="\t")
wax_ce<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/stomata/CE_stomata_wax_PC_Hansey_FS.txt",header=T,sep="\t")

#OL<-c("NP87","SG_18","SG_1533","NC472")# 2 popcorn + 1 flint (genetic outliers); last one is tropical
OL<-c("NP87","SG_18","SG_1533","FC46","N501","NC472")# 2 popcorn + 3 flint (genetic outliers); last one is tropical

data<-wax_ce[-which(wax_ce$MLC_STANDARD %in% OL),]
##### handle missing data
for (i in 1:nrow(data)){
  len<-length(data[i,is.na(data[i,])])
  if (len>0){
    print(paste(i,": ",len,sep=""))
  }
}
data<-data[-c(52,53),-55]#AC_UK3: col 55
#data<-data[,-54]

#### keep all ACs
data0<-data

### Plot: simple regression

wax_name<-read.table("/Users/Meng/Desktop/LabServer/OfficeCmp/GoogleDrive/MLC_AZ_2017/uplift_to_v4/wax_GeneExp/raw_wax/wax_cutin_gradient.txt",header=T,sep="\t")
wax_name<-as.character(wax_name[c(1:18,20:22,24,26,28:35,38:51,55:57,60:61,64:75),1])
wax_name[48]<-"Alicyclic Campesterol"
wax_name
trait_names<-c(wax_name,colnames(data)[c(68:70)])

results<-vector()
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata/plot")
pdf("wax_stomata_ce_naiveCor_noPFT.pdf")
par(mar=c(6,6,4,4))
index<-1
for (i in c(5:66,68:70)){

  fit<-cor.test(data[,i],data$ce_SD18_both_untr,use="completed.obs")
  plot(data[,i],data$ce_SD18_both_untr,xlab=trait_names[index],ylab=expression(italic(g)[c](gh^{-1}*g^{-1})),col=data$tail)

  #legend("right",legend=unique(data$tail),pch=1,col=as.numeric(as.factor(levels(data$tail))))
  fit1<-lm(data$ce_SD18_both_untr~data[,i])
  abline(fit1,col="red",lwd=2)
  legend("topright", bty="n", legend=paste("r =",format(fit$estimate, digits=3)))
  res<-c(fit$estimate,fit$p.value,summary(fit1)$r.squared)
  results<-rbind(results,res)
  index=index+1
}
dev.off()

results<-as.data.frame(results)
results$R2<-results[,1]^2
colnames(results)<-c("r","P-value","R2_regression","r*r")
rownames(results)<-trait_names
write.table(results,"wax_stomata_ce_cor_naive.txt",col.names=T,row.names = T,sep="\t",quote=F)


############################

#################################
# pick mtry=47 based on 5-fold CV*****
#################################

wax_ce<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/stomata/CE_stomata_wax_PC_Hansey_FS.txt",header=T,sep="\t")

#OL<-c("NP87","SG_18","SG_1533","NC472")# 2 popcorn + 1 flint (genetic outliers); last one is tropical
OL<-c("NP87","SG_18","SG_1533","FC46","N501","NC472")# 2 popcorn + 3 flint (genetic outliers); last one is tropical

data<-wax_ce[-which(wax_ce$MLC_STANDARD %in% OL),]
##### handle missing data
for (i in 1:nrow(data)){
  len<-length(data[i,is.na(data[i,])])
  if (len>0){
    print(paste(i,": ",len,sep=""))
  }
}
data<-data[-c(52,53),-55]#AC_UK3: col 55

#### keep all ACs
data0<-data

#####
one_four<-floor(nrow(data)/5)
five=nrow(data)-4*one_four
## to repeat, start from here
sample_tool<-c(rep(1:4,each=one_four),rep(5,five))

SAMPLE_tool<-c(1:length(sample_tool))
set.seed(999)
for (l in 1:50){
  sample_tool<-sample(sample_tool)
  SAMPLE_tool<-cbind.data.frame(SAMPLE_tool,sample_tool)
}
colnames(SAMPLE_tool)<-c("index",paste("iteration",1:50,sep=""))


for(mtry in c(10,20,30,40,50,60)){
  for(ntree in c(500,1000,2000)){
    Fitted.TE<-vector()
    IMP<-colnames(data0)[c(5:66,68:70)]

    for(l in 1:50){

      sample_tool<-SAMPLE_tool[,(l+1)]
      data_wk<-as.data.frame(cbind(data0,sample_tool))

      fitted.TE<-vector()
      Imp<-colnames(data0)[c(5:66,68:70)]
      for (j in 1:5){
        data.te<-data_wk[which(data_wk$sample_tool==j),]
        data.tr<-data_wk[-which(data_wk$sample_tool==j),]

        fml<-paste("cforest(ce_SD18_both_untr~",paste(colnames(data.tr)[c(5:66,68:70)],collapse="+"),
                   ",data=data.tr,controls=cforest_unbiased(ntree=",ntree,",mtry=",mtry,"))",sep="")
        model<-eval(parse(text=fml))
        #model

        #fitted.tr <- predict(model, subset(data.tr,select=c(5:67,69:71)), OOB=TRUE, type= "response")
        fitted.te <- predict(model, newdata= subset(data.te,select=c(5:66,68:70)),type='response')## not correct

        if(l==1){
          combine<-cbind.data.frame(data.te[,c(1,3)],fitted.te)
          fitted.TE<-rbind.data.frame(fitted.TE,combine)
        }else{
          combine<-cbind.data.frame(data.te$MLC_STANDARD,fitted.te)
          colnames(combine)[1]<-"MLC_STANDARD"
          fitted.TE<-rbind.data.frame(fitted.TE,combine)

        }
        imp<-varimp(model)
        Imp<-cbind.data.frame(Imp,imp)

      }# loop of 5-fold cross validation

      ## predicted gc
      if(l==1){
        Fitted.TE<-fitted.TE
      }else{
        Fitted.TE<-merge(Fitted.TE,fitted.TE,by="MLC_STANDARD",all=T) # accumulate fitted values from all five folds
      }

      ## importance score
      for (m in 2:6){
        Imp[,m]<-as.numeric(as.character(Imp[,m]))
      }

      ## mean of importance score in 5-fold CV
      RF.imp.1<-apply(Imp[,-1],1,mean)
      IMP<-cbind.data.frame(IMP,RF.imp.1)

    }
    colnames(Fitted.TE)<-c("MLC_STANDARD","ce_SD18_both_untr",paste("predicted_",1:50,sep=""))
    colnames(IMP)<-c("Variable",paste("iteration_",1:50,sep=""))

    setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
    write.table(Fitted.TE,paste("predicted_gc_5fCV_51lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),row.names=F,sep="\t",quote=F)
    write.table(IMP,paste("importanceScore_gc_5fCV_51lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),row.names=F,sep="\t",quote=F)
  }
}

# calculate RMSE and R2 for each mtry
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")

ACCU<-c(paste("iteration_",1:50,sep=""))
for(mtry in c(10,20,30,40,50,60)){
  for(ntree in c(500,1000,2000)){

    Fitted.TE<-read.table(paste("predicted_gc_5fCV_51lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),header=T,sep="\t")


    RMSE.te<-vector()
    R.te<-vector()

    for (k in 3:ncol(Fitted.TE)){
      rmse.te<-RMSE(Fitted.TE$ce_SD18_both_untr,Fitted.TE[,k])
      RMSE.te<-c(RMSE.te,rmse.te)

      cor.te<-cor(Fitted.TE[,k],Fitted.TE$ce_SD18_both_untr)
      R.te<-c(R.te,cor.te)
    }
    accuracy<-cbind.data.frame(RMSE.te,R.te)
    colnames(accuracy)<-c(paste("RMSE_mtry",mtry,"_ntree",ntree,sep=""),paste("r_mtry",mtry,"_ntree",ntree,sep=""))

    ACCU<-cbind(ACCU,accuracy)

  }
}
ACCU1<-ACCU[,-1]
rownames(ACCU1)<-ACCU[,1]
Mean<-apply(ACCU1,2,mean)
SD<-apply(ACCU1,2,sd)
ACCU1<-rbind.data.frame(ACCU1,Mean,SD)
rownames(ACCU1)[c(51,52)]<-c("mean","sd")
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
write.table(ACCU1,"accuracy_gc_5fCV_51lines_mtry_ntree_combination_cforest.txt",row.names=T,sep="\t",quote=F)

########################

##############################

## Importance
library(reshape2)
library(ggplot2)

setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
wax_name<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/short_vs_complete_wax_names.txt",header=T,sep="\t")

setwd("/Users/zhenghao/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
wax_name<-read.table("/Users/zhenghao/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/short_vs_complete_wax_names.txt",header=T,sep="\t")

for(mtry in c(10,20,30,40,50,60)){
  for(ntree in c(500,1000,2000)){
    IMP<-read.table(paste("importanceScore_gc_5fCV_51lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),header=T,sep="\t")
    IMP<-merge(wax_name,IMP,by.x="short_names",by.y="Variable",all.y=T)
    IMP<-IMP[,-1]
    IMP <- melt(IMP, id.vars=c("Chemicals"),value.name = "Importance")
    IMP$Importance<-as.numeric(as.character(IMP$Importance))
    IMP$Chemicals<-as.factor(IMP$Chemicals)

    imp<-aggregate(Importance~Chemicals,data=IMP,median)
    imp<-imp[order(imp$Importance),]
    wax<-as.character(imp$Chemicals)
    IMP$Chemicals<-factor(IMP$Chemicals,levels=c(paste(wax,sep=",")))

    pdf(paste("Importance_5foldCV_mtry",mtry,"_ntree",ntree,"_noPFT_allAC_cforest.pdf",sep=""),height=5,width=10,family="sans")
    print(ggplot(IMP,aes(Chemicals,Importance))+geom_boxplot()+
            theme_light()+
            labs(x="Features of maize adult leaf cuticular wax composition and structural features",y="Importance score")+
            #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            theme_light()+
            theme(strip.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  strip.text.y = element_text(size=5,colour = 'black'),
                  strip.text.x = element_text(size=5,colour = 'black'),
                  strip.placement = "outside",
                  axis.text.x = element_text(angle = 90, hjust = 1,vjust=0),
                  text = element_text(size=8),
                  #legend.title=element_blank(),
                  legend.text.align = 0,
                  legend.text=element_text(size=rel(0.8)))

    )
    dev.off()

  }
}

##### generate table ###
#i=47
mtry=10

setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
#IMP_0<-read.table(paste("importanceScore_gc_6fCV_54lines_mtry",i,"_cforest.txt",sep=""),head=T,sep="\t")
IMP_0<-read.table(paste("importanceScore_gc_5fCV_51lines_mtry",mtry,"ntree1000_cforest.txt",sep=""),header=T,sep="\t")
#wax_name<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/short_vs_complete_wax_names.txt",header=T,sep="\t")
wax_name<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/short_vs_complete_wax_names.txt",header=T,sep="\t")

IMP_0<-merge(wax_name,IMP_0,by.x="short_names",by.y="Variable",all.y=T)

IMP<-IMP_0
rownames(IMP)<-IMP[,2]
IMP<-IMP[,-c(1:2)]
stat_mean<-apply(IMP,1,mean)
stat_var<-apply(IMP,1,var)
IMP_sum<-cbind.data.frame(IMP_0[,2],stat_mean,stat_var)
colnames(IMP_sum)<-c("Cuticular features","Mean of importance score","Variance of importance score")

#pearson<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/plot/wax_stomata_ce_cor_naive.txt",header=T,sep="\t")
pearson<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/plot/wax_stomata_ce_cor_naive.txt",header=T,sep="\t")
pearson<-cbind.data.frame(rownames(pearson),pearson[,c(1,2)])
colnames(pearson)<-c("Cuticular features","Pearson's correlation","P-value of correlation")
IMP_Pearson<-merge(pearson,IMP_sum,by="Cuticular features",all=T)
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/cforest")
write.table(IMP_Pearson,"TableSx_pearson_cor_importance_score_0903.txt",row.names=F,sep="\t",quote=F)
###### Ends here ########################
