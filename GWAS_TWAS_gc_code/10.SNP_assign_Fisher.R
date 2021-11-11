# ### SNP-gene pairs
#
# chr=1
# ###### just need to run once #######################
#
#
# #gff<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
# #gff<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
# gff<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
# gff$center<-(gff$start+gff$end)/2
# gff$gene_name<-as.character(gff$gene_name)
#
# nearestGenes<-function(snp,gff.sub,nm_gene=5){
#   ch<-snp[1]
#   pos<-snp[2]
#   #snpID<-snp[3]
#   gff.sub$distTOcenter<-abs(gff.sub$center-pos)
#   gff.sub<-gff.sub[order(gff.sub$distTOcenter),]
#   near_genes<-gff.sub[1:nm_gene,]
#   #snps<-snp[rep(seq_len(nrow(snp)),each=nm_gene),]
#   snps<-matrix(rep(snp,nm_gene),nrow=nm_gene,byrow=T)
#   near_genes<-cbind(snps,near_genes)
#   return(near_genes)
# }
#
# #pairs<-nearestGenes(snp,gff.sub,nm_gene=5)
# #chr=1
#
# #for (chr in 2:10){
#   SNP<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/filter_prune_468K/imputedGBS_filtered_chr",chr,"_pos_ID.txt",sep=""),header=F)
#   #SNP<-read.table(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/filter_prune_468K/imputedGBS_filtered_chr",chr,"_pos_ID.txt",sep=""),header=F)
#   colnames(SNP)<-c("chr","pos","ID")
#
#   SNP.1<-SNP[,1:2]
#   SNPid<-cbind.data.frame(SNP[,3],paste(SNP[,1],"_",SNP[,2],sep=""))
#   colnames(SNPid)<-c("ID","chr_pos")
#   SNPid$ID<-as.character(SNPid$ID)
#   SNPid$chr_pos<-as.character(SNPid$chr_pos)
#
#   gff.sub<-gff[which(gff$chr==chr),c(1,7,6)]
#   # for (i in 1:nrow(SNP)){
#   #   snp<-SNP[i,]
#   #   ch<-snp[1,1]
#   #   pos<-snp[1,2]
#   #   snpID<-snp[1,3]
#   #   gff.sub$distTOcenter<-abs(gff.sub$center-pos)
#   #   gff.sub<-gff.sub[order(gff.sub$distTOcenter),]
#   #   near_genes<-gff.sub[1:nm_gene,]
#   #   snps<-snp[rep(seq_len(nrow(snp)),each=nm_gene),]
#   #   near_genes<-cbind(snps,near_genes)
#   # }
#
#   # Pair<-vector()
#   # for(j in 1:2){
#   #   snp=SNP.1[j,]
#   #   pair<-nearestGenes(snp, gff.sub,nm_gene=5)
#   #   print(pair)
#   #   Pair<-rbind(Pair,pair)
#   # }
#   # SNP.test<-SNP.1[1:2,]
#
#   pairs<-apply(SNP.1,1,FUN=function(snp) nearestGenes(snp, gff.sub,nm_gene=5))
#
#   PAIR<-do.call(rbind, pairs)
#
#   PAIR$chr_pos<-paste(PAIR$`1`,"_",PAIR$`2`,sep="")
#   PAIR.final<-merge(SNPid,PAIR,by="chr_pos",all=T)
#   colnames(PAIR.final)[2:4]<-c("SNP","SNP_chr","pos")
#  write.table(PAIR.final,paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_SNPto5Genes.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
#}



#########################
## gwas results
rep<-"BLUP"

#for (trait in c("ce","ft")){
for (trait in c("ce")){
  ## gwas results
  gwas<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/",trait,"_SD18_both_HapMap3_topSNPs01.txt",sep=""),header=T,sep="\t")
  #gwas<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/",trait,"_SD18_both_HapMap3_topSNPs01.txt",sep=""),header=T,sep="\t")

  #for (pairing_set in c("first","second")){
  for (pairing_set in c("first")){
    gwas$SNP<-as.character(gwas$SNP)

    ### combine gwas results with SNP-gene pairs
    gwas_new<-vector()

    for (chr in 1:10){
      gwas_sub<-gwas[which(gwas$Chromosome==chr),]

      ## overlapping gene considered
      if(pairing_set=="first"){
        pair<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_SNP-Gene_merged1_0715.txt",sep=""),header=T,sep="\t")
      }else if(pairing_set=="second"){
        pair<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_SNP-Gene_merged2.txt",sep=""),header=T,sep="\t")
      }

      #pair<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_geneAssignment.txt",sep=""),header=T,sep="\t")
      #pair<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_geneAssignment.txt",sep=""),header=T,sep="\t")
      #pair.1<-melt(pair,id.vars=c("ID"),measure.vars = c("gene", "gene2"))

      gwas_sub<-merge(gwas_sub,pair,by.x="SNP",by.y="ID",all.x=T)
      gwas_new<-rbind.data.frame(gwas_new,gwas_sub)
    }

    gwas_new_1<-gwas_new[,c(1:5,16,18,19)]

    #setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb/")
    setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb/")

    if(pairing_set=="first"){
      write.table(gwas_new_1,paste("Gene_SNP_pairs_forFisher_",trait,"_HMP3_kinLD02_450to310_merged1_0715.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
    }else if(pairing_set=="second"){
      write.table(gwas_new_1,paste("Gene_SNP_pairs_forFisher_",trait,"_HMP3_kinLD02_450to310_merged2.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
    }

  }
}


## the above code just need run once

#################################################
# combine with TWAS results
#################################################
pairing_set<-1 # 1 or 2

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb/")


# 300 kinship,BLUP (OLRM)
rep<-"BLUP"
#traits<-c("ce","ft")
Models<-c("Kin")
#model="Kin" # "PC0","PC5","PC10", "Kin"
K=20

#twas.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D/rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS.txt",header=T)
#twas.all<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D/rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS.txt",header=T)

## TWAS results including plastid and contig genes
twas.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D/rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS_MtCp.txt",header=T)

## convert -logp to pvalue
for(i in 4:ncol(twas.all)){
  twas.all[,i]<-10^((-1)*twas.all[,i])
}


Fisher <- function(x){res <- sumlog(x); return(res$p)}

for (k in 4:ncol(twas.all)){
  twas<-twas.all[,c(1:3,k)]
  t=colnames(twas)[4]
  t=strsplit(t,"_",fixed=T)
  t=t[[1]][1] # trait
  colnames(twas)[4]<-"pvalue"

  #GWAS_new<-read.table(paste("Gene_SNP_pairs_forFisher_",t,"_HMP3_kinLD02_450to310_merged",pairing_set,".txt",sep=""),header=T,sep="\t")
  GWAS_new<-read.table(paste("Gene_SNP_pairs_forFisher_",t,"_HMP3_kinLD02_450to310_merged",pairing_set,"_0715.txt",sep=""),header=T,sep="\t")

  for (model in Models){

    #T_GWAS<-merge(twas,GWAS_new,by.x="gene_name",by.y="value",all.y=T) ## 1 SNP paired with 2 genes
    T_GWAS<-merge(twas,GWAS_new,by="gene_name",all.y=T) # SNP paired with merged gene region
    T_GWAS<-as.data.frame(T_GWAS)
    T_GWAS$pvalue[is.na(T_GWAS$pvalue)]<-1

    ###############################################################
    # Fisher's combined test
    ###############################################################
    #install.packages("metap")
    library(metap)

    df <- cbind(T_GWAS$P.value,T_GWAS$pvalue)
    a=apply(df, 1, Fisher)

    fisher<-cbind.data.frame(T_GWAS$gene_name,T_GWAS$SNP,T_GWAS$Chromosome,T_GWAS$Position,a,T_GWAS$P.value,T_GWAS$pvalue,T_GWAS$strand)
    colnames(fisher)<-c("gene_name","SNP","chr","position","p","pG","pT","strand")
    fisher$DF<-4

    #write.table(fisher,paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_merge",pairing_set,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    #write.table(fisher,paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_MtCp_merge",pairing_set,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(fisher,paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_MtCp_merge",pairing_set,"_0715.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  } # end of model loop
}# end of trait loop

##########################################
# how many SNPs in a gene; how many genes and SNPs were tested
##################################################
library(data.table)

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb")

fisher_res<-fread("FisherRes_ce_Kin_TWAS.BLUP_K20_HMP3.450to310_rlog_MtCp_merge1.txt",data.table=F)
fisher_res$gene_name<-as.character(fisher_res$gene_name)
fisher_res$SNP<-as.character(fisher_res$SNP)

nm_SNP<-length(unique(fisher_res$SNP)) #971508
nm_gene<-length(unique(fisher_res$gene_name)) #31879, 45041 genes in the genome

uniq_gene<-unique(fisher_res$gene_name)
nm_snp<-c()
for (g in uniq_gene){
  sub_fisher<-fisher_res[which(fisher_res$gene_name==g),]
  nm<-nrow(sub_fisher)
  nm_snp<-c(nm_snp,nm)
}

range(nm_snp) #1 to 2271
mean(nm_snp) # 32.55297
median(nm_snp) # 10

nm_snp_inGene<-cbind.data.frame(uniq_gene,nm_snp)
write.table(nm_snp_inGene,"nm_snp_inGene_Fisher_test.txt",row.names=F,sep="\t",quote=F)
pdf("nm_snp_inGene_Fisher_test.pdf",height=4,width = 5,family="serif")
hist(nm_snp,xlab="Number of SNPs paired with a gene")
dev.off()

twas.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Kin_LD02/P3D/rlog_PEER_K20_noP3D_kin_ce_ft_SD18_TWAS_MtCp.txt",header=T)
twas.all$gene_name<-as.character(twas.all$gene_name)
no_twas<-fisher_res[-which(fisher_res$gene_name %in% twas.all$gene_name),]
length(unique(no_twas$gene_name)) #16019

no_twas2<-fisher_res[which(fisher_res$pT==1),]
length(unique(no_twas2$gene_name))

31879-16019

#################################################
################################################
# Manhattan plot
################################################
library(qqman)
library(qvalue)
align<-"hisat"
#act_P=1
K=20
rep="BLUP"
transf="rlog"


# 300 kinship,BLUP (OLRM)
traits<-c("ce","ft")
#Models<-c("PC0","PC5","PC10", "Kin")
Models<-"Kin"
#model<-"PC0" # "PC0","PC1", "Kin"

pairing_set<-1

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb")

#filename<-paste("FisherRes_",model,"_PEERres25_TWAS.",rep,"_HMP3.310IMP.prun09",sep="")
#filename<-paste(m,"_",t,"_",transf,"_",align,"GE_PC",act_P,"_PEER",K,"_TWAS_results_rep",rep,sep="")

for (t in traits[1]){
  for(model in Models){

    #filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_merge",pairing_set,sep="")
    #filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_MtCp_merge",pairing_set,sep="")
    filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_MtCp_merge",pairing_set,"_0715",sep="")
    twas<-read.table(paste(filename,".txt",sep=""),header=T,sep="\t")

    pdf(paste(filename,".pdf",sep=""),width=8,height=4)
    manhattan(twas, chr = "chr", bp = "position", p = "p", snp = "SNP",
              col = c("navy", "darkorange1"), chrlabs = NULL,
              suggestiveline = FALSE,genomewideline =FALSE,
              main=paste(filename))
    dev.off()

    jpeg(paste(filename,"_qq.jpg",sep=""))
    qq(twas$p)
    dev.off()

  }
}

######################
#####################################
# pull out important genes
#####################################
library(pkgcond)

align<-"hisat"
#act_P=2
K=20
rep="BLUP" # 1, 2,"b","PEERresBLUE"
transf="rlog"
traits<-c("ce","ft")

#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb")
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb")

#annotation<-read.delim("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
#annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene_name"
Models<-"Kin"

for (t in traits[1]){
  # if (t=="ce"){
  #   Models=c("PC5", "Kin")
  # }else if(t=="ft"){
  #   Models=c("PC10", "Kin")
  # }

  for (model in Models){
    #setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb")
    setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/Fishers_comb")
    #filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_HMP3.450to310",sep="")
    #filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_merge",pairing_set,sep="")
    #filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_MtCp_merge",pairing_set,sep="")
    filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_K",K,"_HMP3.450to310_rlog_MtCp_merge",pairing_set,"_0715",sep="")

    twas0<-read.table(paste(filename,".txt",sep=""),header=T,sep="\t")
    twas0<-twas0[order(twas0$p),]
    twas0$gene_name<-as.character(twas0$gene_name)
    nm_cand<-ceiling(length(unique(twas0$gene_name))*0.01)

    # F_genes<-twas0[1,]
    # r=2
    # while(nrow(F_genes)<nm_cand){
    #   check<-twas0[r,]
    #   if(check$gene_name %!in% F_genes$gene_name){
    #     F_genes<-rbind.data.frame(F_genes,twas0[r,])
    #   }
    #   r=r+1
    # }

    uniq_genes<-unique(twas0$gene_name)

    # for checking
    #F_genes.2<-uniq_genes[1:nm_cand]
    #which(F_genes$gene_name != F_genes.2)


    F_genes<-twas0[match(uniq_genes[1:nm_cand],twas0$gene_name),]
    # gene_list<-c("Zm00001d029140","Zm00001d025149","Zm00001d043509","Zm00001d044162","Zm00001d051599",
    #              "Zm00001d036765","Zm00001d010426")
    # F_genes[which(F_genes$gene_name %in% gene_list),]

    #which(F_genes$gene_name != F_genes.2$gene_name)
    #which(F_genes$p != F_genes.2$p)
    result_sum<-merge(F_genes,annotation,by="gene_name",all.x=T)

      #setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection")
      #setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection")
      setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection")
      #write.table(F_genes,paste("Cand_",filename,"_",postfix,".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
      write.table(result_sum,paste("Cand_",filename,"_anno.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
      write.table(F_genes,paste("Cand_",filename,".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
  }
}


##### *****should give exactly same results with TWAS that not includes plastid and contig genes

#### check overlapping gene models

## overlap_geneID is from "SNP_gene_pair_0415.R"
overlap_geneID<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/overlapping_gene_list.txt",header=T,sep="\t")
overlap_geneID$gene_name<-as.character(overlap_geneID$gene_name)

fisher_cand<-read.delim("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/Cand_FisherRes_ce_Kin_TWAS.BLUP_K20_HMP3.450to310_rlog_merge1_anno.txt",header=T,sep="\t")
fisher_cand$gene_name<-as.character(fisher_cand$gene_name)
uniq_cand<-unique(fisher_cand$gene_name)
#uniq_fisher<-fisher_cand[match(uniq_cand,fisher_cand$gene_name),]

flag_genes<-c()
for (i in nrow(overlap_geneID)/2){
  ID1<-overlap_geneID$gene_name[(2*i-1)]
  ID2<-overlap_geneID$gene_name[(2*i)]

  if(ID1 %in% uniq_cand & ID2 %in% uniq_cand){
    flag_genes<-c(flag_genes,ID1,ID2)
  }
}
length(flag_genes)
## result is null. Good!
