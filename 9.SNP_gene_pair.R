##############################
# assign SNPs to closest genes
##############################
gff.all<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
gff.all<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
#gff.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_formatted.txt",header=T,sep="\t")
length(unique(gff.all$gene_name))
gff.all$chr<-as.numeric(as.character(gff.all$chr))
gff.all$start<-as.numeric(as.character(gff.all$start))
gff.all$end<-as.numeric(as.character(gff.all$end))

gff.all<-gff.all[order(gff.all$chr,gff.all$start),]

last_bp<-c(307041717,244442276,235667834,246994605,223902240,174033170,182381542,181122637,159769782,150982314) # from .fa.fai file of B73_v4 reference genome, may need to change for other reference genome
## remove reduandent genes
redandent<-c()
for(g in 2:nrow(gff.all)){
  if(gff.all$chr[g]==gff.all$chr[g-1] & gff.all$start[g]==gff.all$start[g-1]& gff.all$end[g]==gff.all$end[g-1] & gff.all$strand[g]==gff.all$strand[g-1]){
    redandent<-c(redandent,g)
  }
}
check<-gff.all[redandent,]
gff.all<-gff.all[-redandent,]

#################################################
# merge "gene regions" for overlapping gene models
#################################################
#gff.all<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/Cand_selection/SNP_gene_pair_test.txt",header=T,sep="\t")

gff.all$start_merge<-gff.all$start
gff.all$end_merge<-gff.all$end

i=2
index1<-2
index<-2
gff.merge<-c()

for (ch in 1:10){
  print(ch)
  gff<-gff.all[which(gff.all$chr==ch),]

  i=2
  while (i <= nrow(gff)){
    if (gff$start[i]<=gff$end[(i-1)]){ ##  overlap

      index1<-i-1 ## the index of the first involved overlapping gene in this merged region
      index<-i ## the index of searching list, will change in each while loop

      ## if the start of gene[index] is smaller than the end of gene[index1]
      ## or if the start of gene[index] is smaller than the end of gene[index-1] (in case there are >= 3 genes)

      #while((gff$end[index1]>=gff$start[index] | gff$start[index]<=gff$end[index-1])&index<=nrow(gff)){ # not use this line
      while(gff$start[index] <= max(gff$end[index1:(index-1)]) &index<=nrow(gff)){ # Xiaowei's edit
        #gff$end_merge[index-1]<-gff$end[index]
        index<-index+1
      }
      gff$start_merge[index1:(index-1)]<-rep(gff$start[index1],(index-index1))
      #gff$end_merge[index1:(index-1)]<-rep(gff$end[(index-1)],(index-index1))
      gff$end_merge[index1:(index-1)] <- rep(max(gff$end[index1:(index-1)]),(index-index1))   ### updated by Xiaowei
      i= index
    }else{
      i=i+1
    }
  }
  gff.merge<-rbind.data.frame(gff.merge,gff)
}

# define extended gene intervals to cover the whole genome
gff.new<-vector()
for (i in 1:10){
  gff.sub<-gff.merge[which(gff.merge$chr==i),]
  gff.sub$ext_end<-last_bp[i]
  nm_gene<-nrow(gff.sub)
  for (j in 1:(nm_gene-1)){ # extended end for the last gene is pre-determined as the last base pair of a chromosome
    gff.sub$ext_end[j]<-(gff.sub$end_merge[j]+gff.sub$start_merge[j+1])/2
  }
  gff.new<-rbind(gff.new,gff.sub)
}
gff.new$gene_name<-as.character(gff.new$gene_name)
#########################
# fix the extended boundary for merged gene regions
#########################
i=nrow(gff.new)-1 ## will check from the end of the list to the start to make sure the correct ext_end for merged genes
while (i>=1){
  if (gff.new$end_merge[i]==gff.new$end_merge[(i+1)]){
    gff.new$ext_end[i]<-gff.new$ext_end[(i+1)]
  }
  i=i-1
}
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/TWAS/")
#write.table(gff.new,"v4.42_gene_pos_overlapMerged_extEnd.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(gff.new,"v4.42_gene_pos_overlapMerged_extEnd_0715.txt",col.names=T,row.names=F,sep="\t",quote=F)



############################################
## start the assignment: closest and 2nd closest
#############################################
chr=1 # 1 to 10

#gff.new<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_overlapMerged_extEnd.txt",header=T,sep="\t")
#gff.new<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_overlapMerged_extEnd.txt",header=T,sep="\t")

gff.new<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_overlapMerged_extEnd_0715.txt",header=T,sep="\t")
gff.new<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_overlapMerged_extEnd_0715.txt",header=T,sep="\t")

### check the number of merged genes due to overlapped gene models
gff.new$interval<-paste(gff.new$start_merge,',',gff.new$ext_end,sep="")
nrow(gff.new)-length(unique(gff.new$interval))
test<-gff.new[duplicated(gff.new$interval),]
test.2<-c()
for(i in 1:nrow(test)){
  interval<-as.character(test$interval[i])
  temp<-gff.new[which(gff.new$interval==interval),]
  test.2<-rbind.data.frame(test.2,temp)
}

nrow(test.2)-(nrow(gff.new)-length(unique(gff.new$interval)))
## test
#snps<-read.table(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection/test_SNP_pos.txt",sep=""),header=T)


#for (chr in 1:10){
  print(chr)
  snps<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/filter_prune_468K/imputedGBS_filtered_chr",chr,"_pos_ID.txt",sep=""),header=F)
  #snps<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/filter_prune_468K/imputedGBS_filtered_chr",i,"_pos_ID.txt",sep=""),header=F)
  colnames(snps)<-c("chr","pos","ID")

  SNP_gene1<-c()

  gff.sub<-gff.new[which(gff.new$chr==chr),]

  j=1
  #while (j <= 100){
  while (j <= nrow(gff.sub)){
    if(j %% 100==0) {
      # Print on the screen some message
      cat(paste0("iteration: ", j, "\n"))
    }
    #start_merge<-gff.sub$start_merge[j]
    end_merge<-gff.sub$end_merge[j]
    index<-which(gff.sub$end_merge==end_merge)
    gene_cluster<-gff.sub[index,]
    if(j>1){
      #snp_group<-snps[which(snps$pos>gff.sub$ext_end[j-1] & snps$pos<=gff.sub$ext_end[j]),]
      snp_group<-snps[which(snps$pos>=gff.sub$ext_end[j-1] & snps$pos<=gff.sub$ext_end[j]),] # Xiaowei's edit
    }else if (j==1){
      snp_group<-snps[which(snps$pos<=gff.sub$ext_end[j]),]
    }
    snp_group_multi<-snp_group[rep(seq_len(nrow(snp_group)), each=nrow(gene_cluster)), ]
    gene_cluster_multi<-gene_cluster[rep(seq_len(nrow(gene_cluster)), nrow(snp_group)), ]
    snp_gene1<-cbind.data.frame(snp_group_multi,gene_cluster_multi)
    SNP_gene1<-rbind.data.frame(SNP_gene1,snp_gene1)
    j<-index[length(index)]+1
  }
  write.table(SNP_gene1,paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_SNP-Gene_merged1_0715.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
#}

### second closest gene:
gff.new<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_overlapMerged_extEnd.txt",header=T,sep="\t")
gff.new<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.42_gene_pos_overlapMerged_extEnd.txt",header=T,sep="\t")
gff.new$gene_name<-as.character(gff.new$gene_name)

chr=1
for(chr in 1:10){
  SNP_gene1<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_SNP-Gene_merged1.txt",sep=""),header=T,sep="\t")
  SNP_gene1$ID<-as.character(SNP_gene1$ID)
  SNP_gene1$gene_name<-as.character(SNP_gene1$gene_name)

  SNP_gene2<-c()

  gff.sub<-gff.new[which(gff.new$chr==chr),]

  j=1
  #while (j <= 100){
  while (j <= nrow(gff.sub)){
    if(j %% 100==0) {
      # Print on the screen some message
      cat(paste0("iteration: ", j, "\n"))
    }

    end_merge<-gff.sub$end_merge[j]
    index<-which(gff.sub$end_merge==end_merge) # index of all genes in a gene cluster
    gene_cluster<-gff.sub[index,]
    gene1<-gene_cluster$gene_name[1]# first gene in a gene cluster

    # pull out SNPs assigned to this gene cluster (using the first gene ID in a cluster)
    snp_group<-SNP_gene1[which(SNP_gene1$gene_name==gene1),]

    gff.sub.update<-gff.sub[-index,] # remove a gene cluster from gff.sub

    # update "extended gene region" for gff.sub.update
    if(nrow(snp_group)>0){
      if (index[1]!=1 & index[length(index)]!=nrow(gff.sub)){
        gff.sub.update$ext_end[(index[1]-1)]<-(gff.sub.update$end_merge[(index[1]-1)]+gff.sub.update$start_merge[index[1]])/2

        ## fix the ext_end for "merged gene regions"
        gff.sub.update$ext_end[which(gff.sub.update$end_merge==gff.sub.update$end_merge[(index[1]-1)])]<-gff.sub.update$ext_end[(index[1]-1)]

        ## pull out gene clusters right before and after the removal
        index_bf<-which(gff.sub.update$end_merge==gff.sub.update$end_merge[(index[1]-1)])
        index_af<-which(gff.sub.update$end_merge==gff.sub.update$end_merge[index[1]])
        gene_cluster_bf<-gff.sub.update[index_bf,]
        gene_cluster_af<-gff.sub.update[index_af,]

        ## re-assign SNPs in snp_group to adjacent gene clusters
        snp_group_bf<-snp_group[which(snp_group$pos<=gff.sub.update$ext_end[(index[1]-1)]),1:3]
        snp_group_af<-snp_group[which(snp_group$pos>gff.sub.update$ext_end[(index[1]-1)]),1:3]

        ## check the amount
        if(nrow(snp_group_bf)+nrow(snp_group_af)!=nrow(snp_group)){
          break
        }

        ## format and aggregate
        snp_group_multi_bf<-snp_group_bf[rep(seq_len(nrow(snp_group_bf)), each=nrow(gene_cluster_bf)), ]
        gene_cluster_multi_bf<-gene_cluster_bf[rep(seq_len(nrow(gene_cluster_bf)), nrow(snp_group_bf)), ]
        snp_group_multi_af<-snp_group_af[rep(seq_len(nrow(snp_group_af)), each=nrow(gene_cluster_af)), ]
        gene_cluster_multi_af<-gene_cluster_af[rep(seq_len(nrow(gene_cluster_af)), nrow(snp_group_af)), ]

        snp_gene2_bf<-cbind.data.frame(snp_group_multi_bf,gene_cluster_multi_bf)
        snp_gene2_af<-cbind.data.frame(snp_group_multi_af,gene_cluster_multi_af)

        SNP_gene2<-rbind.data.frame(SNP_gene2,snp_gene2_bf,snp_gene2_af)
        j<-index[length(index)]+1

      }else if(index[1]==1){
        ## pull out gene clusters right before and after the removal
        index_af<-which(gff.sub.update$end_merge==gff.sub.update$end_merge[1])
        gene_cluster_af<-gff.sub.update[index_af,]

        snp_group_af<-snp_group[,1:3]## all SNPs used in the first gene cluster (on a chr) will be assigned to the second

        snp_group_multi_af<-snp_group_af[rep(seq_len(nrow(snp_group_af)), each=nrow(gene_cluster_af)), ]
        gene_cluster_multi_af<-gene_cluster_af[rep(seq_len(nrow(gene_cluster_af)), nrow(snp_group_af)), ]

        snp_gene2_af<-cbind.data.frame(snp_group_multi_af,gene_cluster_multi_af)

        SNP_gene2<-rbind.data.frame(SNP_gene2,snp_gene2_af)
        j<-index[length(index)]+1


      }else if(index[length(index)]==nrow(gff.sub)){ ## if the closest genes include the last gene on the chr
        #gff.sub.update$ext_end[(index[1]-1)]<-last_bp[chr]

        ## fix the ext_end for "merged gene regions"
        #gff.sub.update$ext_end[which(gff.sub.update$end_merge==gff.sub.update$end_merge[(index[1]-1)])]<-gff.sub.update$ext_end[(index[1]-1)]

        ## pull out gene clusters right before the removal
        index_bf<-which(gff.sub.update$end_merge==gff.sub.update$end_merge[index[1]-1])
        gene_cluster_bf<-gff.sub.update[index_bf,]

        snp_group_bf<-snp_group[,1:3]## all SNPs used in the last gene cluster (on a chr) will be assigned to the last but one

        snp_group_multi_bf<-snp_group_bf[rep(seq_len(nrow(snp_group_bf)), each=nrow(gene_cluster_bf)), ]
        gene_cluster_multi_bf<-gene_cluster_bf[rep(seq_len(nrow(gene_cluster_bf)), nrow(snp_group_bf)), ]

        snp_gene2_bf<-cbind.data.frame(snp_group_multi_bf,gene_cluster_multi_bf)

        SNP_gene2<-rbind.data.frame(SNP_gene2,snp_gene2_bf)
        j<-index[length(index)]+1
      }# end of if.. else if .. else if..
    } else {
      j<-index[length(index)]+1
    } # end of if nrow(snp_group)>0
  }# end of while for gene cluster j

  write.table(SNP_gene2,paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_SNP-Gene_merged2.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

}
