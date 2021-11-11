#thres=0.0002

setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_3_005")## Up-to-date verstion for 0.05% SNP
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_3_01")## Up-to-date verstion for 0.1% SNP
library(data.table)


# set variables -----------------------------------------------------------
break_width=200000 # used in step 1, distance of adjacent SNPs to break peaks
win_size=10 #used in step 2
step_size=5 #used in step 2

#import gwas results
gwas_res=fread('/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_2/ce_SD18_both_HapMap3_topSNPs_forCand_test.txt',data.table = F, stringsAsFactors = F)
gwas_res=gwas_res[,c("SNP","Chromosome","Position","P.value","qvalues","trait")]
#change to numeric
gwas_res<-gwas_res[order(gwas_res$P.value),]
#gwas_res<-gwas_res[1:ceiling(nrow(gwas_res)/2),] ## test 0.05% GWAS SNPs
gwas_res$Chromosome=as.numeric(gwas_res$Chromosome)
gwas_res$Position=as.numeric(gwas_res$Position)
gwas_res$qvalues=as.numeric(gwas_res$qvalues)
gwas_res$P.value=as.numeric(gwas_res$P.value)

gwas_res$logp=-log10(gwas_res$P.value)
traits=unique(gwas_res$trait)

#rm(curr_gwas_peak1.1)

# loop through each traits, assign peak numbering 1.1 ---------------------

for(t in traits){

  print(t)
  #gwas results of current trait
  curr_gwas=gwas_res[which(gwas_res$trait==t),]

  peak_count=0
  uni_chr=sort(unique(curr_gwas$Chromosome))

  #loop through each chromosome
  for (chr in uni_chr){
    # chr <- 1
    print(chr)
    peak_count=peak_count+1
    curr_chr=curr_gwas[which(curr_gwas$Chromosome==chr),]
    #sort by physical position
    curr_chr_sorted=curr_chr[order(curr_chr$Position),]
    #add additional columns to calculate physical position distance between two adjacent SNPs, and assign peak number
    curr_chr_sorted$dis[1]=0
    curr_chr_sorted$peak_num1.1[1]=peak_count
    if (nrow(curr_chr_sorted)>1){#only do this if more than 1 SNP in the peak
      for(i in 2:nrow(curr_chr_sorted)){
        curr_chr_sorted$dis[i]=curr_chr_sorted$Position[i]-curr_chr_sorted$Position[i-1]
        if(curr_chr_sorted$dis[i]<break_width){
          curr_chr_sorted$peak_num1.1[i]=peak_count
        }else{
          peak_count=peak_count+1
          curr_chr_sorted$peak_num1.1[i]=peak_count
        }
      }
      #append to a file containing all chr
      if('curr_gwas_peak1.1' %in% ls()){
        curr_chr_sorted$peak_num1.1.1=curr_chr_sorted$peak_num1.1
        curr_gwas_peak1.1$peak_num1.1.1=curr_gwas_peak1.1$peak_num1.1
        curr_gwas_peak1.1=rbind.data.frame(curr_gwas_peak1.1,curr_chr_sorted)
      }else{
        curr_gwas_peak1.1=curr_chr_sorted
      }
    }else{ ## Meng's edit
      #append to a file containing all chr
      if('curr_gwas_peak1.1' %in% ls()){
        curr_chr_sorted$peak_num1.1.1=curr_chr_sorted$peak_num1.1
        curr_gwas_peak1.1$peak_num1.1.1=curr_gwas_peak1.1$peak_num1.1
        curr_gwas_peak1.1=rbind.data.frame(curr_gwas_peak1.1,curr_chr_sorted)
      }else{
        curr_gwas_peak1.1=curr_chr_sorted
      }
    } ## end of Meng's edit
  }   ### end loop chr


  if ('curr_gwas_peak1.1' %in% ls()){
    #rearrange order to avoid any skipped peak numbers
    c=1
    curr_gwas_peak1.1$peak_num1.1.1=curr_gwas_peak1.1$peak_num1.1
    for (x in unique(curr_gwas_peak1.1$peak_num1.1)){
      ind=which(curr_gwas_peak1.1$peak_num1.1==x)
      curr_gwas_peak1.1$peak_num1.1.1[ind]=rep(c,length(ind))
      c=c+1
    }
    curr_gwas_peak1.1$peak_num1.1=curr_gwas_peak1.1$peak_num1.1.1

    #write out v1.1
    write.csv(curr_gwas_peak1.1[,c(1:9)],paste('gwas_peak_v1.1_',t,'.csv',sep=''),row.names = F,quote = F)
    rm('curr_gwas_peak1.1')
  }
}

# loop through each traits, assign peak numbering 1.2 (sliding window) --------

for(t in traits){

  print(t)
  #peak assignment 1.1 of current trait
  curr_gwas=fread(paste('gwas_peak_v1.1_',t,'.csv',sep=''),data.table = F)

  peak_count=max(curr_gwas$peak_num1.1)+1
  uni_peak=max(curr_gwas$peak_num1.1)

  #loop through each peak
  # pdf(paste('peak_v1.1_',t,'.pdf',sep=''))

  for (i in 1:uni_peak){
    print(i)
    # i <- 1
    curr_peak=curr_gwas[which(curr_gwas$peak_num1.1==i),]
    curr_peak$peak_num1.2= curr_peak$peak_num1.1

    if (nrow(curr_peak)<(win_size+step_size)){ #if peak has less than 10 SNPs, no further splitting
    }else if (nrow(curr_peak)>=(win_size+step_size)){#if peak has more than 10 SNPs, see if there are needs for splitting
      curr_peak$order=1:nrow(curr_peak)
      num_win=ceiling(nrow(curr_peak)/step_size)
      win_mean=matrix(nrow=0,ncol=2)

      for(x in 1:num_win){
        # x <- 1
        curr_SNP=(step_size*(x-1)+1) : (step_size*(x-1)+win_size)
        curr_logp=round(median(curr_peak$logp[which(curr_peak$order %in% curr_SNP)]))   ###??? median
        win_mean=rbind(win_mean,c(x,curr_logp))
      }

      win_mean=data.frame(win_mean)
      names(win_mean)=c('win_num','mean_logp')
      win_mean$diff[1]=0

      for(x in 2:nrow(win_mean)){
        win_mean$diff[x]=win_mean$mean_logp[x]-win_mean$mean_logp[x-1]
      }

      #find all positive and negative values
      ind_neg=which(win_mean$diff<0)
      ind_pos=which(win_mean$diff>0)
      ind=c()

      if (length(ind_neg)>1){ #if more than one negative values
        for (z in 1:(length(ind_neg)-1)){
          # z <- 1
          #find if there are positive values between two negatives   ###???# the crucial step
          ind_curr=ind_pos[which(ind_pos>ind_neg[z])][1]

          if (is.na(ind_curr)){
          }else if (ind_curr<ind_neg[z+1]){
            ind=c(ind,ind_curr)
          }
        }
        #manually check the last one
        ind_curr=ind_pos[which(ind_pos>ind_neg[length(ind_neg)])][1]
        ind=c(ind,ind_curr)
      }else{
        ind_curr=ind_pos[which(ind_pos>ind_neg)][1]
        ind=c(ind,ind_curr)
      }

      #assign additional peak number
      if (length(ind)>0){
        for (y in ind){
          # y <- 8
          #mid point of the window
          break_point=round((step_size*(y-1)+1)+win_size/2)
          curr_peak$peak_num1.2[which(curr_peak$order>break_point)]=peak_count
          peak_count=peak_count+1
        }#end loop y
      }#end if

      #plotting each peaks for visual check
      # plot(curr_peak$logp~curr_peak$Position,main=paste('peak1.1:',i),pch=16)
      #
      # if (length(ind)>0){
      #   ps <- unique(curr_peak$peak_num1.2)
      #   for (i in ps){
      #     # i <- ps[1]
      #     curr=curr_peak[which(curr_peak$peak_num1.2==i),]
      #     s=curr[which(curr$P.value==min(curr$P.value))[1],]
      #     abline(v=s$Position)
      #   }
      # }
    }#end if

    #append to one file
    if(i==1){
      curr_gwas_peak1.2=curr_peak[,1:10]
    }else{
      curr_gwas_peak1.2=rbind.data.frame(curr_gwas_peak1.2,curr_peak[,1:10])
    }

  }#end loop i

  #rearrange order to avoid any skipped peak numbers
  c=1
  curr_gwas_peak1.2$peak_num1.2s=curr_gwas_peak1.2$peak_num1.2
  for (x in unique(curr_gwas_peak1.2$peak_num1.2)){
    ind=which(curr_gwas_peak1.2$peak_num1.2==x)
    curr_gwas_peak1.2$peak_num1.2s[ind]=rep(c,length(ind))
    c=c+1
  }
  curr_gwas_peak1.2$peak_num1.2=curr_gwas_peak1.2$peak_num1.2s
  curr_gwas_peak1.2 <- curr_gwas_peak1.2[,-11]

  #write out v1.2
  write.csv(curr_gwas_peak1.2,paste('gwas_peak_v1.2_sliding_window_',t,'.csv',sep=''),row.names = F,quote = F)
  # dev.off()

  #assign peak SNPs
  peaks=data.frame(matrix(nrow=max(curr_gwas_peak1.2$peak_num1.2),ncol=ncol(curr_gwas_peak1.2)+1))   ###???# adding +1

  for (i in 1:max(curr_gwas_peak1.2$peak_num1.2)){
    # i <- 63
    curr=curr_gwas_peak1.2[which(curr_gwas_peak1.2$peak_num1.2==i),]
    s=curr[which(curr$P.value==min(curr$P.value))[1],]
    s$SNP_count=nrow(curr)
    peaks[i,]=s
  }

  names(peaks)=c(names(curr_gwas_peak1.2),"SNP_count")
  print(nrow(peaks))
  #write out v1.2
  write.table(peaks,paste('peak_snps_v1.2_',t,'.csv',sep=''),sep=',',quote=F,row.names = F)
}

# using LD to merge SNP, assigning peak loci v2 ---------------------------
#setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_2")
#setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/test_Di_code_3")
library(data.table)

system('mkdir SNPs_to_extract_geno')
traits=c("ce")
for(t in traits){

  print(t)
  curr_trait_peaks=fread(paste('peak_snps_v1.2_',t,'.csv',sep=''),data.table = F)

  #sort by chromosome and position
  curr_trait_peaks=curr_trait_peaks[order(curr_trait_peaks$Chromosome,curr_trait_peaks$Position),]

  #loop through each chromosome
  for (chr in unique(curr_trait_peaks$Chromosome)){
    # chr <- 1
    curr_chr=curr_trait_peaks[which(curr_trait_peaks$Chromosome==chr),]

    if (nrow(curr_chr)>=2){
      #calculate distance between peak SNPs
      curr_chr$dis[1]=0
      for(i in 2:nrow(curr_chr)){
        curr_chr$dis[i]=curr_chr$Position[i]-curr_chr$Position[i-1]
      }

      #extract genotype of all peak snps
      write.table(curr_chr[,2:3],paste('SNPs_to_extract_geno/',t,'_chr',chr,'.txt',sep=''),sep='\t',quote=F,row.names = F)
      # system(paste('vcftools --gzvcf /workdir/dw524/ionomics/Hapmap_genotype/post_imputation_filtering/3.filtering/AGPv4_Ames_ionomics',chr,'_maf.vcf.gz --positions SNPs_to_extract_geno/',t,'_chr',chr,'.txt --recode --out SNPs_to_extract_geno/peak_snps_trait_',t,'_chr_',chr,' &',sep=''))
      system(paste('vcftools --gzvcf /workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/filter_prune_468K/MLC_GBSSNP310_468K_v4_chr',chr,'_filter.recode.vcf.gz --positions SNPs_to_extract_geno/',t,'_chr',chr,'.txt --recode --out SNPs_to_extract_geno/peak_snps_trait_',t,'_chr_',chr,' &',sep=''))
    }
  }
}

system('wait')

#loop again to calculate all pairwise LD
for(t in traits){
  print(t)
  curr_trait_peaks=fread(paste('peak_snps_v1.2_',t,'.csv',sep=''),data.table = F)

  #sort by chromsome and position
  curr_trait_peaks=curr_trait_peaks[order(curr_trait_peaks$Chromosome,curr_trait_peaks$Position),]

  #loop through each chromosome
  for (chr in unique(curr_trait_peaks$Chromosome)){
    curr_chr=curr_trait_peaks[which(curr_trait_peaks$Chromosome==chr),]
    if (nrow(curr_chr)>=2){
      #calculate distance between peak SNPs
      curr_chr$dis[1]=0
      for(i in 2:nrow(curr_chr)){
        curr_chr$dis[i]=curr_chr$Position[i]-curr_chr$Position[i-1]
      }

      #extract genotype of all peak snps
      write.table(curr_chr[,2:3],paste('SNPs_to_extract_geno/',t,'_chr',chr,'.txt',sep=''),sep='\t',quote=F,row.names = F)
      system(paste('/programs/tassel-5-standalone/run_pipeline.pl -Xmx350g -vcf SNPs_to_extract_geno/peak_snps_trait_',t,'_chr_',chr,'.recode.vcf -ld -ldType All -export SNPs_to_extract_geno/',t,'_chr_',chr,'_all_pairwise_LD.txt &',sep=''))
    }
  }
}


## getting a list of peak SNP pairs that have LD > 0.2, and filter
library(data.table)
for(t in traits){
  print(t)
  curr_trait_snps=fread(paste('gwas_peak_v1.2_sliding_window_',t,'.csv',sep=''),data.table = F)
  curr_trait_peaks=fread(paste('peak_snps_v1.2_',t,'.csv',sep=''),data.table = F)

  #sort by chromosome and position
  curr_trait_peaks=curr_trait_peaks[order(curr_trait_peaks$Chromosome,curr_trait_peaks$Position),]
  curr_trait_peaks$id=paste(curr_trait_peaks$Chromosome,curr_trait_peaks$Position,sep='_')

  #get chromosomes we calculated LD on
  all_chr=list.files(path = "SNPs_to_extract_geno/", pattern = '_all_pairwise_LD.txt')
  all_chr=all_chr[grep(paste(t,'_',sep=''),all_chr)]
  all_chr=gsub('_all_pairwise_LD.txt','',all_chr)
  all_chr=as.numeric(gsub(paste(t,'_chr_',sep=''),'',all_chr))

  #loop through each chromosome
  peak_snp_remove_list=c()

  for (chr in all_chr){
    # chr <- 1
    LD=read.delim(paste('SNPs_to_extract_geno/',t,'_chr_',chr,'_all_pairwise_LD.txt',sep=''))
    LD_s=LD[,c('Locus1','Position1','Position2','R.2')]
    LD_0.2=LD_s[which(LD_s$R.2>0.2),]
    LD_0.2$snp1=paste(LD_0.2$Locus1,LD_0.2$Position1,sep='_')
    LD_0.2$snp2=paste(LD_0.2$Locus1,LD_0.2$Position2,sep='_')
    write.csv(LD_0.2,paste('SNPs_to_extract_geno/',t,'_chr_',chr,'_LD0.2.csv',sep=''),row.names = F,quote=F)
    if(nrow(LD_0.2)>0){
      for (i in 1:nrow(LD_0.2)){ #remove less significant SNP per pair
        # i <- 1
        p1=curr_trait_peaks$logp[which(curr_trait_peaks$id==LD_0.2$snp1[i])]
        p2=curr_trait_peaks$logp[which(curr_trait_peaks$id==LD_0.2$snp2[i])]
        if(p1>=p2){
          peak_snp_remove_list=c(peak_snp_remove_list,LD_0.2$snp2[i])
        }else{
          peak_snp_remove_list=c(peak_snp_remove_list,LD_0.2$snp1[i])
        }
      }
    }
  }

  peak_snp_remove_list=unique(peak_snp_remove_list)
  curr_trait_peaks_cleaned=curr_trait_peaks[!curr_trait_peaks$id %in% peak_snp_remove_list,]

  curr_trait_peaks_cleaned$peak_num2=1:nrow(curr_trait_peaks_cleaned)
  write.table(curr_trait_peaks_cleaned,paste('peak_snps_v2_',t,'.csv',sep=''),sep=',',quote=F,row.names = F)

  a=merge(curr_trait_snps,curr_trait_peaks_cleaned[,c("peak_num1.2","peak_num2")],by="peak_num1.2",all=T)
  a <- a[,c(2:10,1,11)]
  a <- a[order(a$Chromosome,a$Position),]
  # a <- a[!is.na(a$peak_num2),]
  write.table(a,paste('gwas_snps_v2_',t,'.csv',sep=''),sep=',',quote=F,row.names = F)
}

######manually check all assigned loci, assign peak_num2


# loop through each traits, assign comments (peak filtering) ---------------------
for(t in traits){

  print(t)
  #peak assignment 3 of current trait
  curr_gwas=fread(paste('gwas_snps_v2_',t,'.csv',sep=''),data.table = F)
  curr_gwas$SNP_count=0
  uni_peak=max(curr_gwas$peak_num2,na.rm=TRUE)

  #loop through each peak
  for (i in 1:uni_peak){
    curr_peak=curr_gwas[which(curr_gwas$peak_num2==i),]
    count=nrow(curr_peak)
    curr_gwas$SNP_count[which(curr_gwas$peak_num2==i)]=count
  }

  #remove helicopters
  curr_gwas$comments=''
  ind=which(curr_gwas$SNP_count==1)
  if (length(ind)>0){curr_gwas$comments[ind]='helicopter'}

  #if < 5 SNPs, check if there are pairs <100K in distance
  peak_ind=unique(curr_gwas[which(curr_gwas$SNP_count<5 & curr_gwas$SNP_count>1 ),"peak_num2"])

  # peak_ind=peak_ind[!is.na(peak_ind)]
  for (i in peak_ind){
    check=curr_gwas[which(curr_gwas$peak_num2==i),'dis']
    if (min(check)>100000){
      curr_gwas$comments[which(curr_gwas$peak_num2==i)]='dis_too_large'
    }
  }

  #add annotation for merged peaks
  # curr_gwas$comments[which(is.na(curr_gwas$peak_num2))]='merged_by_LD'

  #examine to see if there are peaks merged incorrectly
  curr_gwas_peak3=curr_gwas[order(curr_gwas$Chromosome,curr_gwas$Position),]
  curr_gwas_peak3_no_na=curr_gwas_peak3[!is.na(curr_gwas_peak3$peak_num2),]

  for (i in 2:nrow(curr_gwas_peak3_no_na)){
    if(curr_gwas_peak3_no_na$peak_num2[i]<curr_gwas_peak3_no_na$peak_num2[i-1]){
      snp=curr_gwas_peak3_no_na$SNP[i]
      curr_gwas_peak3[which(curr_gwas_peak3$SNP==snp),'comments']='merged_wrong_check'
    }
  }

  #find peak SNPs
  peaks=data.frame(matrix(nrow=0,ncol=ncol(curr_gwas_peak3)))
  xx=unique(curr_gwas_peak3$peak_num2)
  if(length(which(is.na(xx)))>0){xx=xx[-which(is.na(xx))]}

  for (i in xx){
    curr=curr_gwas_peak3[which(curr_gwas_peak3$peak_num2==i),]
    s=curr[which(curr$P.value==min(curr$P.value))[1],]
    s$SNP_count=nrow(curr)
    names(peaks)=names(s)
    peaks=rbind.data.frame(peaks,s)
  }

  # names(peaks)=names(curr_gwas_peak3)   ###???# no need
  print(nrow(peaks))
  #write out v2
  write.table(peaks,paste('peak_snps_v3_',t,'.csv',sep=''),sep=',',quote=F,row.names = F)
  write.csv(curr_gwas_peak3,paste('gwas_peak_v3_',t,'.csv',sep=''),row.names = F,quote = F)
}


#plot_function -------------------------------------------------------
plot_LD=function(all,dis,snp_in_region){
  # ind=which(all_snps$`R^2`==1)
  ind=which(all_snps$Position==LD$Position[1])
  all=all_snps[which(all_snps$Position>=(all_snps$Position[ind]-dis*1000000) & all_snps$Position<=(all_snps$Position[ind]+dis*1000000)),]

  a<-range(-log10(all$P.value))[2]+1
  b<-all$snp_pos[1]/10^6 - 0.1 - dis #dis in Mb
  c<-all$snp_pos[1]/10^6 + 0.1 + dis#dis in Mb
  # pos=all$Position[which(all$`R^2`==1)]
  pos <- LD$Position[1]

  par(mar = c(4,5,2,4))
  plot(NULL, xlim=c(b,c), ylim=c(0,a), ylab=expression(-log[10](italic(p))), xlab="Physical position (Mb)",main=paste('trait: ',all$trait[1],', SNP: ',all$SNP[1],sep=''))
  rect(xleft=(pos-100000)/10^6,xright=(pos+100000)/10^6,ybottom=0,ytop=(a+0.8),col='lavender',border = NA)

  for (m in 1:nrow(all)){
    logp=-log10(all$P.value)[m]
    xpos=all$Position[m]/10^6
    segments(x0=xpos,x1=xpos,y0=0,y1=logp,col='gray85',lty=1,ylim=c(0,a),lwd=1)
  }

  #plot again snps with significance
  all_sig=all[which(-log10(all$P.value)>(-log10(all$FDR_thres0.05[1]))),]

  if (nrow(all_sig)>0){
    for (m in 1:nrow(all_sig)){
      logp=-log10(all_sig$P.value)[m]
      xpos=all_sig$Position[m]/10^6
      segments(x0=xpos,x1=xpos,y0=0,y1=logp,col="black",lty=1,ylim=c(0,a),lwd=1.5)
    }
  }

  #plot other peak snps within selected interval
  if (nrow(snp_in_region)>1){
    snp_in_region2=snp_in_region[-which(snp_in_region$Position==pos),]
    for (x in 1:nrow(snp_in_region2)){
      logp=-log10(snp_in_region2$P.value)[x]
      xpos=snp_in_region2$Position[x]/10^6
      segments(x0=xpos,x1=xpos,y0=0,y1=logp,col="orange",lty=4,ylim=c(0,a),lwd=1.5)
    }
  }

  #add significance line
  abline(h=-log10(all$FDR_thres0.05[1]),col="red",lty=2)

  #new plot for LD
  par(new = T)
  plot(all$Position/10^6,all$`R^2`,
       # pch=ifelse(all$Position==pos,24,17),
       pch=17,
       cex=ifelse(all$Position==pos, 0.75, ifelse(all$`R^2`>0.4, 0.6,0.45)),
       col=ifelse(all$Position==pos, "red", ifelse(all$`R^2`>0.4, 'royalblue','powderblue')),
       axes=F, xlab=NA, ylab=NA,xlim=c(b,c),ylim=c(0,1))


  #plot other peak snps within selected interval
  if (nrow(snp_in_region)>1){
    par(new = T)
    snp_ld=all[all$Position %in% snp_in_region2$Position,]
    plot(snp_ld$Position/10^6,snp_ld$`R^2`,
         pch=8,cex=1,
         col='orange',
         axes=F, xlab=NA, ylab=NA,xlim=c(b,c),ylim=c(0,1))
  }

  axis(side = 4)
  mtext(side = 4, line = 3, expression(paste("r"^"2")))

}


#end ---------------------------------------------------------------------

# plot final peak selection -----------------------------------------------
system('mkdir peak_snp_LD')
setwd('./peak_snp_LD')
#calculate LD to peak SNPs +-5M

for(t in traits){
  print(t)
  peak_snp=read.csv(paste('../peak_snps_v3_',t,'.csv',sep=''))

  #get physical position +/- 1M of the peak snp
  peak_snp$minus1M=peak_snp$Position-1000000
  peak_snp$plus1M=peak_snp$Position+1000000

  #get snps from Hapmap and calculate LD
  for (i in 1:nrow(peak_snp)){
    curr_ch=peak_snp$Chr[i]
    pos=peak_snp$Position[i]
    up=peak_snp$plus1M[i]
    low=peak_snp$minus1M[i]
    snp=peak_snp$SNP[i]
    #parse out all SNPs in region
    # system(paste('vcftools --gzvcf /workdir/dw524/ionomics/Hapmap_genotype/post_imputation_filtering/3.filtering/AGPv4_Ames_ionomics',curr_ch,'_maf.vcf.gz --chr ',curr_ch,' --from-bp ',low ,' --to-bp ',up,' --recode --out sel_snps_SNP_',snp,'_',trait,' &',sep=''))
    system(paste('vcftools --gzvcf /workdir/dw524/vitamaize/genotype/3.filtering/AGPv4_Ames_vitamaize',curr_ch,'_maf.vcf.gz --chr ',curr_ch,' --from-bp ',low ,' --to-bp ',up,' --recode --out sel_snps_SNP_',snp,'_',t,' &',sep=''))

  }
}

for(t in traits){
  print(t)
  peak_snp=read.csv(paste('../peak_snps_v3_',t,'.csv',sep=''))

  for (i in 1:nrow(peak_snp)){
    curr_ch=peak_snp$Chr[i]
    pos=peak_snp$Position[i]
    snp=peak_snp$SNP[i]
    #read in SNPs
    # a=fread(paste('sel_snps_SNP_',snp,'_',t,'.recode.vcf',sep=''),colClasses=c(rep(NA,3),rep("NULL",1819)))
    a=fread(paste('sel_snps_SNP_',snp,'_',t,'.recode.vcf',sep=''),colClasses=c(rep(NA,3),rep("NULL",1468)))   # 1468 == 1462 (lines count) + 6

    #calculate LD
    index=grep(pos,a$POS)-1
    system(paste('/programs/tassel-5-standalone/run_pipeline.pl -Xmx350g -vcf sel_snps_SNP_',snp,'_',t,'.recode.vcf -ld -ldType SiteByAll -ldTestSite ',index, ' -export ',snp,'_site_by_all.txt &',sep=''))
  }
}

for(t in traits){

  print(t)
  peak_snp=read.csv(paste('../peak_snps_v3_',t,'.csv',sep=''),stringsAsFactors = F)

  #get GWAS results for selected SNPs of the traits
  pdf(paste('../gwas_peak_snp_assignment_',t,'_.pdf',sep=''),family='serif')
  gwas_res=fread(paste('../../GAPIT_GWAS_Result_all_chr_',t,'.csv',sep=''),data.table = F)
  snps=unique(peak_snp$SNP)

  for (j in snps){
    chr=peak_snp$Chr[which(peak_snp$SNP==j)]
    gwas_res_curr_chrom=gwas_res[which(gwas_res$Chromosome==chr),]
    snp_pos=as.numeric(peak_snp$Position[which(peak_snp$SNP==j & peak_snp$trait==t)])
    LD=fread(paste(j,'_site_by_all.txt',sep=''))
    #format LD
    LD$Position=LD$Position2
    LD[which(LD$Position2==snp_pos),'Position']=LD[which(LD$Position2==snp_pos),'Position1']
    LD=LD[,c("Position","R^2")]
    #add back the peak SNP
    LD=rbind(t(matrix(c(snp_pos,1))),LD,use.names=FALSE)
    names(LD)=c("Position","R^2")

    #merge with GWAS results
    all_snps=merge(LD,gwas_res_curr_chrom[,c('Position','P.value','FDR_thres0.05')])
    all_snps$SNP=j
    all_snps$trait=t
    all_snps$snp_pos=snp_pos

    #find if any other peaks in the interval
    snp_in_region=peak_snp[which(peak_snp$Chromosome==chr & peak_snp$Position> (snp_pos-5000000) & peak_snp$Position< (snp_pos+5000000)),]

    #plot LD
    plot_LD(all_snps,1,snp_in_region)
  }
  dev.off()
}



# get candidate gene list -------------------------------------------------
library(qdapRegex)
library(reshape2)

### gene position
gff=fread('/workdir/dw524/AGPv4_genes_Gramene/Zea_mays.B73_RefGen_v4.59_anno.csv',data.table = F)
# names(gff)=c('chr','source','type','start','end','score','strand','phase','attribute')

gff$chr=as.numeric(as.character(gff$chr))
gff$start=as.numeric(as.character(gff$start))
gff$end=as.numeric(as.character(gff$end))
# gff=gff[which(gff$type=='gene'),] #retain only gene models

# id_sp=colsplit(unlist(gff[,'attribute']),';',c('ID','other'))#split to retain only gene id
# gff=cbind.data.frame(gff,id_sp)

for(t in traits){
  print(t)

  ### import GWAS peaks
  peak_snp=read.csv(paste('../peak_snps_v3_',t,'.csv',sep=''))

  #annotation file
  snp_annotated=as.data.frame(matrix(NA,nrow=0,ncol=(ncol(peak_snp)+ncol(gff))))
  names(snp_annotated)=c(names(peak_snp),names(gff))

  for (i in 1:nrow(peak_snp)){
    print(i)
    snp_pos=peak_snp$Position[i]
    snp_chr=peak_snp$Chromosome[i]
    #search region left bound
    left=min(c(peak_snp$r2_02_l[i],peak_snp$r2_05_l[i],snp_pos-100000))
    #search region right bound
    right=max(c(peak_snp$r2_02_r[i],peak_snp$r2_05_r[i],snp_pos+100000))

    #pull out genes within search region
    index=which(gff$chr==snp_chr & gff$start<=right & gff$end>=left) #partial overlap allowed
    if (length(index)!=0){
      curr=as.data.frame(c(peak_snp[i,],gff[index,]))
    }else{
      curr=as.data.frame(c(peak_snp[i,],rep(NA,ncol(gff))))
    }
    names(curr)=names(snp_annotated)
    snp_annotated=rbind.data.frame(snp_annotated,curr)
  }

  #add annotation and gene-to-snp distance
  snp_annotated$rel_pos=NA
  for (i in 1:nrow(snp_annotated)){
    pos=snp_annotated$Position[i]
    # snp_annotated[i,'description']=rm_between(snp_annotated[i,'attribute'], "description=", ";gene_id=", extract=TRUE)
    # snp_annotated$ID=gsub('ID=gene:','',snp_annotated$ID)
    snp_annotated$'distance to ORF_start in bp'=snp_annotated$Position-snp_annotated$start
    snp_annotated$'distance to ORF_end in bp'=snp_annotated$Position-snp_annotated$end
    if(!is.na(snp_annotated$"distance to ORF_start in bp"[i])){
      if (snp_annotated$'distance to ORF_start in bp'[i]>0 & snp_annotated$'distance to ORF_end in bp'[i]>0){snp_annotated$rel_pos[i]='downstream'}
      if (snp_annotated$'distance to ORF_start in bp'[i]>0 & snp_annotated$'distance to ORF_end in bp'[i]<0){snp_annotated$rel_pos[i]='within'}
      if (snp_annotated$'distance to ORF_start in bp'[i]<0 & snp_annotated$'distance to ORF_end in bp'[i]<0){snp_annotated$rel_pos[i]='upstream'}
    }
  }
  # #r2 0.2 genes
  # snp_annotated$LD02='N'
  # snp_annotated$LD02[which(snp_annotated$start <= snp_annotated$r2_02_r & snp_annotated$end >= snp_annotated$r2_02_l)]='Y'#partial overlap is allowed
  #
  # snp_annotated$LD02[snp_annotated$r2_02_r==snp_annotated$Position & snp_annotated$start>snp_annotated$Position]=0
  # snp_annotated$LD02[snp_annotated$r2_02_l==snp_annotated$Position & snp_annotated$end<snp_annotated$Position]=0
  #
  # #r2 0.5 genes
  # snp_annotated$LD05='N'
  # snp_annotated$LD05[which(snp_annotated$start <= snp_annotated$r2_05_r & snp_annotated$end >= snp_annotated$r2_05_l)]='Y'
  #
  # snp_annotated$LD05[snp_annotated$r2_05_r==snp_annotated$Position & snp_annotated$start>snp_annotated$Position]=0
  # snp_annotated$LD05[snp_annotated$r2_05_l==snp_annotated$Position & snp_annotated$end<snp_annotated$Position]=0
  #
  #r2 100K
  snp_annotated$`100K`='N'
  snp_annotated$`100K`[which(snp_annotated$start <= (snp_annotated$Position + 100000) & snp_annotated$end >= (snp_annotated$Position - 100000))]='Y'

  names(snp_annotated)[1]='snp_number'
  write.csv(snp_annotated, paste('../candidate_genes_peak_snps_v3_',t,'.csv',sep=''),row.names=F)
}


# r2 of SNP to gene -------------------------------------------------------
for(t in traits){
  print(t)
  cand=fread(paste('../candidate_genes_peak_snps_v3_',t,'.csv',sep=''),data.table = F)
  gwas_res=fread(paste('../../GAPIT_GWAS_Result_all_chr_',t,'.csv',sep=''),data.table = F)
  cand_no_missing=cand[!is.na(cand$ID),]#skip SNPs without candidate genes in the interval
  snps=unique(cand_no_missing$snp_number)
  for(j in snps){
    print(j)
    curr=cand_no_missing[which(cand_no_missing$snp_number==j),]
    snp_pos=cand[which(cand$snp_number==j),'Position'][1]
    LD=fread(paste('../peak_snp_LD/',j,'_site_by_all.txt',sep=''), data.table=F)
    #format LD
    LD$Position=LD$Position2
    LD[which(LD$Position2==snp_pos),'Position']=LD[which(LD$Position2==snp_pos),'Position1']
    LD=LD[,c("Position","R^2")]
    #merge with GWAS result
    ch=curr$chr[1]
    LD=merge(LD,gwas_res[which(gwas_res$Chromosome==ch),c("Position","P.value")],by='Position',all.x=T)
    LD=LD[order(LD$Position),]
    for(x in 1:nrow(curr)){
      gene=curr$ID[x]
      st=curr$start[x]
      end=curr$end[x]
      sub=LD[which(LD$Position<=end & LD$Position>=st),]
      cand[which(cand$snp_number==j & cand$ID==gene),'num_SNP_in_gene']=nrow(sub)
      cand[which(cand$snp_number==j & cand$ID==gene),'r2_min']=min(sub$`R^2`)
      cand[which(cand$snp_number==j & cand$ID==gene),'r2_max']=max(sub$`R^2`)
      cand[which(cand$snp_number==j & cand$ID==gene),'r2_mean']=mean(sub$`R^2`)
      cand[which(cand$snp_number==j & cand$ID==gene),'r2_median']=median(sub$`R^2`)
      cand[which(cand$snp_number==j & cand$ID==gene),'p_min']=min(sub$`P.value`)
      cand[which(cand$snp_number==j & cand$ID==gene),'p_max']=max(sub$`P.value`)
      cand[which(cand$snp_number==j & cand$ID==gene),'p_mean']=mean(sub$`P.value`)
      cand[which(cand$snp_number==j & cand$ID==gene),'p_median']=median(sub$`P.value`)
      if(nrow(sub)==0){
        if(snp_pos>end){
          sub_alt=LD[which(LD$Position>end),][1:10,]
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_min']=min(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_max']=max(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_mean']=mean(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_median']=median(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_min']=min(sub_alt$`P.value`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_max']=max(sub_alt$`P.value`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_mean']=mean(sub_alt$`P.value`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_median']=median(sub_alt$`P.value`)
        }else if (snp_pos<st){
          sub_alt=LD[which(LD$Position<st),]
          sub_alt=sub_alt[(nrow(sub_alt)-9):nrow(sub_alt),]
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_min']=min(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_max']=max(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_mean']=mean(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'r2_median']=median(sub_alt$`R^2`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_min']=min(sub_alt$`P.value`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_max']=max(sub_alt$`P.value`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_mean']=mean(sub_alt$`P.value`)
          cand[which(cand$snp_number==j & cand$ID==gene),'p_median']=median(sub_alt$`P.value`)
        }
      }
    }
  }
  fwrite(cand,paste('../candidate_genes_peak_snps_v3_',t,'_manually_curated_r2_added.csv',sep=''))
}
