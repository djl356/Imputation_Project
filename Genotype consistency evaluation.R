## required packages ## 
options(stringsAsFactors = F)

require(data.table)
require(tidyr)
require(vcd)
require(BEDMatrix)
require(snow)

## Ethnicity classification for UKB SNP array data
load("/home/biostat/DiskA/resource/UKB/ukb20220804.RData")
## based on "ethnic" field ##
White = ukb[grep("^1",ukb$ethnic),]
Asian = ukb[grep("^3",ukb$ethnic),]
Black = ukb[grep("^4",ukb$ethnic),]

### Use plink to extract the genetic data of the three ethnicities
### in the impute data (480k) and calculate the MAC respectively

## Ethnicity classification for WGS data
WGS_all = fread("WGS_150k.sample")
colnames(WGS_all) = c("ID1","ID2")
WGS_white = WGS_all[WGS_all$ID2 %in% White$ID,]
WGS_asian = WGS_all[WGS_all$ID2 %in% Asian$ID,]
WGS_black = WGS_all[WGS_all$ID2 %in% Black$ID,]

### Use plink to extract the genetic data of the three ethnicities
### in the WGS data (150k) and calculate the MAC respectively

## Divide the MAC interval (conducted separately for the three ethnicities)
## Randomly select 5000 variants in each interval of each chromosome
## Note:should change the path for differnt ethnic groups
for (i in 1:22) {
  wgs_count = fread(paste0("/home/biostat/DiskB/jldai/White/WGS_counts/chr",i,".acount"))
  wgs_count = wgs_count[wgs_count$ALT_CTS>0,]

    type = cut(wgs_count$ALT_CTS,breaks=c(0,5,10,20,50,3000),
             labels=c("(0,5]","(5,10]","(10,20]","(20,50]",">50"))
  wgs_countType = cbind(wgs_count,type)
  
  wgs_rare1 = na.omit(wgs_count[wgs_countType$type=="(0,5]",2])
  wgs_rare2 = na.omit(wgs_count[wgs_countType$type=="(5,10]",2])
  wgs_rare3 = na.omit(wgs_count[wgs_countType$type=="(10,20]",2])
  wgs_rare4 = na.omit(wgs_count[wgs_countType$type=="(20,50]",2])
  wgs_rare5 = na.omit(wgs_count[wgs_countType$type==">50 and MAF<0.001",2])
  
  
  wgs_rareSample1 = sample(wgs_rare1$ID,5000)
  wgs_rareSample2 = sample(wgs_rare2$ID,5000)
  wgs_rareSample3 = sample(wgs_rare3$ID,5000)
  wgs_rareSample4 = sample(wgs_rare4$ID,5000)
  wgs_rareSample5 = sample(wgs_rare5$ID,5000)
  
  write.table(wgs_rareSample1,paste0("rare_sample1/chr",i,"_WGS.txt"), quote=F,col.names=F,row.names=F,sep="\t")
  write.table(wgs_rareSample2,paste0("rare_sample2/chr",i,"_WGS.txt"), quote=F,col.names=F,row.names=F,sep="\t")
  write.table(wgs_rareSample3,paste0("rare_sample3/chr",i,"_WGS.txt"), quote=F,col.names=F,row.names=F,sep="\t")
  write.table(wgs_rareSample4,paste0("rare_sample4/chr",i,"_WGS.txt"), quote=F,col.names=F,row.names=F,sep="\t")
  write.table(wgs_rareSample5,paste0("rare_sample5/chr",i,"_WGS.txt"), quote=F,col.names=F,row.names=F,sep="\t")
}

### Use plink to extract the genetic data of each MAC interval 
### for the three ethnicity in the WGS data based on the extracted variants above

## Extract variants in the TOPMed-imputed data that match WGS variants that have been extracted before 
## This following step is used to generate variants files for subsequent plink extraction
## Note:should change the path for differnt ethnic groups
## Take White as a example here
setwd("/home/biostat/DiskB/jldai/White")
for (i in 1:22) {
  TopMed = fread(paste0("/home/biostat/DiskD/TopMed_part1/ukb_topmed_c",i,"_qc.bim"))
  for (j in 1:5) {
    rareSample = read.table(paste0("rare_sample/rare_sample",j,"/chr",i,"_WGS.txt"),header = F)
    rareSample = sapply(strsplit(rareSample$V1,":"),"[",2)   
    rare_TopMed = TopMed[TopMed$V4 %in% rareSample,V2]
    write.table(rare_TopMed,paste0("rare_sample/rare_sample",j,"/chr",i,"_TopMed.txt"),quote=F,col.names=F,row.names=F,sep="\t")
  }
}

### Use plink to extract the genetic data of each MAC interval 
### for the three ethnicity in the TOPMed-imputed data based on the extracted variants above

## Extract variants in the HRC+UK10K-imputed data that match WGS variants that have been extracted before 
## This step is used to generate variants files for subsequent plink extraction
## Note:should change the path for differnt ethnic groups
setwd("/home/biostat/DiskB/jldai/White")
for (i in 1:22) {
  GWAS = fread(paste0("/home/biostat/DiskA/resource/UKB/GWAS/chr",i,"_qc.bim"))
  for (j in 1:5) {
    rareSample = read.table(paste0("rare_sample/rare_sample",j,"/chr",i,"_WGS.txt"),header = F)
    rareSample = sapply(strsplit(rareSample$V1,":"),"[",2)   
    rare_GWAS = GWAS[GWAS$V4 %in% rareSample,V2]
    write.table(rare_GWAS,paste0("rare_sample/rare_sample",j,"/chr",i,"_GWAS.txt"), quote=F,col.names=F,row.names=F,sep="\t")
  }
}

### Use plink to extract the genetic data of each MAC interval 
### for the three ethnicity in the TOPMed-imputed data based on the extracted variants above

################################## corr_eval function ######################################################
## This function is used to calculate the correlation between the variant of the imputed data 
## and the corresponding variant of the WGS data
## Parameter 'p':The position of the variant to be tested
## Parameter 'a':The bed file of the current MAC interval of the current chromosome of WGS data
## Parameter 'b':The bed file of the current MAC interval of the current chromosome of TOPMed-imputed data
## Parameter 'c':The bed file of the current MAC interval of the current chromosome of HRC+UK10K-imputed data
corr_eval = function(p,a,b,c){
  require(vcd)
  pos_now = p
  
  ## find index of position
  id_wgs = grep(pos_now,colnames(a));if(length(id_wgs)>1) id_wgs=id_wgs[1]
  id_top = grep(pos_now,colnames(b));if(length(id_top)>1) id_top=id_top[1]
  id_gwas = grep(pos_now,colnames(c));if(length(id_gwas)>1) id_gwas=id_gwas[1]
  
  ## match imputed data samples to WGS samples
  wgs_top = a[match(rownames(b),rownames(a)),]
  wgs_gwas = a[match(rownames(c),rownames(a)),]
  
  ## convert to numeric vector
  x = as.numeric(a[,id_wgs])
  x1 = as.numeric(wgs_top[,id_wgs])
  x2 = as.numeric(wgs_gwas[,id_wgs])
  y = as.numeric(b[,id_top])
  z = as.numeric(c[,id_gwas])
  
  ## caculate MAF and MAC
  MAF = 0.5*mean(x,na.rm=T)
  MAC = sum(x>0,na.rm=T)
  temp = c(pos_now,MAF,MAC)
  
  ## TopMed vs WGS
  t_l1 = length(unique(x1)[!is.na(unique(x1))])
  t_l2 = length(unique(y)[!is.na(unique(y))])
  
  if(t_l1==0 | t_l2==0){
    t1=1;t2=1
  }else t1=t_l1;t2=t_l2
  
  if(length(y)>0 & (t1!=1 | t2!=1)){
    stat = assocstats(table(x1,y))
    t_cramer = round(stat$cramer,digits = 2)
    temp = c(temp,t_cramer)
  }else temp = c(temp,t_cramer=NA)
  
  ## HRC+UK10K vs WGS
  g_l1 = length(unique(x2)[!is.na(unique(x2))])
  g_l2 = length(unique(z)[!is.na(unique(z))])
  
  if(g_l1==0 | g_l2==0){
    g1=1;g2=1
  }else g1=g_l1;g2=g_l2
  
  if(length(z)>0 & (g1!=1 | g2!=1)){
    stat = assocstats(table(x2,z))
    g_cramer = round(stat$cramer,digits = 2)
    temp = c(temp,g_cramer)
  }else temp = c(temp,g_cramer=NA)
  
  return(temp)
}
############################################################################################################

## Loop of the consistence analysis
## Outer loop: 22 chromosomes ;Inner loop: 5 MAC intervals
## Note:should change the path for differnt ethnic groups
## Take White as a example here
setwd("/home/sshen/EPYC/DiskB/jldai/White")
cl = makeCluster(15)
for (i in 1:22) {
  for (j in 1:5) {
    ## WGS 
    wgs_bed = as.matrix(BEDMatrix(paste0("interval",j,"/WGS/chr",i)),simple_names=T)
    wgs_col = colSums(wgs_bed,na.rm=T)
    rownames(wgs_bed) = gsub("0_","",rownames(wgs_bed))
    
    ## filter variant 
    wgs_col = subset(wgs_col,wgs_col!=nrow(wgs_bed) & wgs_col!= 2*nrow(wgs_bed))
    wgs_bed = wgs_bed[,wgs_col>0]
    
    ## TopMed
    top_bed = as.matrix(BEDMatrix(paste0("interval",j,"/TopMed/chr",i),simple_names=T))
    raw_top = rownames(top_bed)
    rownames(top_bed) = paste0(raw_top,"_",raw_top)
    top_col = colSums(top_bed,na.rm=T)
    
    ## filter variant 
    top_col = subset(top_col,top_col!=nrow(top_bed) & top_col!= 2*nrow(top_bed))
    Tcol_name = names(top_col[top_col>0])
    top_bed = top_bed[,Tcol_name]
    
    ## GWAS
    gwas_bed = as.matrix(BEDMatrix(paste0("interval",j,"/GWAS/chr",i),simple_names=T))
    raw_gwas = rownames(gwas_bed)
    rownames(gwas_bed) = paste0(raw_gwas,"_",raw_gwas)  
    gwas_col = colSums(gwas_bed,na.rm=T)
    
    ## filter variant 
    gwas_col = subset(gwas_col,gwas_col!=nrow(gwas_bed) & gwas_col!= 2*nrow(gwas_bed))
    Gcol_name = names(gwas_col[gwas_col>0])
    gwas_bed = gwas_bed[,Gcol_name]
    
    ## get the position of variants to be tested 
    wgs_pos = unlist(lapply(colnames(wgs_bed),function(x) strsplit(x,":")[[1]][2]))
    top_pos = unlist(lapply(colnames(top_bed),function(x) strsplit(x,":")[[1]][2]))
    gwas_pos = unlist(lapply(colnames(gwas_bed),function(x) strsplit(x,":")[[1]][2]))
    ## rename
    colnames(wgs_bed) = wgs_pos
    colnames(top_bed) = top_pos
    colnames(gwas_bed) = gwas_pos
    
    ## get intersection of imputed data and wgs data
    wgs_top = intersect(wgs_pos,top_pos)
    wgs_gwas = intersect(wgs_pos,gwas_pos)
    top_bed = top_bed[,wgs_top]
    gwas_bed = gwas_bed[,wgs_gwas]
    
    ## keep the variants that intersect with WGS data in each imputed data respectively
    top_if = wgs_pos %in% wgs_top
    gwas_if = wgs_pos %in% wgs_gwas
    wgs_pos = wgs_pos[top_if | gwas_if]
    wgs_bed = wgs_bed[,wgs_pos]
    
    ## parallel computing
    corr_result = parLapply(cl,wgs_pos,corr_eval,wgs_bed,top_bed,gwas_bed)
    
    ## process result
    corr_result = t(data.frame(corr_result))
    colnames(corr_result) = c("pos","MAF","MAC",
                              "t_cramer","g_cramer")
    corr_result = as.data.frame(corr_result,row.names = FALSE)
    corr_result$pos = paste0("chr",i,":",corr_result$pos)

      for (k in 2:ncol(corr_result)) {
      corr_result[,k]= as.numeric(corr_result[,k])
    }
    write.table(corr_result, file = paste0("interval",j,"/corr_result/result_chr",i,".csv"),
                col.names = T, row.names = F,sep = ",")
  }
}
stopCluster(cl)

## Integrate the result of 22 chromosomes
## Note:should change the path for differnt ethnic groups
## Take White as a example here
setwd("/home/biostat/DiskB/jldai/White")
for (i in 1:5) {
  corr_all = data.frame()
  for (j in 1:22) {
    temp = read.csv(paste0("interval",i,"/corr_result/result_chr",j,".csv"),header = T,sep = ",")
    corr_all = rbind(corr_all,temp)
  }
  write.table(corr_all, file = paste0("interval",i,"/corr_result/result_all.csv"),
              col.names = T, row.names = F,sep = ",")
}


