options(stringsAsFactors = F)

require(data.table)
require(tidyr)
require(vcd)
require(BEDMatrix)
require(snow)
require(dplyr)
require(gtools)

######### all snv in WGS  ##############
setwd("/home/sshen/mount/t7920/DiskC/plink/zstd")
for (i in 1:22) {
  chr = fread(paste0("chr",i,".afreq"))
  chr$maf = ifelse(chr$ALT_FREQS <= 0.5,chr$ALT_FREQS,1-chr$ALT_FREQS)
  chr$pos = unlist(lapply(chr$ID,function(x) strsplit(x,":")[[1]][2]))
  chr = unite(chr,id,c("#CHROM","pos","REF","ALT"),sep = ":",remove = T)
  chr = chr[chr$maf<0.01,]
  write.table(chr$id,file=paste0("/home/biostat/DiskB/jldai/GPB_revise/White/false_positive/WGS_varID/stat0.01/chr",i,".txt"),
              col.names = F,row.names = F,quote = F)
}

######## negative snv in WGS #######
setwd("/home/biostat/DiskA/resource/UKB/BED/biochemistry_WGS")
## p：5e-8
threshold = -log10(5e-08)
nagetive_table = NULL
for (i in 1:30) {
  pheno = paste0("c",i)
  now_sum = 0
  for (j in 1:22) {
    chr_biochemistry = fread(paste0("chr",j,"_C",i,".regenie"))
    chr_biochemistry = select(chr_biochemistry,ID,ALLELE0,ALLELE1,A1FREQ,LOG10P)
    chr_biochemistry$maf = ifelse(chr_biochemistry$A1FREQ>=0.5,chr_biochemistry$A1FREQ,1-chr_biochemistry$A1FREQ)
    chr_biochemistry_keep = chr_biochemistry[chr_biochemistry$LOG10P<=threshold & chr_biochemistry$maf<0.01,]
    now_sum = now_sum+nrow(chr_biochemistry_keep)
  }
  temp = data.frame(pheno = pheno,nagetive_sum = now_sum)
  nagetive_table = rbind(nagetive_table,temp)
}
write.table(nagetive_table,file = "/home/biostat/DiskB/jldai/GPB_revise/White/false_positive/keep_inWGS/WGS_negative_stat.txt",col.names = T, row.names = F)

######## significant snv of WGS ###########
setwd("/home/biostat/DiskB/jldai/GPB_revise/White/false_positive")
wgs = fread("WGS_chemi_single.txt")

######## significant snv of HRC #################
setwd("/home/biostat/DiskB/jldai/GPB_revise/White/false_positive")
hrc = fread("HRC_chemi_single.txt")
raw_id = hrc$ID
keep_id = NULL
for (i in 1:22) {
  now = fread(paste0("WGS_varID/chr",i,".txt"),header = F)
  wgs_id = now$V1
  temp_id = raw_id[raw_id %in% wgs_id]
  keep_id = c(keep_id,temp_id)
}
write.table(keep_id,file="keep_inWGS/HRC_keepID.txt",col.names = F,row.names = F,quote = F)

###### stat #####
setwd("/home/biostat/DiskB/jldai/GPB_revise/White/false_positive")
wgs = fread("WGS_chemi_single.txt")
wgs = select(wgs,ID,TRAIT,panel)
hrc = fread("HRC_chemi_single.txt")
hrc = select(hrc,ID,TRAIT,panel)
keep = fread("keep_inWGS/HRC_keepID.txt",header = F)
keep_id = keep$V1
hrc = hrc[hrc$ID %in% keep_id,]

all = NULL
for (i in 1:30) {
  trait = paste0("c",i)
  now_wgs = wgs[wgs$TRAIT==trait,c("ID")]
  wgs_id = now_wgs$ID
  now_hrc = hrc[hrc$TRAIT==trait,c("ID")]
  hrc_id = now_hrc$ID
  
  if (length(hrc_id)>0 & length(wgs_id)>0) {
    wgs_hrc = intersect(wgs_id,hrc_id)
    hrc_false = length(hrc_id)-length(wgs_hrc)
    false_ratio = hrc_false /length(hrc_id)
    temp = data.frame(trait = trait,hrc_false = hrc_false,
                      hrc_find = length(hrc_id),wgs_find = length(wgs_id),false_positive = round(false_ratio,4))
    temp$panel = "HRC+UK10K"
    all = rbind(all,temp)
  }
}

label = fread("/home/biostat/DiskA/resource/UKB/BED/total_result/biochemistry/biochemistry_datafields.csv")
label = select(label,label,Description)
label$label = tolower(label$label)
label = label[mixedorder(label),]
colnames(label) = c("trait","label")
all = inner_join(all,label,by="trait")

temp = data.frame(trait="average",
                  hrc_false = NA,
                  hrc_find = NA,
                  wgs_find = NA,
                  false_positive=round(mean(all$false_positive,na.rm=T),4),
                  panel = "HRC+UK10K",
                  label = "")
all = rbind(all,temp)
all = select(all,trait,label,hrc_false,hrc_find,wgs_find,false_positive,panel)

write.table(all,file = "keep_inWGS/(abs)HRC_WGS_stat-single.csv",sep = ",",
            col.names = T,row.names = F,quote = F)


######## significant snv of TOPMed #################
setwd("/home/biostat/DiskB/jldai/GPB_revise/White/false_positive")
top = fread("TOP_chemi_single.txt")
raw_id = top$ID
keep_id = NULL
for (i in 1:22) {
  now = fread(paste0("WGS_varID/chr",i,".txt"),header = F)
  wgs_id = paste0("chr",now$V1)
  temp_id = raw_id[raw_id %in% wgs_id]
  keep_id = c(keep_id,temp_id)
}
write.table(keep_id,file="keep_inWGS/TOP_keepID.txt",col.names = F,row.names = F,quote = F)

##### stat ##########
setwd("/home/biostat/DiskB/jldai/GPB_revise/White/false_positive")
wgs = fread("WGS_chemi_single.txt")
wgs = select(wgs,ID,TRAIT,panel)
top = fread("TOP_chemi_single.txt")
top = select(top,ID,TRAIT,panel)
keep = fread("keep_inWGS/TOP_keepID.txt",header = F)
keep_id = keep$V1
top = top[top$ID %in% keep_id,]

all = NULL
for (i in 1:30) {
  trait = paste0("c",i)
  now_wgs = wgs[wgs$TRAIT==trait,c("ID")]
  wgs_id = now_wgs$ID
  wgs_id = paste0("chr",wgs_id)
  now_top = top[top$TRAIT==trait,c("ID")]
  top_id = now_top$ID
  
  if (length(top_id)>0 & length(wgs_id)>0) {
    wgs_top = intersect(wgs_id,top_id)
    top_false = length(top_id)-length(wgs_top)
    false_ratio = top_false /length(top_id)
    temp = data.frame(trait = trait,top_false = top_false,
                      top_find = length(top_id),wgs_find = length(wgs_id),false_positive = round(false_ratio,4))
    temp$panel = "TOPMed"
    all = rbind(all,temp)
  }
}

label = fread("/home/biostat/DiskA/resource/UKB/BED/total_result/biochemistry/biochemistry_datafields.csv")
label = select(label,label,Description)
label$label = tolower(label$label)
label = label[mixedorder(label),]
colnames(label) = c("trait","label")
all = inner_join(all,label,by="trait")

temp = data.frame(trait="average",
                  top_false = NA,
                  top_find = NA,
                  wgs_find = NA,
                  false_positive=round(mean(all$false_positive,na.rm=T),4),
                  panel = "TOPMed",
                  label = "")
all = rbind(all,temp)
all = select(all,trait,label,top_false,top_find,wgs_find,false_positive,panel)

write.table(all,file = "keep_inWGS/(abs)TOP_WGS_stat-single.csv",sep = ",",
            col.names = T,row.names = F,quote = F)


