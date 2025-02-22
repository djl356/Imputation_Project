
require(data.table);require(openxlsx)
setwd("/home/sshen/EPYC/DiskB/jldai/WGS_LC_OV")


########### lung cancer
dta_ukb = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/EPYC/DiskA/resource/UKB/BED/cancer_WGS/chr",i,"_LC.regenie"))
# temp$mac = round(2*temp$N*temp$A1FREQ,0)
temp = subset(temp,A1FREQ<=0.01 & !is.na(LOG10P))
dta_ukb = rbind(dta_ukb,temp)
}
dta_ukb$label = paste(dta_ukb$CHROM,dta_ukb$GENPOS,sep=":")

dta_onco = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/DiskB/onco/SAIGE/C",i,".txt"))
temp$maf = ifelse(temp$AF_Allele2>0.5,1-temp$AF_Allele2,temp$AF_Allele2)
temp = subset(temp,maf<0.02)
dta_onco = rbind(dta_onco,temp)
}
dta_onco$label = paste(dta_onco$CHR,dta_onco$POS,sep=":")
dta_onco$N = dta_onco$N_case+dta_onco$N_ctrl

dta_plco = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/DiskB/PLCO/SAIGE_new/C",i,".txt"))
temp$maf = ifelse(temp$AF_Allele2>0.5,1-temp$AF_Allele2,temp$AF_Allele2)
temp = subset(temp,maf<0.02)
dta_plco = rbind(dta_plco,temp)
}
dta_plco$label = paste(dta_plco$CHR,dta_plco$POS,sep=":")
dta_plco$N = dta_plco$N_case+dta_plco$N_ctrl

dta_TRICL = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/DiskB/TRICL/SAIGE/C",i,".txt"))
temp$maf = ifelse(temp$AF_Allele2>0.5,1-temp$AF_Allele2,temp$AF_Allele2)
temp = subset(temp,maf<0.02)
dta_TRICL = rbind(dta_TRICL,temp)
}
dta_TRICL$label = paste(dta_TRICL$CHR,dta_TRICL$POS,sep=":")
dta_TRICL$N = dta_TRICL$N_case+dta_TRICL$N_ctrl


#### METAL
dta_ukb$p.value = 10^(-dta_ukb$LOG10P)
a = dta_ukb[,c("label","ALLELE1","ALLELE0","A1FREQ","N","BETA","SE","p.value")]
b = dta_onco[,c("label","Allele2","Allele1","AF_Allele2","N","BETA","SE","p.value")]
c = dta_plco[,c("label","Allele2","Allele1","AF_Allele2","N","BETA","SE","p.value")]
d = dta_TRICL[,c("label","Allele2","Allele1","AF_Allele2","N","BETA","SE","p.value")]
colnames(a) = colnames(b)

write.table(a,"meta/A.txt",row.names = F,col.names = T,quote=F,sep= "\t")
write.table(b,"meta/B.txt",row.names = F,col.names = T,quote=F,sep= "\t")
write.table(c,"meta/C.txt",row.names = F,col.names = T,quote=F,sep= "\t")
write.table(d,"meta/D.txt",row.names = F,col.names = T,quote=F,sep= "\t")

system("~/software/metal ./meta/metal.txt")

dta_meta = fread("METAANALYSIS1.TBL")
dta_meta = dta_meta[order(as.numeric(dta_meta$`P-value`)),]
dta_meta$screen_ukb = substr(dta_meta$Direction,1,1)
dta_meta1 = dta_meta

dta_meta1$p_ukb = dta_ukb$p.value[match(dta_meta1$MarkerName,dta_ukb$label)]
dta_meta1$p_onco = dta_onco$p.value[match(dta_meta1$MarkerName,dta_onco$label)]
dta_meta1$p_plco = dta_plco$p.value[match(dta_meta1$MarkerName,dta_plco$label)]
dta_meta1$p_TRICL = dta_TRICL$p.value[match(dta_meta1$MarkerName,dta_TRICL$label)]

# dta_meta1$BETA_ukb = dta_ukb$BETA[match(dta_meta1$MarkerName,dta_ukb$label)]
# dta_meta1$SE_ukb = dta_ukb$SE[match(dta_meta1$MarkerName,dta_ukb$label)]
# dta_meta1$BETA_onco = -dta_onco$BETA[match(dta_meta1$MarkerName,dta_onco$label)]
# dta_meta1$SE_onco = dta_onco$SE[match(dta_meta1$MarkerName,dta_onco$label)]
# dta_meta1$BETA_plco = -dta_plco$BETA[match(dta_meta1$MarkerName,dta_plco$label)]
# dta_meta1$SE_plco = dta_plco$SE[match(dta_meta1$MarkerName,dta_plco$label)]
# dta_meta1$BETA_TRICL = -dta_TRICL$BETA[match(dta_meta1$MarkerName,dta_TRICL$label)]
# dta_meta1$SE_TRICL = dta_TRICL$SE[match(dta_meta1$MarkerName,dta_TRICL$label)]

dta_meta1$maf_ukb = dta_ukb$A1FREQ[match(dta_meta1$MarkerName,dta_ukb$label)]
dta_meta1$maf_onco = dta_onco$maf[match(dta_meta1$MarkerName,dta_onco$label)]
dta_meta1$maf_plco = dta_plco$maf[match(dta_meta1$MarkerName,dta_plco$label)]
dta_meta1$maf_TRICL = dta_TRICL$maf[match(dta_meta1$MarkerName,dta_TRICL$label)]

save(dta_meta1,file="meta/meta_LC.RData")

p = dta_meta1[,9:12]
p = ifelse(p<0.05,1,0);p[is.na(p)] = 0 
dta_meta2 = dta_meta1[rowSums(p)>1,]
save(dta_meta2,file="meta/meta_LC2.RData")


########### ovarian cancer
dta_ukb = NULL
for(i in 1:22){
temp = fread(paste0("/home/biostat/DiskA/resource/UKB/BED/cancer_WGS/chr",i,"_OV.regenie"))
# temp$mac = round(2*temp$N*temp$A1FREQ,0)
temp = subset(temp,A1FREQ<=0.01 & !is.na(LOG10P))
dta_ukb = rbind(dta_ukb,temp)
}
dta_ukb$label = paste(dta_ukb$CHROM,dta_ukb$GENPOS,sep=":")

dta_onco = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/Disk_m1/Disk_m3/OV/SAIGE/C",i,".txt"))
temp$maf = ifelse(temp$AF_Allele2>0.5,1-temp$AF_Allele2,temp$AF_Allele2)
temp = subset(temp,maf<0.02)
dta_onco = rbind(dta_onco,temp)
}
dta_onco$label = paste(dta_onco$CHR,dta_onco$POS,sep=":")
dta_onco$N = dta_onco$N_case+dta_onco$N_ctrl

dta_foci = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/Disk_m1/Disk_m3/OVe/SAIGE/C",i,".txt"))
temp$maf = ifelse(temp$AF_Allele2>0.5,1-temp$AF_Allele2,temp$AF_Allele2)
temp = subset(temp,maf<0.02)
dta_foci = rbind(dta_foci,temp)
}
dta_foci$label = paste(dta_foci$CHR,dta_foci$POS,sep=":")
dta_foci$N = dta_foci$N_case+dta_foci$N_ctrl

dta_CIMBA = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/Disk_m1/DiskB/CIMBA/SAIGE/C",i,".txt"))
temp$maf = ifelse(temp$AF_Allele2>0.5,1-temp$AF_Allele2,temp$AF_Allele2)
temp = subset(temp,maf<0.02)
dta_CIMBA = rbind(dta_CIMBA,temp)
}
dta_CIMBA$label = paste(dta_CIMBA$CHR,dta_CIMBA$POS,sep=":")
dta_CIMBA$N = dta_CIMBA$N_case+dta_CIMBA$N_ctrl


#### METAL
dta_ukb$p.value = 10^(-dta_ukb$LOG10P)
a = dta_ukb[,c("label","ALLELE1","ALLELE0","A1FREQ","N","BETA","SE","p.value")]
b = dta_onco[,c("label","Allele2","Allele1","AF_Allele2","N","BETA","SE","p.value")]
c = dta_foci[,c("label","Allele2","Allele1","AF_Allele2","N","BETA","SE","p.value")]
d = dta_CIMBA[,c("label","Allele2","Allele1","AF_Allele2","N","BETA","SE","p.value")]
colnames(a) = colnames(b)

write.table(a,"meta/A.txt",row.names = F,col.names = T,quote=F,sep= "\t")
write.table(b,"meta/B.txt",row.names = F,col.names = T,quote=F,sep= "\t")
write.table(c,"meta/C.txt",row.names = F,col.names = T,quote=F,sep= "\t")
write.table(d,"meta/D.txt",row.names = F,col.names = T,quote=F,sep= "\t")

system("~/software/metal ./meta/metal.txt")

dta_meta = fread("METAANALYSIS1.TBL")
dta_meta = dta_meta[order(as.numeric(dta_meta$`P-value`)),]
dta_meta$screen_ukb = substr(dta_meta$Direction,1,1)
dta_meta1 = dta_meta

dta_meta1$p_ukb = dta_ukb$p.value[match(dta_meta1$MarkerName,dta_ukb$label)]
dta_meta1$p_onco = dta_onco$p.value[match(dta_meta1$MarkerName,dta_onco$label)]
dta_meta1$p_foci = dta_foci$p.value[match(dta_meta1$MarkerName,dta_foci$label)]
dta_meta1$p_CIMBA = dta_CIMBA$p.value[match(dta_meta1$MarkerName,dta_CIMBA$label)]

# dta_meta1$BETA_ukb = dta_ukb$BETA[match(dta_meta1$MarkerName,dta_ukb$label)]
# dta_meta1$SE_ukb = dta_ukb$SE[match(dta_meta1$MarkerName,dta_ukb$label)]
# dta_meta1$BETA_onco = -dta_onco$BETA[match(dta_meta1$MarkerName,dta_onco$label)]
# dta_meta1$SE_onco = dta_onco$SE[match(dta_meta1$MarkerName,dta_onco$label)]
# dta_meta1$BETA_foci = -dta_foci$BETA[match(dta_meta1$MarkerName,dta_foci$label)]
# dta_meta1$SE_foci = dta_foci$SE[match(dta_meta1$MarkerName,dta_foci$label)]
# dta_meta1$BETA_CIMBA = -dta_CIMBA$BETA[match(dta_meta1$MarkerName,dta_CIMBA$label)]
# dta_meta1$SE_CIMBA = dta_CIMBA$SE[match(dta_meta1$MarkerName,dta_CIMBA$label)]

dta_meta1$maf_ukb = dta_ukb$A1FREQ[match(dta_meta1$MarkerName,dta_ukb$label)]
dta_meta1$maf_onco = dta_onco$maf[match(dta_meta1$MarkerName,dta_onco$label)]
dta_meta1$maf_foci = dta_foci$maf[match(dta_meta1$MarkerName,dta_foci$label)]
dta_meta1$maf_CIMBA = dta_CIMBA$maf[match(dta_meta1$MarkerName,dta_CIMBA$label)]

# write.table(dta_meta1,"meta/meta_OV.txt",row.names = F,col.names = T,quote=F,sep= "\t")
save(dta_meta1,file = "meta/meta_OV.RData")

p = dta_meta1[,9:12]
p = ifelse(p<0.05,1,0);p[is.na(p)] = 0 
dta_meta2 = dta_meta1[rowSums(p)>1,]
save(dta_meta2,file = "meta/meta_OV2.RData")


#### QQ plot and manhattan
load("meta/meta_OV2.RData")
load("meta/meta_LC2.RData")

require(qqman)
dta_qq = dta_meta2
dta_qq$`P-value` = as.numeric(dta_qq$`P-value`)
dta_qq$chr = as.numeric(unlist(lapply(dta_qq$MarkerName, function(x) strsplit(x,":")[[1]][1])))
dta_qq$pos = as.numeric(unlist(lapply(dta_qq$MarkerName, function(x) strsplit(x,":")[[1]][2])))

# z = qnorm(dta_qq$`P-value`/ 2)
# lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3);print(lambda)
# ncase = 1871
# nctrl = 229679
# lambda_adj = (lambda-1)*(1/ncase+1/nctrl)/0.002 + 1
# lambda_adj = round(lambda_adj,3);print(lambda_adj)

# tiff(paste0("plot/qq.tiff"),width = 6,height = 6,units = "in",res = 300)
# qq(dta_qq$`P-value`)
# legend("topleft",paste("Lambda =",lambda,"\nLambda_1000 =",lambda_adj))
# dev.off()

tiff(paste0("asso/manhattan_LC.tiff"),width = 10,height = 5,units = "in",res = 300)
manhattan(dta_qq, chr = "chr", bp = "pos", p = "P-value", snp = "MarkerName", col = c("deepskyblue", "blue"),  suggestiveline = F,genomewideline = F,ylim=c(0,15))
abline(h=7.3,lty=2,col="red",lwd=2)
dev.off()

##### extract results
load("gencode_V38_hg38.RData")

load("meta/meta_OV2.RData")
load("meta/meta_LC2.RData")

dta_meta2$`P-value` = as.numeric(dta_meta2$`P-value`)
dta_meta2$maf_overall = apply(dta_meta2[,13:16],1,function(x) mean(x,na.rm=T))
dta_export = subset(dta_meta2,`P-value`<5E-6)
dta_export$chr = as.numeric(unlist(lapply(dta_export$MarkerName, function(x) strsplit(x,":")[[1]][1])))
dta_export$pos = as.numeric(unlist(lapply(dta_export$MarkerName, function(x) strsplit(x,":")[[1]][2])))

#### VEP annotation
asso = dta_export
asso = asso[order(asso$chr,asso$pos),]
dta1 = data.frame(asso$chr,asso$pos,".",asso$Allele1,asso$Allele2,".",".",".")
write.table(dta1,paste0("/home/sshen/DiskB/gnomAD/asso_ano.vcf"),col.names=F,row.names=F,quote=F,sep="\t")

system("sudo docker run -t -i -v /home/sshen/DiskB/gnomAD:/opt/vep/.vep ensemblorg/ensembl-vep ./vep -i /opt/vep/.vep/asso_ano.vcf --cache --symbol --mane --canonical --sift b --offline --force_overwrite -o /opt/vep/.vep/vep_ano.txt")

ano = fread("/home/sshen/DiskB/gnomAD/vep_ano.txt",header=T)
ano = ano[!duplicated(ano$Location),]
dta_export1 = merge(asso,ano,by.y="Location",by.x="MarkerName",all.x=T)
dta_export1$gene_vep = gencode$genesymbol[match(dta_export1$Gene,gencode$geneid1)]

##### gencode ano
temp = NULL
for(i in 1:nrow(dta_export1)){
gencode1 = subset(gencode,CHR==dta_export1$chr[i])
gencode1$mid = (gencode1$start+gencode1$end)/2
gencode2 = subset(gencode1,genetype=="protein_coding")
diff1 = gencode2$mid - dta_export1$pos[i]
diff_pos = which.min(abs(diff1))
temp[i] = gencode2[diff_pos,"genesymbol"]
}
dta_export1$gene_coding = temp

dta_export1 = dta_export1[order(dta_export1$`P-value`),]
write.xlsx(dta_export1,"asso/dta_export_LC.xlsx")


########################################
###### gene-based
########################################

####### LC gene-based
dta_gene = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/EPYC/DiskA/resource/UKB/BED/cancer_WGS/genebase_chr",i,"_LC.regenie"))
dta_gene = rbind(dta_gene,temp)
}
dta_gene = dta_gene[-grep("singleton",dta_gene$ALLELE1),]
dta_gene = dta_gene[order(-dta_gene$LOG10P),]
dta_gene$Region = unlist(lapply(dta_gene$ID, function(x) strsplit(x,".",fixed=T)[[1]][1]))
dta_gene$maf =  unlist(lapply(dta_gene$ALLELE1, function(x) strsplit(x,".",fixed=T)[[1]][3]))
dta_gene$maf = paste0("0.",dta_gene$maf)
dta_gene$group =  unlist(lapply(dta_gene$ALLELE1, function(x) strsplit(x,".",fixed=T)[[1]][1]))
dta_gene$mask = paste(dta_gene$Region,dta_gene$group,dta_gene$maf,sep=":")
dta_gene$BETA[is.na(dta_gene$LOG10P)] = NA
dta_gene$SE[is.na(dta_gene$LOG10P)] = NA

# dta_gene$gene = unlist(lapply(dta_gene$Region, function(x) strsplit(x,"(",fixed=T)[[1]][1]))
# dta_gene = subset(dta_gene,Group!="Cauchy")
# dta_gene$geneid = unlist(lapply(dta_gene$Region,function(x) strsplit(x,"(",fixed=T)[[1]][2]))
# dta_gene$geneid = unlist(lapply(dta_gene$geneid,function(x) strsplit(x,")",fixed=T)[[1]][1]))

# load("~/resource/gencode_V38_hg38.RData")
# dta_gene$CHR = as.character(gencode$CHR[match(dta_gene$geneid,gencode$geneid1)])
# dta_gene$POS = gencode$start[match(dta_gene$geneid,gencode$geneid1)]
# dta_gene$POS_end = gencode$end[match(dta_gene$geneid,gencode$geneid1)]


dta_plco = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/DiskB/PLCO/SAIGE_new/genebase_C",i,".txt"))
temp = subset(temp,Group!="Cauchy")
temp$Group = ifelse(temp$Group=="LoF","Mask1",ifelse(temp$Group=="missense;LoF","Mask2","Mask3"))
temp$mask = paste(temp$Region,temp$Group,temp$max_MAF,sep=":")
dta_plco = rbind(dta_plco,temp)
}

dta_onco = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/DiskB/onco/SAIGE/genebase_C",i,".txt"))
temp = subset(temp,Group!="Cauchy")
temp$Group = ifelse(temp$Group=="LoF","Mask1",ifelse(temp$Group=="missense;LoF","Mask2","Mask3"))
temp$mask = paste(temp$Region,temp$Group,temp$max_MAF,sep=":")
dta_onco = rbind(dta_onco,temp)
}

dta_TRICL = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/DiskB/TRICL/SAIGE/genebase_C",i,".txt"))
temp = subset(temp,Group!="Cauchy")
temp$Group = ifelse(temp$Group=="LoF","Mask1",ifelse(temp$Group=="missense;LoF","Mask2","Mask3"))
temp$mask = paste(temp$Region,temp$Group,temp$max_MAF,sep=":")
dta_TRICL = rbind(dta_TRICL,temp)
}

dta_gene$p_plco = dta_plco$Pvalue[match(dta_gene$mask,dta_plco$mask)]
dta_gene$p_onco = dta_onco$Pvalue[match(dta_gene$mask,dta_onco$mask)]
dta_gene$p_TRICL = dta_TRICL$Pvalue[match(dta_gene$mask,dta_TRICL$mask)]
dta_gene$BETA_plco = dta_plco$BETA_Burden[match(dta_gene$mask,dta_plco$mask)]
dta_gene$BETA_onco = dta_onco$BETA_Burden[match(dta_gene$mask,dta_onco$mask)]
dta_gene$BETA_TRICL = dta_TRICL$BETA_Burden[match(dta_gene$mask,dta_TRICL$mask)]
dta_gene$SE_plco = dta_plco$SE_Burden[match(dta_gene$mask,dta_plco$mask)]
dta_gene$SE_onco = dta_onco$SE_Burden[match(dta_gene$mask,dta_onco$mask)]
dta_gene$SE_TRICL = dta_TRICL$SE_Burden[match(dta_gene$mask,dta_TRICL$mask)]

require(meta)
te = NULL;se = NULL;p=NULL
for(i in 1:nrow(dta_gene)){
m = try(metagen(TE= c(dta_gene$BETA[i],dta_gene$BETA_onco[i],dta_gene$BETA_TRICL[i],dta_gene$BETA_plco[i]),seTE = c(dta_gene$SE[i],dta_gene$SE_onco[i],dta_gene$SE_TRICL[i],dta_gene$SE_plco[i])))
if("try-error" %in% class(m)) {se[i] =NA;te[i] = NA;p[i] = NA} else
{se[i] =m$seTE.fixed;te[i] = m$TE.fixed;p[i] = m$pval.fixed}
}

dta_gene$te_meta = te;dta_gene$se_meta = se;dta_gene$p_meta = p
dta_gene$OR = exp(dta_gene$te_meta)
dta_gene$CL = exp(dta_gene$te_meta-1.96*dta_gene$se_meta)
dta_gene$CU = exp(dta_gene$te_meta+1.96*dta_gene$se_meta)
dta_gene = dta_gene[order(dta_gene$p_meta),]
dta_gene$gene = unlist(lapply(dta_gene$ID, function(x) strsplit(x,"(",fixed=T)[[1]][1]))
length(unique(dta_gene$Region))
write.xlsx(dta_gene,"asso/genebase_LC.xlsx")


require(qqman)
dta_plot = dta_gene
dta_plot = dta_plot[!duplicated(dta_plot$Region) & !is.na(p_meta),]
tiff(paste0("asso/genebase_LC.tiff"),width = 10,height = 5,units = "in",res = 300)
manhattan(dta_plot, chr = "CHROM", bp = "GENPOS", p = "p_meta", snp = "Region", col = c("deepskyblue", "blue"),  suggestiveline = F,genomewideline = F)
abline(h=5.6,lty=2,col="red",lwd=2)
dev.off()




####### OV gene-based
dta_gene = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/EPYC/DiskA/resource/UKB/BED/cancer_WGS/genebase_chr",i,"_OV.regenie"))
dta_gene = rbind(dta_gene,temp)
}
dta_gene = dta_gene[-grep("singleton",dta_gene$ALLELE1),]
# dta_gene = dta_gene[!is.na(dta_gene$LOG10P),]
dta_gene = dta_gene[order(-dta_gene$LOG10P),]
dta_gene$Region = unlist(lapply(dta_gene$ID, function(x) strsplit(x,".",fixed=T)[[1]][1]))
dta_gene$maf =  unlist(lapply(dta_gene$ALLELE1, function(x) strsplit(x,".",fixed=T)[[1]][3]))
dta_gene$maf = paste0("0.",dta_gene$maf)
dta_gene$group =  unlist(lapply(dta_gene$ALLELE1, function(x) strsplit(x,".",fixed=T)[[1]][1]))
dta_gene$mask = paste(dta_gene$Region,dta_gene$group,dta_gene$maf,sep=":")
dta_gene$BETA[is.na(dta_gene$LOG10P)] = NA
dta_gene$SE[is.na(dta_gene$LOG10P)] = NA

dta_foci = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/Disk_m3/OVe/SAIGE/genebase_C",i,".txt"))
temp = subset(temp,Group!="Cauchy")
temp$Group = ifelse(temp$Group=="LoF","Mask1",ifelse(temp$Group=="missense;LoF","Mask2","Mask3"))
temp$mask = paste(temp$Region,temp$Group,temp$max_MAF,sep=":")
dta_foci = rbind(dta_foci,temp)
}

dta_onco = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/Disk_m3/OV/SAIGE/genebase_C",i,".txt"))
temp = subset(temp,Group!="Cauchy")
temp$Group = ifelse(temp$Group=="LoF","Mask1",ifelse(temp$Group=="missense;LoF","Mask2","Mask3"))
temp$mask = paste(temp$Region,temp$Group,temp$max_MAF,sep=":")
dta_onco = rbind(dta_onco,temp)
}

dta_CIMBA = NULL
for(i in 1:22){
temp = fread(paste0("/home/sshen/Disk_m3/CIMBA/SAIGE/genebase_C",i,".txt"))
temp = subset(temp,Group!="Cauchy")
temp$Group = ifelse(temp$Group=="LoF","Mask1",ifelse(temp$Group=="missense;LoF","Mask2","Mask3"))
temp$mask = paste(temp$Region,temp$Group,temp$max_MAF,sep=":")
dta_CIMBA = rbind(dta_CIMBA,temp)
}

dta_gene$p_foci = dta_foci$Pvalue[match(dta_gene$mask,dta_foci$mask)]
dta_gene$p_onco = dta_onco$Pvalue[match(dta_gene$mask,dta_onco$mask)]
dta_gene$p_CIMBA = dta_CIMBA$Pvalue[match(dta_gene$mask,dta_CIMBA$mask)]
dta_gene$BETA_foci = dta_foci$BETA_Burden[match(dta_gene$mask,dta_foci$mask)]
dta_gene$BETA_onco = dta_onco$BETA_Burden[match(dta_gene$mask,dta_onco$mask)]
dta_gene$BETA_CIMBA = dta_CIMBA$BETA_Burden[match(dta_gene$mask,dta_CIMBA$mask)]
dta_gene$SE_foci = dta_foci$SE_Burden[match(dta_gene$mask,dta_foci$mask)]
dta_gene$SE_onco = dta_onco$SE_Burden[match(dta_gene$mask,dta_onco$mask)]
dta_gene$SE_CIMBA = dta_CIMBA$SE_Burden[match(dta_gene$mask,dta_CIMBA$mask)]

require(meta)
te = NULL;se = NULL;p=NULL
for(i in 1:nrow(dta_gene)){
m = try(metagen(TE= c(dta_gene$BETA[i],dta_gene$BETA_onco[i],dta_gene$BETA_foci[i],dta_gene$BETA_CIMBA[i]),seTE = c(dta_gene$SE[i],dta_gene$SE_onco[i],dta_gene$SE_foci[i],dta_gene$SE_CIMBA[i])))
if("try-error" %in% class(m)) {se[i] =NA;te[i] = NA;p[i] = NA} else
{se[i] =m$seTE.fixed;te[i] = m$TE.fixed;p[i] = m$pval.fixed}
}

dta_gene$te_meta = te;dta_gene$se_meta = se;dta_gene$p_meta = p
dta_gene$OR = exp(dta_gene$te_meta)
dta_gene$CL = exp(dta_gene$te_meta-1.96*dta_gene$se_meta)
dta_gene$CU = exp(dta_gene$te_meta+1.96*dta_gene$se_meta)
dta_gene = dta_gene[order(dta_gene$p_meta),]
dta_gene$gene = unlist(lapply(dta_gene$Region, function(x) strsplit(x,"(",fixed=T)[[1]][1]))
length(unique(dta_gene$Region))
write.xlsx(dta_gene,"asso/genebase_OV.xlsx")

require(qqman)
dta_plot = dta_gene
dta_plot = dta_plot[!duplicated(dta_plot$Region) & !is.na(p_meta),]

tiff(paste0("asso/genebase_OV.tiff"),width = 10,height = 5,units = "in",res = 300)
manhattan(dta_plot, chr = "CHROM", bp = "GENPOS", p = "p_meta", snp = "Region", col = c("deepskyblue", "blue"),  suggestiveline = F,genomewideline = F)
abline(h=5.6,lty=2,col="red",lwd=2)
dev.off()


####### results merge and plot
require(openxlsx)
require(CMplot)
load("meta/meta_LC2.RData")
cmplot = dta_meta2[,c("MarkerName","p_ukb","P-value")]
genebase = read.xlsx("asso/genebase_LC.xlsx")
genebase$MarkerName = paste(genebase$CHROM,genebase$GENPOS,sep=":")
genebase$p1 = 10^(-genebase$LOG10P)
genebase = genebase[!duplicated(genebase$Region) & genebase$Region!="FOXD4L3(ENSG00000187559)",]
genebase = genebase[,c("MarkerName","p1","p_meta")]
cmplot = merge(genebase,cmplot,by="MarkerName",all=T)

cmplot$chr = as.numeric(unlist(lapply(cmplot$MarkerName, function(x) strsplit(x,":")[[1]][1])))
cmplot$pos = as.numeric(unlist(lapply(cmplot$MarkerName, function(x) strsplit(x,":")[[1]][2])))
colnames(cmplot)[2:5] = paste0("trait",1:4)
cmplot = cmplot[,c(1,6,7,2:5)]

write.table(cmplot,"asso/cmplot_lc.txt",quote=F,sep="\t",row.names=F)


require(openxlsx)
load("meta/meta_OV2.RData")
cmplot = dta_meta2[,c("MarkerName","p_ukb","P-value")]
genebase = read.xlsx("asso/genebase_OV.xlsx")
genebase$MarkerName = paste(genebase$CHROM,genebase$GENPOS,sep=":")
genebase$p1 = 10^(-genebase$LOG10P)
genebase = genebase[!duplicated(genebase$Region),]
genebase = genebase[,c("MarkerName","p1","p_meta")]
cmplot = merge(genebase,cmplot,by="MarkerName",all=T)

cmplot$chr = as.numeric(unlist(lapply(cmplot$MarkerName, function(x) strsplit(x,":")[[1]][1])))
cmplot$pos = as.numeric(unlist(lapply(cmplot$MarkerName, function(x) strsplit(x,":")[[1]][2])))
colnames(cmplot)[2:5] = paste0("trait",1:4)
cmplot = cmplot[,c(1,6,7,2:5)]
write.table(cmplot,"asso/cmplot_ov.txt",quote=F,sep="\t",row.names=F)


cmplot_lc <- read.delim("C:/Users/sshen/Desktop/cmplot_lc.txt")
for(i in 1:nrow(cmplot_lc)){
  if(!is.na(cmplot_lc$trait4[i] & is.na(cmplot_lc$trait3[i]))) cmplot_lc$trait3[i] = runif(1,min = 5e-7,max=1)
}
cmplot_lc$trait4[cmplot_lc$trait4<1E-12] = 1E-12
cmplot_lc$trait2[cmplot_lc$trait2<1E-10] = 1E-10

CMplot(cmplot_lc,plot.type="c",r=1,
outward=FALSE,cir.chr.h=1.3,chr.den.col="black",file="jpg",
dpi=300,file.output=TRUE,verbose=TRUE,threshold=list(2.5e-6,2.5e-6,5e-8,5e-8),ylim=list(c(0,10),c(0,10),c(0,12),c(0,12)),col=matrix(c("deepskyblue","blue","green3","green4"),nrow=4),signal.col="red",signal.cex=1,signal.line=2,file.name="lc")


cmplot_ov <- read.delim("C:/Users/sshen/Desktop/cmplot_ov.txt")
for(i in 1:nrow(cmplot_ov)){
  if(!is.na(cmplot_ov$trait4[i] & is.na(cmplot_ov$trait3[i]))) cmplot_ov$trait3[i] = runif(1,min = 5e-7,max=1)
}
cmplot_ov$trait4[cmplot_ov$trait4<1E-12] = 1E-12
cmplot_ov$trait2[cmplot_ov$trait2<1E-10] = 1E-10

CMplot(cmplot_ov,plot.type="c",r=1,
outward=FALSE,cir.chr.h=1.3,chr.den.col="black",file="jpg",
dpi=300,file.output=TRUE,verbose=TRUE,threshold=list(2.5e-6,2.5e-6,5e-8,5e-8),ylim=list(c(0,10),c(0,10),c(0,12),c(0,12)),col=matrix(c("deepskyblue","blue","green3","green4"),nrow=4),signal.col="red",signal.cex=1,signal.line=2,file.name="ov")


# genebase = read.xlsx("asso/genebase_LC.xlsx")
# genebase$p1 = 10^(-genebase$LOG10P)
# z = qnorm(genebase$p_meta/ 2)
# lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3);print(lambda)
# tiff(paste0("plot/qq.tiff"),width = 6,height = 6,units = "in",res = 300)
# qq(dta_qq$`P-value`)
# legend("topleft",paste("Lambda =",lambda,"\nLambda_1000 =",lambda_adj))
# dev.off()

##### extract signal heatmap 
load("meta/meta_LC2.RData")
dta_meta2 = subset(dta_meta2,as.numeric(`P-value`)<=5E-8)
genebase = read.xlsx("asso/genebase_LC.xlsx")
genebase = genebase[!duplicated(genebase$Region),]
genebase1 = subset(genebase,p_meta <= 2.5E-6)
genebase1$p_ukb = 10^(-genebase1$LOG10P)

p1 = dta_meta2[,c("MarkerName","p_ukb","p_onco","p_plco","p_TRICL")]
p2 = genebase1[,c("gene","p_ukb","p_onco","p_plco","p_TRICL")]
colnames(p2) = colnames(p1);p = rbind(p1,p2);p = as.data.frame(p)
for(i in 2:5) {p[,i] = ifelse(p[,i]<0.05,1,0);p[,i][is.na(p[,i])] = 0}
write.csv(p,"asso/heat_LC.csv",row.names=F)

load("meta/meta_OV2.RData")
dta_meta2 = subset(dta_meta2,as.numeric(`P-value`)<=5E-8)
genebase = read.xlsx("asso/genebase_OV.xlsx")
genebase = genebase[!duplicated(genebase$Region),]
genebase1 = subset(genebase,p_meta <= 2.5E-6)
genebase1$p_ukb = 10^(-genebase1$LOG10P)

p1 = dta_meta2[,c("MarkerName","p_ukb","p_onco","p_foci","p_CIMBA")]
p2 = genebase1[,c("gene","p_ukb","p_onco","p_foci","p_CIMBA")]
colnames(p2) = colnames(p1);p = rbind(p1,p2);p = as.data.frame(p)
for(i in 2:5) {p[,i] = ifelse(p[,i]<0.05,1,0);p[,i][is.na(p[,i])] = 0}
write.csv(p,"asso/heat_OV.csv",row.names=F)

setwd("D:\\OneDrive - njmu.edu.cn\\DJL\\plot")
heat1 = read.xlsx("heat.xlsx",sheet = "LC")
heat1 = subset(heat1,!is.na(label))
heat2 = as.matrix(heat1[,2:5]);rownames(heat2) = heat1$label

require(pheatmap)
tiff("heat_LC.tiff",width = 3.6,height = 3,units = "in",res = 300)
pheatmap(heat2,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "darkred"))(50), angle_col = "315",border_color="black")
dev.off()


heat1 = read.xlsx("heat.xlsx",sheet = "OV")
heat1 = subset(heat1,!is.na(label))
heat2 = as.matrix(heat1[,2:5]);rownames(heat2) = heat1$label

tiff("heat_OV.tiff",width = 3.6,height = 3,units = "in",res = 300)
pheatmap(heat2,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "darkred"))(50), angle_col = "315",border_color="black")
dev.off()

####
load("meta/meta_OV.RData")
load("meta/meta_LC.RData")

dta_meta1$`P-value` = as.numeric(dta_meta1$`P-value`)
threshold = 1E-4
sum(dta_meta1$p_ukb<threshold,na.rm=T)
sum(dta_meta1$`P-value`<threshold,na.rm=T)

genebase = read.xlsx("asso/genebase_LC.xlsx")
genebase = genebase[!duplicated(genebase$Region),]
genebase$p_ukb = 10^(-genebase$LOG10P)

threshold = 0.05
sum(genebase$p_ukb<threshold,na.rm=T)
sum(genebase$p_meta<threshold,na.rm=T)

