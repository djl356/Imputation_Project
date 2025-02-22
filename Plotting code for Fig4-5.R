############# Plotting code of Figure 4 ###########
#### Fig 4A ####
library(gtools)
library(ggplot2)
library(ggtext)
library(grid)
library(ggpubr)

## The following three files are obtained after processing the results of association analysis
## The code of processing association result see the 'Association tests for rare variants.R'

setwd("/home/biostat/DiskA/resource/UKB/BED")

hrc = fread("biochemistry_HRC/result/0.01summary/single0.01_sum.csv")
topmed = fread("biochemistry_topmed/result/0.01summary/single0.01_sum.csv")
wgs = fread("biochemistry_WGS/result/0.01summary/single0.01_sum.csv")

hrc$panel = "HRC+UK10K"
topmed$panel = "TOPMed"
wgs$panel = "WGS"
panel_pheno = rbind(topmed,hrc,wgs)

# get the label file
lab_file = fread("biochemistry_datafields.csv",header = T,sep = ",")
labels = lab_file[,c("Description","label","classification")]
labels$label = tolower(labels$label)
labels = labels[mixedorder(labels$label),]

pheno_level = labels$label
pheno_lab = labels$Description

panel_pheno$type = factor(panel_pheno$type,levels = pheno_level,labels = pheno_lab)
panel_pheno$panel = factor(panel_pheno$panel,levels = c("TOPMed","HRC+UK10K","WGS"))

options(scipen=200)

panel_col = ggplot(panel_pheno, aes(fill=panel, y=sum, x=type)) + 
  geom_col(position="identity",width = 0.8) +
  labs(x = "",y = "Numbers of significant 
rare variants") +
  scale_fill_manual(name="",
                    labels = c("TOPMed","HRC+UK10K","WGS"),
                    values=c("salmon","darkturquoise","royalblue3")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 19,color = "black",angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size = 19,color = "black"),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.5, colour="black", linetype="solid"),
        axis.ticks.length = unit(2, "mm"),
        axis.title.x = element_text(size = 21,vjust=-1,hjust=0.5),
        axis.title.y = element_text(size = 21,vjust=1.5,hjust = 0.5),
        legend.text = element_text(size = 21),
        legend.key.size = unit(24,"pt"),
        legend.spacing.x = unit(0.5,"cm"),
        legend.position = "top")

#### Fig 4B (take single-variant tests as a example) ####
## The following three files are obtained after processing the results of association analysis
## The code of processing association result see the 'Association tests for rare variants.R'
## Note:replace the 'single0.01_sum.csv' file in each path with 'gene_unique_sum.csv' files to draw plot for gene-based tests

hrc = fread("biochemistry_HRC/result/0.01summary/single0.01_sum.csv")
topmed = fread("biochemistry_topmed/result/0.01summary/single0.01_sum.csv")
wgs = fread("biochemistry_WGS/result/0.01summary/single0.01_sum.csv")

hrc$panel = "HRC+UK10K"
topmed$panel = "TOPMed"
wgs$panel = "WGS"
panel_pheno = rbind(topmed,hrc,wgs)

library(scales)
library(dplyr)
single0.01 = merge(panel_pheno,labels,by.x="type",by.y="label")
single0.01 = single0.01[,c("Description","sum","panel")]

single0.01 = spread(single0.01,panel,sum)

## calculate the improvement ratio of each filled data relative to WGS data

single0.01$hrc_per = (single0.01$`HRC+UK10K`-single0.01$WGS)/single0.01$WGS
single0.01$top_per = (single0.01$`TOPMed`-single0.01$WGS)/single0.01$WGS

per_data = single0.01[,c("Description","hrc_per","top_per")]

## negative ratio are recorded as 0 (equivalent to only calculating the improvement ratio)
per_data$hrc_per[per_data$hrc_per<0] = 0 
per_data$top_per[per_data$top_per<0] = 0 

panel = c("HRC+UK10K","TOPMed")
hrc_per = mean(per_data$hrc_per)
top_per = mean(per_data$top_per)

per_mean = c(hrc_per,top_per)
per_data = data.table(panel,per_mean)

per_data$panel = factor(per_data$panel,levels = c("TOPMed","HRC+UK10K"))
per_data = per_data[order(per_data$panel),]
per_data$text = percent(per_data$per_mean,accuracy=0.01)

library(ggplot2)
library(ggtext)
library(grid)
library(ggpubr)

options(scipen=200)

per_col = ggplot(per_data, aes(fill=panel, y=per_mean, x=panel)) + 
  geom_col(width=0.5) +
  labs(x = "Reference panel",y = "Percentage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(name="",
                    labels = c("TOPMed","HRC+UK10K"),
                    values=c("salmon","darkturquoise")) +
  geom_text(aes(label=per_data$text),size=8,color = "black",vjust=-0.5) +
  geom_hline(yintercept = 1,size = 1,
             linetype="dashed", color="royalblue3") +
  annotate('text',label = "WGS",
           x = 0.56,y = 1.03,
           size = 7,color = "royalblue3") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 22,color = "black"),
        axis.text.y = element_text(size = 22,color = "black"),
        axis.ticks.x=element_line(size=0.5, colour="black", linetype="solid"),
        axis.ticks.y=element_line(size=0.5, colour="black", linetype="solid"),
        axis.ticks.length = unit(2, "mm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(16,"pt"),
        legend.spacing.x = unit(0.5,"cm"),
        legend.position = "none")

#### Fig 4C (take HRC+UK10K-imputed data as example)####
## Note:replace the 'hrc_sactter.csv' file in each path with 'topmed_scatter.csv' files to draw plot for TOPMed-imputed data
scatter_data = fread("/home/biostat/DiskB/jldai/chisq/hrc_scatter.csv")

cor_index = round(cor(scatter_data$CHISQ_HRC,scatter_data$CHISQ,method = "pearson"),2)

library(ggplot2)
library(ggtext)
library(ggpubr)

chisq_plot = ggplot(scatter_data, aes(x=CHISQ_HRC, y=CHISQ)) + 
  geom_point(shape = 21, size = 3.5,color = "darkturquoise") +
  labs(x=expression("Genotype imputed from HRC into n=150K"),y= "Gneotype from n=150K WGS") +
  geom_abline(slope = 1,intercept = 0,lty="solid",color = "gray") +
  annotate('text',label = paste0("Pearson's r = ",cor_index),
           x = 275,y = 100,
           size = 8,color = "black") +
  scale_x_continuous(limits = c(0,350),breaks = c(0,50,100,150,200,250,300,350)) +
  scale_y_continuous(?limits = c(0,350),breaks = c(0,50,100,150,200,250,300,350)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 22,color = "black"),
        axis.text.y = element_text(size = 22,color = "black"),
        axis.ticks.x=element_line(size=0.5, colour="black", linetype="solid"),
        axis.ticks.y=element_line(size=0.5, colour="black", linetype="solid"),
        axis.ticks.length = unit(2, "mm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
############# Plotting code of Figure 5 ###########
#### Manhattan Plot (single-variants tests) ####
library(ggplot2)
library(ggtext)
library(ggpubr)

## Note:should change the file for different dataset
## Here take the HRC+UK10K-imputed data as a example
plot_data = fread("/home/biostat/DiskB/jldai/common+cancer/hrc/hrc_plot_data.csv")
pheno = c("E10","E11","E66","F32","I10",
          "I25","I50","J44","J45","K80",
          "BLCA","BRCA","NHL","PRAD","SKCM")
# change some labels
plot_data$label[plot_data$label=="Non-insulin-dependent diabetes mell"] = "Non-insulin-dependent diabetes mellitus"
plot_data$label[plot_data$label=="Essential (primary) hypertension"] = "Hypertension"
plot_data$label[plot_data$label=="Other chronic obstructive pulmonary"] = "COPD"
plot_data$label[plot_data$label=="Bladder Urothelial Carcinoma"] = "Bladder carcinoma"
plot_data$label[plot_data$label=="Breast invasive carcinoma"] = "Breast carcinoma"
plot_data$label[plot_data$label=="Prostate adenocarcinoma"] = "Prostate cancer"
plot_data$label[plot_data$label=="Skin Cutaneous Melanoma"] = "Melanoma"


# convert to character
plot_data$chrom = as.character(plot_data$chrom)
col_data = plot_data[,c("chrom","col_now")]
col_data = unique(col_data)

plot_data$CANCER = factor(plot_data$CANCER,levels = pheno)
plot_data$chrom = factor(plot_data$chrom,levels = as.character(c(1:330)))
plot_data = plot_data[order(plot_data$CANCER),]

color_data = unique(plot_data[,c("chrom","col_now")])

# set threshold and suggestive line
threshold = -log10(5e-08)
threshold_suggest = -log10(5e-06)

manhattan_plot = ggplot(plot_data,aes(x=chrom,y=LOG10P)) +
  geom_jitter(aes(color=chrom,shape=chrom),size=4.5) +
  scale_shape_manual(values = rep(20,times = 330)) +
  scale_color_manual(values = color_data$col_now)+
  geom_hline(yintercept = threshold,size = 2,
             linetype="dashed", color="red") +
  annotate('text',label = expression("P<5x10"^-8),
           x = 24,y = 8.5,
           size = 9.5,color = "red") +
  geom_hline(yintercept = threshold_suggest,size = 2,
             linetype="dashed", color="gray58") +
  scale_y_continuous(expand = c(0, 0),limits = c(0,25)) +
  scale_x_discrete(expand = expansion(mult = c(0.01,0)),labels = unique(plot_data$label),
                   breaks = as.character(seq(11,330,22))) +
  labs(title = "",x="",y=expression("-log"[10]*"(Pvalue)")) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 24),
        axis.line.y= element_line(size = 1,colour = "black"),
        axis.ticks.y=element_line(size=1, colour="black", linetype="solid"),
        axis.ticks.length.y = unit(2.5, "mm"),
        axis.ticks.x= element_blank(),
        axis.title.x = element_text(size = 24,vjust=-0.5,hjust=0.5),
        axis.title.y = element_text(size = 34,vjust=1.5,hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=32,color="black"),
        legend.position = "none")

#### Manhattan Plot (gene-based tests) ####
library(ggplot2)
library(ggtext)
library(ggrepel)
library(ggpubr)

# creat 15 color (same as the single-variant tests)
col = c("#FB8072","#F4CAE4","cadetblue3","#FDC086","#E7298A","#BC80BD",
        "#8DD3C7","goldenrod1","#4DAF4A","#A6761D","#7570B3","#E78AC3",
        "#B3CDE3","slateblue1","#66A61E")

## Note:should change the file for different dataset
## Here take the HRC+UK10K-imputed data as a example

genebase_plot_data = fread("/home/biostat/DiskB/jldai/common+cancer/hrc/hrc_genebase_plot_data.csv")
genebase_plot_data$CANCER = factor(genebase_plot_data$CANCER,levels = c("E10","E11","E66","F32","I10",
                                                                        "I25","I50","J44","J45","K80",
                                                                        "BLCA","BRCA","NHL","PRAD","SKCM"))
genebase_plot_dataYES = genebase_plot_data[genebase_plot_data$is_annotate=="yes" & genebase_plot_data$UPPER=="0.01",]
genebase_plot_dataLEFT = genebase_plot_data[genebase_plot_data$is_annotate=="no",]

genebase_plot_data = rbind(genebase_plot_dataYES,genebase_plot_dataLEFT)

genebase_plot_data = genebase_plot_data %>%
  arrange(desc(LOG10P)) %>%
  group_by(CANCER,GENE) %>%
  filter(row_number()==1)

genebase_plot_data = genebase_plot_data[order(genebase_plot_data$CANCER),]

manhattan_plot = ggplot(genebase_plot_data,aes(x=CANCER,y=LOG10P)) +
  geom_jitter(aes(color=CANCER,shape=CANCER),size=4.5) +
  scale_shape_manual(values = rep(17,times = 15)) +
  geom_hline(yintercept = -log10(2.5e-06),size = 2,
             linetype="dashed", color="red") +
  annotate('text',label = expression("P<2.5x10"^-6),
           x = 13.5,y = 6.4,
           size = 9.5,color = "red") +
  geom_label_repel(data=subset(genebase_plot_data,is_annotate=="yes"),aes(label=GENE,color=CANCER),
                   size=8,min.segment.length = Inf,
                   box.padding=unit(1, "lines"), point.padding=unit(0.2, "lines"),
                   max.overlaps = getOption("ggrepel.max.overlaps", 50)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,20),
                     breaks = c(0,5,10,15,20)) +
  scale_x_discrete(labels = unique(genebase_plot_data$label)) +
  scale_color_manual(values = col) +
  labs(title="",x="",
       y=expression("-log"[10]*"(Pvalue)")) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.line.y= element_line(size = 1,colour = "black"),
        axis.ticks.y=element_line(size=1, colour="black", linetype="solid"),
        axis.ticks.length.y = unit(2.5, "mm"),
        axis.ticks.x=element_blank(),
        axis.title.x = element_text(size = 24,vjust=-0.5,hjust=0.5),
        axis.title.y = element_text(size = 34,vjust=1.5,hjust = 0.5),
        axis.text.y = element_text(size=32,color="black"),
        axis.text.x = element_blank(),
        legend.position = "none")