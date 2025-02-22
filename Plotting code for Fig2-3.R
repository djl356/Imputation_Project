############# Plotting code of Figure 2 ###########
#### Fig 2B ####
setwd("/home/biostat/DiskA/resource/UKB/BED/total_result/MAF_stat")
hrc = fread("HRC/hrc_stat4.csv")
topmed = fread("TopMed/topmed_stat4.csv")
wgs = fread("WGS/wgs_stat4.csv")

# create panel name
hrc$panel = "HRC+UK10K"
topmed$panel = "TOPMed"
wgs$panel = "WGS"
plot_data = rbind(hrc,topmed,wgs)
plot_data = plot_data[type!="Total",]

# create data labels
plot_data$count_chr = as.character(plot_data$count)
plot_data$text = paste0(sub('......$','',plot_data$count_chr),"M")

plot_data$type[plot_data$type=="singleton"] = "Singleton"
plot_data$type[plot_data$type=="doubleton"] = "Doubleton"

plot_data$type = factor(plot_data$type,levels = c("Singleton","Doubleton","Ultra-rare variation","Rare variation","Common variation"),
                        labels = c("Singleton","Doubleton","Ultra-rare variant","Rare variant","Common variant"))
plot_data$panel = factor(plot_data$panel,levels = c("HRC+UK10K","TOPMed","WGS"))

# required packages
library(ggplot2)
library(ggtext)
library(grid)
library(ggpubr)

snv_type = ggplot(plot_data, aes(fill = panel, y = count, x = type)) +
  geom_col(width = 0.9,position = "dodge") +
  labs( title = "",
        x = "Variant type",
        y = "Numbers of variant") +
  scale_fill_manual(labels = c("HRC+UK10K","TOPMed","WGS"),
                    values = c("darkturquoise","salmon","royalblue3")) +
  scale_y_continuous(name = "Numbers of variant",
                     breaks = c(0,50000000,100000000,150000000,200000000,250000000),
                     labels = c("0M","50M","100M","150M","200M","250M"))+  
  geom_text(aes(label=plot_data$text),size=4,color = "grey35",vjust=-0.5,position = position_dodge(0.8)) +
  theme_bw() +
  theme(plot.margin=unit(rep(1,4),'cm'),
        line = element_line(size = 1,linetype = 2),
        panel.grid = element_blank(),   			## 去除背景所有网格线	
        panel.border = element_blank(),
        plot.title = element_text(size = 15,color = "black",vjust = 2),
        axis.line = element_line(size = 0.2,linetype = "solid",color = "black"),
        axis.text = element_text(face="plain",color="black",size=19),
        axis.text.x = element_text(angle = 30,vjust = 0.80,hjust = 0.75),
        axis.ticks=element_line(size=0.5, colour="black", linetype="solid"),
        axis.ticks.length = unit(1, "mm"),
        axis.title.x = element_text(size=21,vjust=-2,hjust=0.5),
        axis.title.y = element_text(size=21,vjust=5),
        legend.key.size = unit(15,"pt"),
        legend.text=element_text(size=12),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.position = "none") +
  guides(fill = guide_legend(byrow = TRUE))

#### Figure 2C ####
setwd("/home/biostat/DiskA/resource/UKB/BED/total_result")
result = fread("corr_result_stat.txt",header=T,sep="\t")
result = result[result$type!="other",]
result$type = factor(result$type,levels = c("(0,5]","(5,10]","(10,20]",
                                            "(20,50]",">50"),
                     labels = c("(0,5]","(5,10]","(10,20]","(20,50]",">50"))
result

library(ggplot2)

scaleFUN <- function(x) sprintf("%.2f", x)

per_plot = ggplot(result,aes(x=type)) +
  geom_line(aes(y=t_proportion,color="sienna1"),size=0.9,group=1) +
  geom_point(aes(y=t_proportion,color="sienna1"),shape=16,size=8) +
  geom_line(aes(y=g_proportion,color="skyblue3"),size=0.9,group=1) +
  geom_point(aes(y=g_proportion,color="skyblue3"),shape=16,size=8) +
  scale_color_manual(labels = c("TOPMed", "HRC+UK10K"),
                     values=c("salmon","darkturquoise")) +
  scale_y_continuous(limits=c(0,1),labels = scaleFUN) +
  labs(x = "MAC intervals",
       y = "Proportion of  of variants") +
  theme_bw() +
  theme(plot.margin=unit(rep(1,4),'cm'),
        line = element_line(size = 1,linetype = 2),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=13,vjust=1),
        axis.line = element_line(size = 0.2,linetype = "solid",color = "black"),
        axis.text = element_text(color="black",size=19),
        axis.text.x = element_text(angle = 30,vjust = 0.60,hjust = 0.75),
        axis.ticks=element_line(size=0.5, colour="black", linetype="solid"),
        axis.ticks.length = unit(2, "mm"),
        axis.title.x = element_text(size=21,vjust=-3,hjust=0.5),
        axis.title.y = element_text(size=21,vjust=5),
        legend.position = "none",
        legend.key = element_rect(color = NULL,fill=rgb(1,1,1,alpha=0.001)),
        legend.key.size = unit(15, "pt"),
        legend.background = element_rect(fill=rgb(1,1,1,alpha=0.001),color=NULL),
        legend.text = element_text(size=12,color="black"),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))

#### Fig 2D ####
load("/home/biostat/DiskB/jldai/RSQ_stat/imput_stat.RData")
imput_stat$type = factor(imput_stat$type,levels = c("(0,5]","(5,10]","(10,20]","(20,50]",">50 and MAF<0.001"),
                         labels = c("(0,5]","(5,10]","(10,20]","(20,50]",">50"))
# required packages
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(ggtext)
library(grid)
library(ggpubr)

myColors = brewer.pal(4,"YlGnBu")[1:4] 
options(scipen=200)

RSQ_col = ggplot(imput_stat, aes(fill=imput_stat$R2_type, y=N, x=type)) + 
  geom_col(width = 0.7, position = "stack") +
  facet_wrap(~panel,scales = "free_x",ncol = 2)+
  labs(x = "MAC intervals",y = "Numbers of imputed SNVs") +
  scale_fill_manual(name = "",
                    labels = c("0.3-0.4","0.4- 0.6","0.6-0.8","0.8-1"),
                    values = myColors) +
  scale_y_continuous(name = "Numbers of imputed variants",
                     breaks = c(0,20000000,40000000,60000000),
                     labels = c("0M","20M","40M","60M"))+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 19,angle = 30,vjust = 0.70,hjust = 0.75),
        axis.text.y = element_text(size = 19),
        axis.ticks=element_line(size=0.5, colour="gray35", linewidth = 2),
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(size = 0.2,linetype = "solid",color = "black"),
        axis.title.x = element_text(size = 21,vjust=-0.5,hjust=0.5),
        axis.title.y = element_text(size = 21,vjust=1.5,hjust = 0.5),
        strip.background = element_rect(fill = "lightsteelblue3",colour = "white"),
        strip.text = element_text(size = 21,colour = "black"),
        legend.title = element_text(size = 21),
        legend.key.size = unit(15,"pt"),
        legend.text=element_text(size=13),
        legend.position = "none")


############# Plotting code of Figure 3 ###########
#### Fig 3A ####
## Note:should change the files for different ethnicities
library(ggplot2)
library(grid)
library(ggpubr)

# White (shown here as a example) #
result_mean = read.csv("white_result_mean.csv",header = T,sep = ",")

White_vplot = ggplot(result_mean,aes(x=type)) +
  geom_line(aes(y=top_Vmean,color="sienna1"),size=1) +
  geom_point(aes(y=top_Vmean,color="sienna1"),shape=19,size=9.5) +
  geom_line(aes(y=gwas_Vmean,color="skyblue3"),size=1) +
  geom_point(aes(y=gwas_Vmean,color="skyblue3"),shape=19,size=9.5) +
  scale_color_manual(labels = c("TopMed", "HRC/UK10K"),
                     values=c("salmon","darkturquoise")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(breaks = result_mean$type,labels = c("(0,5]","(5,10]","(10,20]",
                                                          "(20,50]",">50")) +
  labs(x = "",
       y = "") +
  theme_bw() +
  theme(plot.margin=unit(rep(1,4),'cm'),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        line = element_line(size = 1,linetype = 2),
        panel.grid.major = element_blank(),
        plot.title = element_text(size=15,vjust=1),
        axis.line = element_line(size = 0.2,linetype = "solid",color = "black"),
        axis.text = element_text(color="black",size=25),
        axis.text.x = element_text(angle = 30,vjust = 0.75,hjust = 0.75),
        axis.ticks=element_line(size=0.3, colour="black", linetype="solid"),
        axis.ticks.length = unit(2, "mm"),
        axis.title.x = element_text(size=21,vjust=-3,hjust=0.5),
        axis.title.y = element_text(size=21,vjust=5),
        legend.position = "none")

#### Fig 3B ####
## Note:should change the files for different ethnicities and different dataset
library(ggplot2)
library(ggtext)
library(ggpubr)
library(ggsci)
library(scico)

# White:hrc (shown here as a example ) #
hrc = fread("white_hrc_rsq.txt")

p2 = ggplot(hrc, aes(x=type, y=g_cramer, fill=type)) + 
  stat_boxplot(geom = "errorbar",width=0.7,aes(x=type,y=g_cramer,group=type)) +
  geom_boxplot(outlier.colour = NA,fill = color) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        plot.title = element_text(size = 12,color = "black"),
        axis.line = element_line(size = 0.2,linetype = "solid",color = "black"),
        axis.text.x = element_text(size = 25,angle = 30,vjust = 0.75,hjust = 0.75),
        axis.text.y = element_text(size = 25,color = "black"),
        axis.ticks.x=element_line(size=0.3, colour="black", linetype="solid"),
        axis.ticks.length.x = unit(2, "mm"),
        axis.ticks.y=element_line(size=0.3, colour="black", linetype="solid"),
        axis.ticks.length = unit(2, "mm"),
        axis.title.y = element_text(size = 8,vjust=2.5,hjust = 0.5),
        legend.position = "none")

#### Fig 3C ####
## Note:should change the files for different ethnicities and differenet dataset
color = scico(7, palette = 'berlin',direction = 1)
color = brewer.pal(7,"YlGnBu")[1:7] 

scaleFUN <- function(x) sprintf("%.2f", x)

# White:hrc (shown here as a example) #
hrc = fread("white_hrc_rsq.txt")
hrc = as.data.frame(hrc %>%
                      group_by(interval,type) %>%
                      summarise(cramer = mean(g_cramer),.groups = "keep"))
hrc$type = factor(hrc$type,levels=c("(0.3,0.4]","(0.4,0.5]","(0.5,0.6]",
                                    "(0.6,0.7]","(0.7,0.8]","(0.8,0.9]","(0.9,1]"))

w_hrc = ggplot(hrc, aes(fill = type, y = cramer, x = interval)) + 
  geom_col(position = "dodge",width = 0.85) +   
  labs( title = "",
        x = "",
        y = "") +
  scale_fill_manual(name = "INFO/RSQ",labels = c("(0.3,0.4]","(0.4,0.5]","(0.5,0.6]",
                                                 "(0.6,0.7]","(0.7,0.8]","(0.8,0.9]","(0.9,1]"),
                    values = color,guide = guide_legend(title.position = "top", nrow = 1)) +
  scale_y_continuous(limits=c(0,1),labels = scaleFUN) +
  scale_x_continuous(breaks = c(1,2,3,4,5),labels = c("(0,5]","(5,10]","(10,20]","(20,50]",">50")) +  
  theme_bw() +
  theme(plot.margin=unit(rep(1,4),'cm'),
        line = element_line(size = 1,linetype = 2),
        panel.grid = element_blank(),   		
        panel.border = element_blank(),
        plot.title = element_text(size = 15,color = "black",vjust = 2),
        axis.line = element_line(size = 0.2,linetype = "solid",color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.box = "horizontal",
        legend.key.size = unit(9,"pt"),
        legend.text=element_text(size=9),
        legend.title = element_text(size = 10),
        legend.position = "top")
