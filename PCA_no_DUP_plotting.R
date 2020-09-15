setwd("I:/Quinoa/Results/SNP_Calling_V2/Final_data_analysis/PCA")
library(tidyverse)
library(ggplot2)

PCA<-read.csv("pca_all_no_dup_GWAS.csv")

country<-read.csv("country.csv")

head(PCA)

highland<-read.table("highland.txt",sep="\t",header = F)

head(highland)

highland %>% mutate(pop = paste("highland"))->highland

lowland<-read.table("lowland.txt",sep="\t",header = F)

head(lowland)

lowland %>% mutate(pop = paste("lowland"))-> lowland

rbind(highland,lowland)-> pop_info

colnames(pop_info)<-c("id","pop")

left_join(PCA,pop_info,by="id") %>% left_join(.,country,by="id")->pca

 write.csv(pca,"PCA_no_duplicates_2.csv")

#pca<-read.csv("PCA_no_duplicates.csv")
 
 pca2<-read.csv("PCA_no_duplicates_2.csv")

head(pca)

ggplot(pca,aes(PC1,PC2,color=country))+geom_point()+
  scale_x_continuous(name="PC1 (23.35%)")+
  scale_y_continuous(name="PC2 (9.45%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme( axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))


ggplot(pca2,aes(PC1,PC2,color=Country))+geom_point(size = 2)+
  scale_x_continuous(name="PC1 (23.35%)")+
  scale_y_continuous(name="PC2 (9.45%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=16, color = "black"))+
  theme(axis.text.x = element_text(colour = "black", size = 12))+
  theme(axis.text.y = element_text(colour = "black", size = 12))+
  theme( axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))+
  scale_color_manual(values = c("#00E7E2","#FA0E0E","#08EE1E",
                                "#0F46F5", "#660066","#1C6223","#F961F9",
                                "#EEFC36","#C6AF28","#FF9900","#666666"))

#Argentina - #00E7E2
#Australia - #FA0E0E
#Bolivia - #08EE1E
#Chile - #0F46F5
#Denmark - #660066
#Ecudor - #1C6223
#Peru - #F961F9
#Switzerland - #EEFC36 (Tree	#FFA500)
#United Kingdom - #C6AF28 
#USA - #FF9900 
#Unknown - #666666




