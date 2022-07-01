library(edgeR)
library(plyr)
library(tidyverse)
library(ggrepel)




remove(list = ls())
load("starting_data.RData")


counts.allele<-merge(x = uni_Dmel[,-1],y = uni_Dsim[,-1],by = 0,suffixes = c(".m",".s"))
rownames(counts.allele)<-counts.allele[,1]
counts.allele<-counts.allele[,-1]
counts.allele<-counts.allele[,!grepl("_M...s",colnames(counts.allele))]
counts.allele<-counts.allele[,!grepl("_S...m",colnames(counts.allele))]
save(counts.allele,file = "starting.counts.ASDE.RData")

uni_Dmel[,grepl("_S",names(uni_Dmel))] <- 0
uni_Dsim[,grepl("_M",names(uni_Dsim))] <- 0
counts.geno <- bind_rows(uni_Dmel%>%`row.names<-`(.,NULL),
                         uni_Dsim%>%`row.names<-`(.,NULL))%>%group_by(GeneSymbol)%>%summarise_all(.funs = sum)
save(counts.geno,file = "starting.counts.DE.RData")

