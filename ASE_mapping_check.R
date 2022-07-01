library(plyr)
library(tidyverse)
library(ggrepel)


remove(list = ls())
load("starting.counts.ASDE.RData")


counts.libwise <- counts.allele%>%select(-GeneSymbol)%>%group_by(allele)%>%summarise_all(.funs = sum)
counts.libwise <- as.data.frame(t(counts.libwise))
colnames(counts.libwise)<-counts.libwise[1,]
counts.libwise<-counts.libwise[-1,]
counts.libwise[,1]<-as.numeric(counts.libwise[,1])
counts.libwise[,2]<-as.numeric(counts.libwise[,2])
counts.libwise$lib <- substr(x = rownames(counts.libwise),start = 4,stop = 9) 
counts.libwise$species <- substr(x = rownames(counts.libwise),start = 7,stop = 7)
counts.libwise$species <- gsub(pattern = "H",replacement = "Hybrid",x = counts.libwise$species)
counts.libwise$species <- gsub(pattern = "S",replacement = "*D. simulans*",x = counts.libwise$species)
counts.libwise$species <- gsub(pattern = "M",replacement = "*D. melanogaster*",x = counts.libwise$species)


ggplot(data = counts.libwise,
       mapping = aes(x = log2(`*D. melanogaster*`),y = log2(`*D. simulans*`)))+
  geom_point(aes(colour=species))+
  geom_abline(slope = 1,intercept = 0)
