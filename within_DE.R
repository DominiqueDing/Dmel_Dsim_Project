library(plyr)
library(tidyverse)
library(ggrepel)
library(ggtext)
options(ggrepel.max.overlaps = Inf)
library(ggpubr)
library(scales)


remove(list = ls())
load("starting.counts.DE.RData") # wrongly mapped reads removed in parental libs

counts <- as.data.frame(counts.geno)
rownames(counts) <- counts[,1]
counts<-counts[,-1]
counts<-counts%>%mutate_all(.funs = as.numeric)




###################################################################
####################################### creating DGEList ##########
###################################################################
y <- DGEList(counts = counts,group = substr(x=names(counts),start = 7,stop = 9))
# keep <- row.names(y$counts)%in%c(genes$gene)
# table(keep)
# keep
# FALSE  TRUE 
# 2355  4969 
# y<-y[keep, , keep.lib.sizes=FALSE]
# y<-calcNormFactors(y)

#need min 3 samples with cpm >5
keep<-rowSums(cpm(y)>4) >= 3
table(keep)
# keep
# FALSE  TRUE 
# 3646  9169
y<-y[keep, , keep.lib.sizes=FALSE]
y<-calcNormFactors(y)
y$samples
included_Genes<-rownames(y$counts)

### MD plot: column=sample(i)
plotMD(cpm(y, log=TRUE), column=2)
abline(h=0, col="red", lty=2, lwd=2)

###################################################################
################# design ##########################################
###################################################################
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
### Estimating the dispersion:
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
# [1] 0.02269442
plotBCV(y)
### fit quasi-likelihood model
fit<-glmQLFit(y,design = design,robust=TRUE)
plotQLDisp(fit)


#############################
#####  ################
#############################
### set contrast:
con.DE <- makeContrasts(
  con.DE.sim = (Sec+Sml+Spr)/3-Sus,
  con.DE.mel = (Mec+Mml+Mpr)/3-Mus,
  con.DE.hyb = (Hec+Hml+Hpr)/3-Hus,
  
  # con.DE.hyb_bc.pr = (Hec+Hml)/2-Hpr,
  # con.DE.hyb_ec.pr = Hec-Hpr,
  # con.DE.hyb_ml.pr = Hml-Hpr,
  
  con.DE.hyb.mel.0 = Hus-Mus,
  con.DE.hyb.sim.0 = Hus-Sus,
  con.DE.hyb.mel.1 = (Hec+Hml+Hpr)/3-(Mec+Mml+Mpr)/3,
  con.DE.hyb.sim.1 = (Hec+Hml+Hpr)/3-(Sec+Sml+Spr)/3,
  con.DE.hyb.mel.r = ((Hec+Hml+Hpr)/3-Hus)-((Mec+Mml+Mpr)/3-Mus),
  con.DE.hyb.sim.r = ((Hec+Hml+Hpr)/3-Hus)-((Sec+Sml+Spr)/3-Sus),
  
  con.DE.mel.sim.r = ((Mec+Mml+Mpr)/3-Mus)-((Sec+Sml+Spr)/3-Sus),
  
  levels = design
)
### QL F-test:
for (con in colnames(con.DE)) {
  .GlobalEnv[[paste0('QLFT_DE.',substr(x = con,start = 8,stop = nchar(con)))]] <- glmQLFTest(fit,contrast = con.DE[,con])
  .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]] <- topTags(object = .GlobalEnv[[paste0('QLFT_DE.',substr(x = con,start = 8,stop = nchar(con)))]],n = Inf)$table
  .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]<-.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]%>%
    mutate(dmel_gene = rownames(.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]))
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM23",replacement = "BomBc1",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM14",replacement = "Dso2",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM1",replacement = "BomS1",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM2",replacement = "BomS2",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "CG15065",replacement = "BomS5",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "CG16836",replacement = "BomT2",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM4",replacement = "Dso1",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
 .GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IMPPP",replacement = "BaraA2",x =.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]]$dmel_gene)
  rownames(.GlobalEnv[[paste0('df.result_DE.',substr(x = con,start = 8,stop = nchar(con)))]])<-c()
}


### arrnage results:

result.within.DE <- bind_rows(df.result_DE.sim%>%mutate(species='*D. simulans*'),
                              df.result_DE.mel%>%mutate(species='*D. melanogaster*'),
                              df.result_DE.hyb%>%mutate(species='Hybrid'))
result.within.DE$dmel_gene<-gsub(pattern = "IM23",replacement = "BomBc1",x = result.within.DE$dmel_gene)
result.within.DE$dmel_gene<-gsub(pattern = "IM14",replacement = "Dso2",x = result.within.DE$dmel_gene)
result.within.DE$dmel_gene<-gsub(pattern = "IM1",replacement = "BomS1",x = result.within.DE$dmel_gene)
result.within.DE$dmel_gene<-gsub(pattern = "IM2",replacement = "BomS2",x = result.within.DE$dmel_gene)
result.within.DE$dmel_gene<-gsub(pattern = "CG15065",replacement = "BomS5",x = result.within.DE$dmel_gene)
result.within.DE$dmel_gene<-gsub(pattern = "CG16836",replacement = "BomT2",x = result.within.DE$dmel_gene)
result.within.DE$dmel_gene<-gsub(pattern = "IM4",replacement = "Dso1",x = result.within.DE$dmel_gene)
result.within.DE$dmel_gene<-gsub(pattern = "IMPPP",replacement = "BaraA2",x = result.within.DE$dmel_gene)


result.within.DE$DE_trt.vs.ctrl[result.within.DE$FDR<0.05&result.within.DE$logFC>0] <- "Up regulated"
result.within.DE$DE_trt.vs.ctrl[result.within.DE$FDR<0.05&result.within.DE$logFC >= 1] <- "Up regulated 1"
# result.within.DE$DE_trt.vs.ctrl[result.within.DE$FDR<0.05&result.within.DE$logFC >= 1.5] <- "Up regulated 1.5"
result.within.DE$DE_trt.vs.ctrl[result.within.DE$FDR<0.05&result.within.DE$logFC<0] <- "Down regulated"
result.within.DE$DE_trt.vs.ctrl[result.within.DE$FDR<0.05&result.within.DE$logFC<= -1] <- "Down regulated 1"
# result.within.DE$DE_trt.vs.ctrl[result.within.DE$FDR<0.05&result.within.DE$logFC<= -1.5] <- "Down regulated 1.5"
result.within.DE<-result.within.DE%>%
  mutate(species=factor(species,levels = c("*D. melanogaster*","Hybrid","*D. simulans*")))
# write.table(x = unique(c(result.within.DE$dmel_gene[result.within.DE$DE_trt.vs.ctrl=="Up regulated 1"],
#                          result.within.DE$dmel_gene[result.within.DE$DE_trt.vs.ctrl=="Down regulated 1"])),
#             file = "genes.sig.in.any.logfc1.txt",quote = F,row.names = F,col.names = F)
# write.table(x = unique(result.within.DE$dmel_gene[!is.na(result.within.DE$DE_trt.vs.ctrl)]),
# file = "immune.responsive.genes.txt",quote = F,row.names = F,col.names = F)
# write.table(x = unique((result.within.DE%>%filter(result.within.DE$FDR<0.05,abs(result.within.DE$logFC) >= 1.5))$dmel_gene),
#             file = "genes.sig.in.any.logfc1.5.txt",quote = F,row.names = F,col.names = F)








plot.colours<-c(
  "Cis only"="#000000",
  "Cis*"="#333333",
  "Cis**"="#707070",
  "Trans only"="#C7050F",
  "Trans*"="#FF5470",
  "Trans**"="#F88379",
  "Cis+Trans"="#8C12E2",
  "CisxTrans"="#2CA03D",
  "Compensatory"="#FFAB00",
  "Conserved"="gray90",
  "Ambiguous"="#0DD9E7",
  
  
  
  "Down regulated"="dodgerblue4",
  "Down regulated 1"="dodgerblue",
  "Down regulated 1.5"="dodgerblue",
  
  "Up regulated"="firebrick",
  "Up regulated 1"="firebrick1",
  "Up regulated 1.5"="firebrick1",
  
  "No significant difference "="black",
  "No significant difference"="black"
)






####################################
######### Within Species volcano ###
####################################
pdf(file = "Within.species.DE.trt.vs.ctrl.volcano.plot.pdf",width = 10,height = 5)
ggplot(data = result.within.DE,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     colour=DE_trt.vs.ctrl))+
  geom_point()+
  # geom_text_repel(aes(label=ifelse(test = abs(logFC)>2,yes = dmel_gene,no = '')))+
  geom_vline(xintercept = 1.5,linetype=3,colour="gray10")+
  geom_vline(xintercept = -1.5,linetype=3,colour="gray10")+
  facet_grid(~species)+
  scale_colour_manual(values = plot.colours)+
  labs(x="Log<sub>2</sub> Fold Change of Gene Expression Against Immune Challenge",
       y="-Log<sub>10</sub>(P value)")+
  # xlim(-15,15)+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        legend.position = 'none',
        strip.background = element_rect(fill='NA',colour = 'black'),
        panel.grid = element_blank(),
        strip.text = element_markdown(size = 20),
        axis.title.x = element_markdown(size = 20),
        axis.title.y = element_markdown(size = 20),
        axis.text = element_text(size = 20))
dev.off()





ggplot(data = result.within.DE%>%filter(species=="Hybrid",!is.na(DE_trt.vs.ctrl)))+
  geom_histogram(aes(x = abs(logFC)),bins = 100)


pdf(file = "Within.species.DE.trt.vs.ctrl.volcano.plot_scratch.pdf",width = 10,height = 10)
ggplot(data = result.within.DE,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     colour=DE_trt.vs.ctrl))+
  geom_point()+
  geom_text_repel(aes(label=ifelse(test = abs(logFC)>2,yes = dmel_gene,no = '')))+
  facet_grid(~species)+
  scale_colour_manual(values = plot.colours)+
  labs(x="Log<sub>2</sub> Fold Change of Gene Expression Against Immune Challenge",
       y="-Log<sub>10</sub>(P value)")+
  # xlim(-15,15)+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        legend.position = 'none',
        strip.background = element_rect(fill='NA',colour = 'black'),
        panel.grid = element_blank(),
        strip.text = element_markdown(size = 20),
        axis.title.x = element_markdown(size = 20),
        axis.title.y = element_markdown(size = 20),
        axis.text = element_text(size = 20))
dev.off()


######################################################
#### Correlation between specific immune response: ###
######################################################
genes.sig.in.f1 <- result.within.DE$dmel_gene[result.within.DE$species=="Hybrid"&!is.na(result.within.DE$DE_trt.vs.ctrl)]
genes.sig.in.mel <- result.within.DE$dmel_gene[result.within.DE$species=="*D. melanogaster*"&!is.na(result.within.DE$DE_trt.vs.ctrl)]
genes.sig.in.sim <- result.within.DE$dmel_gene[result.within.DE$species=="*D. simulans*"&!is.na(result.within.DE$DE_trt.vs.ctrl)]
genes.sig.in.any <- unique(c(genes.sig.in.f1,genes.sig.in.mel,genes.sig.in.sim))
genes.sig.in.any.logfc1 <- unique(c(result.within.DE$dmel_gene[result.within.DE$DE_trt.vs.ctrl=="Up regulated 1"],result.within.DE$dmel_gene[result.within.DE$DE_trt.vs.ctrl=="Down regulated 1"]))
genes.sig.in.any.logfc1 <- genes.sig.in.any.logfc1[!is.na(genes.sig.in.any.logfc1)]

genes.up.in.f1 <- unique(result.within.DE$dmel_gene[result.within.DE$species=="Hybrid"&grepl("Up regulated",result.within.DE$DE_trt.vs.ctrl)])
genes.up.in.mel <- unique(result.within.DE$dmel_gene[result.within.DE$species=="*D. melanogaster*"&grepl("Up regulated",result.within.DE$DE_trt.vs.ctrl)])
genes.up.in.sim <- unique(result.within.DE$dmel_gene[result.within.DE$species=="*D. simulans*"&grepl("Up regulated",result.within.DE$DE_trt.vs.ctrl)])
genes.dw.in.f1 <- unique(result.within.DE$dmel_gene[result.within.DE$species=="Hybrid"&grepl("Down regulated",result.within.DE$DE_trt.vs.ctrl)])
genes.dw.in.mel <- unique(result.within.DE$dmel_gene[result.within.DE$species=="*D. melanogaster*"&grepl("Down regulated",result.within.DE$DE_trt.vs.ctrl)])
genes.dw.in.sim <- unique(result.within.DE$dmel_gene[result.within.DE$species=="*D. simulans*"&grepl("Down regulated",result.within.DE$DE_trt.vs.ctrl)])



hist(result.within.DE$logFC[result.within.DE$dmel_gene%in%c(genes.sig.in.f1,genes.sig.in.mel)&result.within.DE$species=="Hybrid"])
hist(result.within.DE$logFC[result.within.DE$dmel_gene%in%c(genes.sig.in.f1,genes.sig.in.mel)&result.within.DE$species=="*D. melanogaster*"])






plot.dmel.dsim.data <- merge(result.within.DE%>%filter(species=="*D. melanogaster*")%>%select(dmel_gene,logFC,DE_trt.vs.ctrl),
                            result.within.DE%>%filter(species=="*D. simulans*")%>%select(dmel_gene,logFC,DE_trt.vs.ctrl),
                            by ="dmel_gene",all = T,suffixes = c("(*D. melanogaster*)","(*D. simulans*)"))%>%
  filter(dmel_gene%in%genes.sig.in.any)%>%mutate(col = "")
plot.dmel.dsim.data$col[plot.dmel.dsim.data$dmel_gene%in%genes.sig.in.any]<-"any3"
# plot.dmel.dsim.data$col[plot.dmel.dsim.data$dmel_gene%in%c(genes.sig.in.mel,genes.sig.in.sim)]<-"anyPa"
plot.dmel.dsim.data$col[plot.dmel.dsim.data$dmel_gene%in%c(genes.sig.in.any.logfc1)]<-"any3.logfc1"

cor.test(plot.dmel.dsim.data$`logFC(*D. melanogaster*)`[plot.dmel.dsim.data$col=="any3"],
         plot.dmel.dsim.data$`logFC(*D. simulans*)`[plot.dmel.dsim.data$col=="any3"])
cor.test(plot.dmel.dsim.data$`logFC(*D. melanogaster*)`[plot.dmel.dsim.data$col=="any3.logfc1"],
         plot.dmel.dsim.data$`logFC(*D. simulans*)`[plot.dmel.dsim.data$col=="any3.logfc1"])

plot.dmel.dsim<-
  ggplot()+
  # geom_point(data=plot.dmel.dsim.data%>%filter(col=="any3"),aes(x = `logFC(*D. melanogaster*)`,y = `logFC(*D. simulans*)`))+
  geom_point(data=plot.dmel.dsim.data%>%filter(col=="any3"),aes(x = `logFC(*D. melanogaster*)`,y = `logFC(*D. simulans*)`),colour="gray")+
  geom_point(data=plot.dmel.dsim.data%>%filter(col=="any3.logfc1"),aes(x = `logFC(*D. melanogaster*)`,y = `logFC(*D. simulans*)`),colour="black")+
  # scale_colour_manual(values = list("any3"="gray","anyPa"="gray","any3.logfc1"="black"))+
  geom_abline(slope = 1,intercept = 0,linetype=3,colour="dark grey")+
  geom_vline(xintercept = 0,linetype=3,colour="dark grey")+
  geom_hline(yintercept = 0,linetype=3,colour="dark grey")+
  annotate(geom = "text",x=-4,y=4.8,label="italic(r)==0.48",parse=T, colour="gray",size=7)+
  annotate(geom = "text",x=-4,y=4.0,label="italic(r)==0.54",parse=T,colour="black",size=7)+
  # geom_text_repel(aes(label=if_else(condition = (abs(`logFC(*D. melanogaster*)`)>1.2)|(abs(`logFC(*D. simulans*)`)>1.2),
  # true = dmel_gene,false = '')))+
  xlim(-5.1,5.1)+ylim(-5.1,5.1)+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        legend.position = 'none',
        strip.background = element_rect(fill='NA',colour = 'black'),
        panel.grid = element_blank(),
        strip.text = element_markdown(size = 20),
        axis.title.x = element_markdown(size = 20),
        axis.title.y = element_markdown(size = 20),
        axis.text = element_text(size = 20))








plot.dmel.hyb.data <- merge(result.within.DE%>%filter(species=="*D. melanogaster*")%>%select(dmel_gene,logFC,DE_trt.vs.ctrl),
                            result.within.DE%>%filter(species=="Hybrid")%>%select(dmel_gene,logFC,DE_trt.vs.ctrl),
                            by ="dmel_gene",all = T,suffixes = c("(*D. melanogaster*)","(Hybrid)"))%>%
  filter(dmel_gene%in%genes.sig.in.any)%>%mutate(col = "")
plot.dmel.hyb.data$col[plot.dmel.hyb.data$dmel_gene%in%genes.sig.in.any]<-"any3"
# plot.dmel.hyb.data$col[plot.dmel.hyb.data$dmel_gene%in%c(genes.sig.in.mel,genes.sig.in.sim)]<-"anyPa"
plot.dmel.hyb.data$col[plot.dmel.hyb.data$dmel_gene%in%c(genes.sig.in.any.logfc1)]<-"any3.logfc1"

cor.test(plot.dmel.hyb.data$`logFC(*D. melanogaster*)`[plot.dmel.hyb.data$col=="any3"],
         plot.dmel.hyb.data$`logFC(Hybrid)`[plot.dmel.hyb.data$col=="any3"])
cor.test(plot.dmel.hyb.data$`logFC(*D. melanogaster*)`[plot.dmel.hyb.data$col=="any3.logfc1"],
         plot.dmel.hyb.data$`logFC(Hybrid)`[plot.dmel.hyb.data$col=="any3.logfc1"])

plot.dmel.hyb<-
  ggplot()+
  # geom_point(data=plot.dmel.hyb.data%>%filter(col=="any3"),aes(x = `logFC(*D. melanogaster*)`,y = `logFC(Hybrid)`))+
  geom_point(data=plot.dmel.hyb.data%>%filter(col=="any3"),aes(x = `logFC(*D. melanogaster*)`,y = `logFC(Hybrid)`),colour="gray")+
  geom_point(data=plot.dmel.hyb.data%>%filter(col=="any3.logfc1"),aes(x = `logFC(*D. melanogaster*)`,y = `logFC(Hybrid)`),colour="black")+
  # scale_colour_manual(values = list("any3"="gray","anyPa"="black"))+
  geom_abline(slope = 1,intercept = 0,linetype=3,colour="dark grey")+
  geom_vline(xintercept = 0,linetype=3,colour="dark grey")+
  geom_hline(yintercept = 0,linetype=3,colour="dark grey")+
  annotate(geom = "text",x=-4,y=4.8,label="italic(r)==0.19",parse=T, colour="gray",size=7)+
  annotate(geom = "text",x=-4,y=4.0,label="italic(r)==0.53",parse=T,colour="black",size=7)+
  # geom_text_repel(aes(label=if_else(condition = (abs(`logFC(*D. melanogaster*)`)>1.2)|(abs(`logFC(Hybrid)`)>1.2),
  # true = dmel_gene,false = '')))+
  xlim(-5.1,5.1)+ylim(-5.1,5.1)+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        legend.position = 'none',
        strip.background = element_rect(fill='NA',colour = 'black'),
        panel.grid = element_blank(),
        strip.text = element_markdown(size = 20),
        axis.title.x = element_markdown(size = 20),
        axis.title.y = element_markdown(size = 20),
        axis.text = element_text(size = 20))




plot.dsim.hyb.data <- merge(result.within.DE%>%filter(species=="*D. simulans*")%>%select(dmel_gene,logFC,DE_trt.vs.ctrl),
                            result.within.DE%>%filter(species=="Hybrid")%>%select(dmel_gene,logFC,DE_trt.vs.ctrl),
                            by ="dmel_gene",all = T,suffixes = c("(*D. simulans*)","(Hybrid)"))%>%
  filter(dmel_gene%in%genes.sig.in.any)%>%mutate(col = "")
plot.dsim.hyb.data$col[plot.dsim.hyb.data$dmel_gene%in%genes.sig.in.any]<-"any3"
plot.dsim.hyb.data$col[plot.dsim.hyb.data$dmel_gene%in%genes.sig.in.any.logfc1]<-"any3.logfc1"
# plot.dsim.hyb.data$col[plot.dsim.hyb.data$dmel_gene%in%c(genes.sig.in.mel,genes.sig.in.sim)]<-"anyPa"

cor.test(plot.dsim.hyb.data$`logFC(*D. simulans*)`[plot.dsim.hyb.data$col=="any3"],
         plot.dsim.hyb.data$`logFC(Hybrid)`[plot.dsim.hyb.data$col=="any3"])
cor.test(plot.dsim.hyb.data$`logFC(*D. simulans*)`[plot.dsim.hyb.data$col=="any3.logfc1"],
         plot.dsim.hyb.data$`logFC(Hybrid)`[plot.dsim.hyb.data$col=="any3.logfc1"])

plot.dsim.hyb<-
  ggplot()+
  # geom_point(data=plot.dsim.hyb.data%>%filter(col=="any3"),aes(x = `logFC(*D. simulans*)`,y = `logFC(Hybrid)`))+
  geom_point(data=plot.dsim.hyb.data%>%filter(col=="any3"),aes(x = `logFC(*D. simulans*)`,y = `logFC(Hybrid)`),colour="gray")+
  geom_point(data=plot.dsim.hyb.data%>%filter(col=="any3.logfc1"),aes(x = `logFC(*D. simulans*)`,y = `logFC(Hybrid)`),colour="black")+
  # scale_colour_manual(values = list("any3"="gray","anyPa"="black"))+
  geom_abline(slope = 1,intercept = 0,linetype=3,colour="dark grey")+
  geom_vline(xintercept = 0,linetype=3,colour="dark grey")+
  geom_hline(yintercept = 0,linetype=3,colour="dark grey")+
  annotate(geom = "text",x=-4,y=4.8,label="italic(r)==0.39",parse=T, colour="gray",size=7)+
  annotate(geom = "text",x=-4,y=4.0,label="italic(r)==0.70",parse=T,colour="black",size=7)+
  # geom_text_repel(aes(label=if_else(condition = (abs(`logFC(*D. simulans*)`)>1.2)|(abs(`logFC(Hybrid)`)>1.2),
                                    # true = dmel_gene,false = '')))+
  xlim(-5.1,5.1)+ylim(-5.1,5.1)+
  theme_bw()+
  theme(panel.spacing = unit(0,'lines'),
        legend.position = 'none',
        strip.background = element_rect(fill='NA',colour = 'black'),
        panel.grid = element_blank(),
        strip.text = element_markdown(size = 20),
        axis.title.x = element_markdown(size = 20),
        axis.title.y = element_markdown(size = 20),
        axis.text = element_text(size = 20))










pdf(file = "within.speciesDE.correlation.pdf",width = 16,height = 5)
ggarrange(plot.dmel.hyb,plot.dsim.hyb,plot.dmel.dsim,ncol = 3,nrow = 1)
dev.off()










df.result_DE.mel.sim.r.logfc1 <- df.result_DE.mel.sim.r%>%filter(dmel_gene%in%genes.sig.in.any.logfc1)
ggplot(data = df.result_DE.mel.sim.r.logfc1)+
  geom_histogram(aes(x=PValue),bins=40)
ggplot(data = df.result_DE.mel.sim.r.logfc1)+
  geom_point(aes(x = logFC,y = -log10(PValue)))



