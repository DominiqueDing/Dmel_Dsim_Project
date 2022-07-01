library(edgeR)
library(plyr)
library(tidyverse)
library(ggrepel)
library(ggtext)
options(ggrepel.max.overlaps = Inf)
library(ggpubr)
library(cowplot)
library(scales)


remove(list = ls())
load("starting.counts.ASDE.RData") # wrongly mapped reads removed in parental libs
starting.genes<-read.table("starting.genes.txt")
counts<-as.matrix(counts.allele)

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
  
  "*Cis* divergent"="#000000",
  "*Trans* divergent"="#C7050F",
  "Both *Cis* and *Trans* divergent"="gold",
  "*Cis*"="#000000",
  "*Trans*"="#C7050F",
  
  
  
  "Down regulated "="dodgerblue",
  "Up regulated "="firebrick"

)





###################################################################
####################################### creating DGEList ##########
###################################################################
y <- DGEList(counts = counts,group = substr(x=colnames(counts),start = 7,stop = 11))
keep <- row.names(y$counts)%in%c(starting.genes$V1)
table(keep)
# keep
# FALSE  TRUE 
# 3646  9169  
y<-y[keep, , keep.lib.sizes=FALSE]
y<-calcNormFactors(y)

# #need min 3 samples with cpm >5
# keep<-rowSums(cpm(y)>5) >= 3
# table(keep)
# # keep
# # FALSE  TRUE 
# # 3870  8945
# y<-y[keep, , keep.lib.sizes=FALSE]
# y<-calcNormFactors(y)



### MD plot: column=sample(i)
plotMD(cpm(y, log=TRUE), column=2)
abline(h=0, col="red", lty=2, lwd=2)






 



y.parents <- DGEList(counts = counts[,!grepl(pattern = "_H",x = colnames(counts))],
                     group = substr(x=colnames(counts[,!grepl(pattern = "_H",x = colnames(counts))]),start = 7,stop = 11))
keep <- row.names(y.parents$counts)%in%c(starting.genes$V1)
y.parents<-y.parents[keep, , keep.lib.sizes=FALSE]
y.parents<-calcNormFactors(y.parents)
col.mds <- gsub(pattern = "H",replacement = "#089911",x = substr(y.parents$samples$group,1,1))
col.mds <- gsub(pattern = "M",replacement = "#4D9DE0",x = col.mds)
col.mds <- gsub(pattern = "S",replacement = "#f8766d",x = col.mds)
lab.mds <- gsub(pattern = "pr",replacement = "Stabbed",x = substr(y.parents$samples$group,2,3))
lab.mds <- gsub(pattern = "us",replacement = "Control",x = lab.mds)
lab.mds <- gsub(pattern = "ml",replacement = "M. luteus",x = lab.mds)
lab.mds <- gsub(pattern = "ec",replacement = "E. coli",x = lab.mds)
pdf("mds_parents.pdf",width = 4,height = 4)
par(mar=c(4,4,1,1))
plotMDS(y.parents,
            cex = 2,lwd=2,
            pch = as.numeric(as.factor(col.mds)),
            col = as.numeric(as.factor(lab.mds)),
            xlim=c(-3,2.5),
            las=1)
dev.off()


y.hyb <- DGEList(counts = counts[,grepl(pattern = "_H",x = colnames(counts))],
                     group = substr(x=colnames(counts[,grepl(pattern = "_H",x = colnames(counts))]),start = 7,stop = 11))
keep <- row.names(y.hyb$counts)%in%c(starting.genes$V1)
y.hyb<-y.hyb[keep, , keep.lib.sizes=FALSE]
y.hyb<-calcNormFactors(y.hyb)
col.mds <- gsub(pattern = "m",replacement = "#4D9DE0",x = substr(y.hyb$samples$group,5,5))
col.mds <- gsub(pattern = "s",replacement = "#f8766d",x = col.mds)
lab.mds <- gsub(pattern = "pr",replacement = "Stabbed",x = substr(y.hyb$samples$group,2,3))
lab.mds <- gsub(pattern = "us",replacement = "Control",x = lab.mds)
lab.mds <- gsub(pattern = "ml",replacement = "M. luteus",x = lab.mds)
lab.mds <- gsub(pattern = "ec",replacement = "E. coli",x = lab.mds)

pdf("mds_f1.pdf",width = 5.5,height = 4)
par(mar=c(4,4,1,8.7),xpd=TRUE)
plotMDS(y.hyb,
            cex = 2,lwd=2,
            pch = as.numeric(as.factor(col.mds)),
            col = as.numeric(as.factor(lab.mds)),
            xlim=c(-3,2.5),
            las=1)
legend("bottomright", inset=c(-0.6,0),box.lwd = 0,box.lty = 0,pt.cex = 1.5,pt.lwd = 2,
       legend=c((lab.mds)[c(2,3,1,4)],expression(italic(D.~melanogaster)),expression(italic(D.~simulans))),
       pch=c(rep(15,4),1,2),
       col=c(as.numeric(factor(lab.mds))[c(2,3,1,4)],"grey","grey"),
       ncol = 1,cex = 1)
dev.off()



###################################################################
################# design ##########################################
###################################################################
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
### Estimating the dispersion:
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
# [1] 0.0252498
plotBCV(y)
### fit quasi-likelihood model
fit<-glmQLFit(y,design = design,robust=TRUE)
plotQLDisp(fit)


#############################
#####  ################
#############################
### set contrast:
con.ASDE <- makeContrasts(
  con.ctrl.cis = Hus.m - Hus.s,
  con.ctrl.total = Mus.m - Sus.s,
  con.ctrl.trans = (Mus.m - Sus.s)-(Hus.m - Hus.s),
  
  con.trt.cis = (Hec.m+Hus.m+Hpr.m)/3 - (Hec.s+Hus.s+Hpr.s)/3,
  con.trt.total = (Mec.m+Mus.m+Mpr.m)/3 - (Sec.s+Sus.s+Spr.s)/3,
  con.trt.trans = ((Mec.m+Mus.m+Mpr.m)/3 - (Sec.s+Sus.s+Spr.s)/3)-((Hec.m+Hus.m+Hpr.m)/3 - (Hec.s+Hus.s+Hpr.s)/3),
  
  con.reac.cis = ((Hec.m+Hus.m+Hpr.m)/3 - (Hec.s+Hus.s+Hpr.s)/3)-(Hus.m - Hus.s),
  con.reac.total = ((Mec.m+Mus.m+Mpr.m)/3 - (Sec.s+Sus.s+Spr.s)/3)-(Mus.m - Sus.s),
  con.reac.trans = (((Mec.m+Mus.m+Mpr.m)/3 - (Sec.s+Sus.s+Spr.s)/3)-(Mus.m - Sus.s)) - (((Hec.m+Hus.m+Hpr.m)/3 - (Hec.s+Hus.s+Hpr.s)/3)-(Hus.m - Hus.s)),
  
  
  levels = design
)
### QL F-test:
for (con in colnames(con.ASDE)) {
  .GlobalEnv[[paste0('QLFT_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]] <- glmQLFTest(fit,contrast = con.ASDE[,con])
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]] <- topTags(object = .GlobalEnv[[paste0('QLFT_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]],n = Inf)$table
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]<-.GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]%>%
    mutate(dmel_gene = rownames(.GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]))
  rownames(.GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]])<-c()
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM23",replacement = "BomBc1",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM14",replacement = "Dso2",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM1",replacement = "BomS1",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM2",replacement = "BomS2",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "CG15065",replacement = "BomS5",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "CG16836",replacement = "BomT2",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IM4",replacement = "Dso1",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
  .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene<-gsub(pattern = "IMPPP",replacement = "BaraA2",x = .GlobalEnv[[paste0('df.result_ASDE.',substr(x = con,start = 5,stop = nchar(con)))]]$dmel_gene)
}


# ### arrange results:
# result.ASDE.ctrl <- merge(
#   x = merge(x = df.result_ASDE.ctrl.cis,y = df.result_ASDE.ctrl.total,by = "dmel_gene",all = T,suffixes = c(".cis",".total"),),
#   y = df.result_ASDE.ctrl.trans%>%`colnames<-`(paste0(colnames(.),".trans")),
#   by.x = "dmel_gene",by.y = "dmel_gene.trans",all = T
# )%>%mutate(treatment="Control")
# 
# result.ASDE.trt <- merge(
#   x = merge(x = df.result_ASDE.trt.cis,y = df.result_ASDE.trt.total,by = "dmel_gene",all = T,suffixes = c(".cis",".total"),),
#   y = df.result_ASDE.trt.trans%>%`colnames<-`(paste0(colnames(.),".trans")),
#   by.x = "dmel_gene",by.y = "dmel_gene.trans",all = T
# )%>%mutate(treatment="Wounded")
# 
# result.ASDE.reac <- merge(
#   x = merge(x = df.result_ASDE.reac.cis,y = df.result_ASDE.reac.total,by = "dmel_gene",all = T,suffixes = c(".cis",".total"),),
#   y = df.result_ASDE.reac.trans%>%`colnames<-`(paste0(colnames(.),".trans")),
#   by.x = "dmel_gene",by.y = "dmel_gene.trans",all = T
# )%>%mutate(treatment="Response")
# 
# 
# for (df in c("result.ASDE.ctrl","result.ASDE.reac","result.ASDE.trt")) {
#   .GlobalEnv[[df]]$Category<-""
#   .GlobalEnv[[df]]$Category.factor<-""
#   
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis<0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Cis only"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Cis only"]<-1
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis<0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Ambiguous_cis.but.no.total"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_cis.but.no.total"]<-1.1
#   
#   # .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
#   #                           &(.GlobalEnv[[df]]$FDR.cis>0.05)
#   #                           &(abs(.GlobalEnv[[df]]$logFC.cis-.GlobalEnv[[df]]$logFC)-
#   #                               abs(.GlobalEnv[[df]]$logFC.trans-.GlobalEnv[[df]]$logFC)<=-1)
#   #                           &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Ambiguous_total.but.no.cis.no.trans.yet.larger.cis.contribution"
#   # .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.cis.no.trans.yet.large.cisFC"]<-1.2 # large FC
#   
#   
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis>0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans<0.05)]<-"Trans only"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Trans only"]<-2
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis>0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans<0.05)]<-"Ambiguous_trans.but.no.total"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_trans.but.no.total"]<-2.1
#   
#   # .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
#   #                           &(.GlobalEnv[[df]]$FDR.cis>0.05)
#   #                           &(.GlobalEnv[[df]]$FDR.trans>0.05)
#   #                           &(abs(.GlobalEnv[[df]]$logFC.trans-.GlobalEnv[[df]]$logFC.total)-
#   #                               abs(.GlobalEnv[[df]]$logFC.cis-.GlobalEnv[[df]]$logFC.total)<=-1)]<-"Ambiguous_total.but.no.trans.no.cis.yet.larger.trans.contribution"
#   # .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.trans.no.cis.yet.large.transFC"]<-2.2 # large FC
#   
#   
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis<0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans<0.05)
#                             &(.GlobalEnv[[df]]$logFC.cis*.GlobalEnv[[df]]$logFC.trans<0)]<-"CisxTrans"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="CisxTrans"]<-3
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis<0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans<0.05)
#                             &(.GlobalEnv[[df]]$logFC.cis*.GlobalEnv[[df]]$logFC.trans>0)]<-"Cis+Trans"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Cis+Trans"]<-4
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis<0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans<0.05)]<-"Compensatory"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Compensatory"]<-5
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis>0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Conserved"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Conserved"]<-7
#   
#   .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
#                             &(.GlobalEnv[[df]]$FDR.cis>0.05)
#                             &(.GlobalEnv[[df]]$FDR.trans>0.05)
#                             &(.GlobalEnv[[df]]$Category=='')]<-"Ambiguous_total.but.no.cis.or.trans"
#   .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.cis.or.trans"]<-6 # large FC
#   
#   
#   .GlobalEnv[[df]]$Category.factor<-as.numeric(.GlobalEnv[[df]]$Category.factor)
#   .GlobalEnv[[df]]$Cat.main.factor <- round(.GlobalEnv[[df]]$Category.factor)
#   .GlobalEnv[[df]]$Cat.main <- .GlobalEnv[[df]]$Category
#   .GlobalEnv[[df]]$Cat.main[.GlobalEnv[[df]]$Cat.main.factor == 1] <- '*Cis* related'
#   .GlobalEnv[[df]]$Cat.main[.GlobalEnv[[df]]$Cat.main.factor == 2] <- '*Trans* related'
#   .GlobalEnv[[df]]$Cat.main[.GlobalEnv[[df]]$Cat.main.factor == 6] <- 'Ambiguous'
#   .GlobalEnv[[df]]$sub.cat <- .GlobalEnv[[df]]$Category
#   .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==1.1] <- "Cis*"
#   .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==1.2] <- "Cis**"
#   .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==2.1] <- "Trans*"
#   .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==2.2] <- "Trans**"
#   .GlobalEnv[[df]]$sub.cat[.GlobalEnv[[df]]$Category.factor==6] <- "Ambiguous"
#   .GlobalEnv[[df]]$cat2 <- NA
#   .GlobalEnv[[df]]$cat2[.GlobalEnv[[df]]$Category.factor%in%c(1,1.1,1.2,3,4,5)] <- "*Cis*"
#   .GlobalEnv[[df]] <- bind_rows(.GlobalEnv[[df]],(.GlobalEnv[[df]])%>%filter(Category.factor%in%c(3,4,5))%>%mutate(cat2="*Trans*"))
#   .GlobalEnv[[df]]$cat2[.GlobalEnv[[df]]$Category.factor%in%c(2,2.1,2.2)] <- "*Trans*"
#   .GlobalEnv[[df]]$cat2[.GlobalEnv[[df]]$Category.factor==6] <- "Ambiguous"
#   .GlobalEnv[[df]]$cat2.col <- NA
#   .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor%in%c(1,1.1,1.2)] <- "*Cis* divergent"
#   .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor%in%c(2,2.1,2.2)] <- "*Trans* divergent"
#   .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor%in%c(3,4,5)] <- "Both *Cis* and *Trans* divergent"
#   .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor==6] <- "Ambiguous"
#   
#   # .GlobalEnv[[df]]$dmel_gene<-factor(x = .GlobalEnv[[df]]$dmel_gene,
#   #                                    levels = .GlobalEnv[[df]]$dmel_gene[order(.GlobalEnv[[df]]$Category.factor,decreasing = F)])
# }
# 
# 
# 
# 
# plot.data <- bind_rows(result.ASDE.ctrl,result.ASDE.trt,result.ASDE.reac)



immune.genes<-read.table("immune.responsive.genes.txt")























































#############################################################
#### hybrid allele test
#############################################################
# hyb.s.gene <- read.table("genes.sig.in.any.logfc1.5.txt")
hyb.s.gene<-data.frame(V1=immune.genes.logfc1)
hyb.gene <- unique(plot.data$dmel_gene[!plot.data$dmel_gene%in%hyb.s.gene$V1])
f1.test.data <- plot.data%>%filter(dmel_gene%in%hyb.s.gene$V1,treatment=="Control")
hyb.s.gene<-merge(x = hyb.s.gene,y = plot.data%>%filter(treatment=="Control")%>%select(dmel_gene,logFC.cis),by.y = "dmel_gene",by.x = "V1",all.x = T )%>%distinct()
binom.test(x = matrix(c(length(hyb.s.gene$V1[hyb.s.gene$logFC.cis>0]),length(hyb.s.gene$V1[hyb.s.gene$logFC.cis<0])),
                      nrow = 1,ncol = 2),p = 0.5)

hyb.test <- data.frame(dmel=c(length(hyb.s.gene$V1[hyb.s.gene$logFC.cis>0]),length((result.ASDE.reac%>%filter(dmel_gene%in%hyb.gene,logFC.cis>0))$dmel_gene)),
                       dsim=c(length(hyb.s.gene$V1[hyb.s.gene$logFC.cis<0]),length((result.ASDE.reac%>%filter(dmel_gene%in%hyb.gene,logFC.cis<0))$dmel_gene)))%>%`rownames<-`(.,c("s.genes","whole"))
chisq.test(hyb.test)
fisher.test(hyb.test)

hyb.test.1<-bind_rows(
  tibble(set=rep("s.gene",hyb.test$dmel[1]),allele=rep("dmel",hyb.test$dmel[1])),
  tibble(set=rep("s.gene",hyb.test$dsim[1]),allele=rep("dsim",hyb.test$dsim[1])),
  tibble(set=rep("whole",hyb.test$dmel[2]),allele=rep("dmel",hyb.test$dmel[2])),
  tibble(set=rep("whole",hyb.test$dsim[2]),allele=rep("dsim",hyb.test$dsim[2]))
)
glm.fit <- glm(formula = allele ~ set,data = hyb.test.1%>%mutate(allele=factor(allele,levels = c("dsim","dmel")),set=factor(set,levels = c("whole","s.gene"))),family = "binomial")
summary(glm.fit)
confint(glm.fit)
glm.fit.1 <- glm(formula = allele ~ set,data = hyb.test.1%>%mutate(allele=factor(allele,levels = c("dsim","dmel")),set=factor(set,levels = c("s.gene","whole"))),family = "binomial")
summary(glm.fit.1)
confint(glm.fit.1)

# f1.est <- read.csv("f1.test.est_1.5.csv")
f1.est <- read.csv("f1.test.est_1.csv",header = T)

pdf("f1.test.logit.pdf",width = 7,height = 5)
ggplot(data = f1.est%>%mutate(treatment=factor(treatment,levels=c("Control","Wounded")))%>%filter(treatment!="Response"),aes(x = X,y = est))+
  geom_hline(yintercept = 0,linetype=3,colour="red")+
  geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.1)+
  facet_grid(cols = vars(treatment))+
  ylab("Logit(*D. melanogaster* allele)")+xlab("")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0,'lines'),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position=c(0.92, 0.9),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text.x = element_text(size=20,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size=20),
        strip.background = element_rect(fill=F,colour='black'),
        strip.text.x = element_markdown(size=20))
dev.off()

pdf("f1.test.hist.pdf",width = 7,height = 5)
ggplot(data = plot.data%>%filter(dmel_gene%in%hyb.s.gene$V1,treatment!="Response")%>%mutate(treatment=factor(treatment,levels=c("Control","Wounded"))))+
  geom_histogram(aes(x=logFC.cis),bins=17)+
  geom_vline(xintercept = 0,linetype = 5,colour="black")+
  facet_grid(cols = vars(treatment),scales = 'free')+
  labs(y='Counts',
       x='Log<sub>2</sub>(*D. melanogaster* / *D. simulans*) in hybrids')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0,'lines'),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position=c(0.92, 0.9),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size=20),
        strip.background = element_rect(fill=F,colour='black'),
        strip.text.x = element_markdown(size=20))
dev.off()



