plot.colours<-c(
  "*Cis* diverged"="#000000",
  "*Trans* diverged"="#C7050F",
  "Both *Cis* and *Trans* diverged"="gold",
  "Ambiguous"="#0DD9E7",
  "Conserved"="gray90"
)

### arrange results:
result.ASDE.ctrl <- merge(
  x = merge(x = df.result_ASDE.ctrl.cis,y = df.result_ASDE.ctrl.total,by = "dmel_gene",all = T,suffixes = c(".cis",".total"),),
  y = df.result_ASDE.ctrl.trans%>%`colnames<-`(paste0(colnames(.),".trans")),
  by.x = "dmel_gene",by.y = "dmel_gene.trans",all = T
)%>%mutate(treatment="Control")

result.ASDE.trt <- merge(
  x = merge(x = df.result_ASDE.trt.cis,y = df.result_ASDE.trt.total,by = "dmel_gene",all = T,suffixes = c(".cis",".total"),),
  y = df.result_ASDE.trt.trans%>%`colnames<-`(paste0(colnames(.),".trans")),
  by.x = "dmel_gene",by.y = "dmel_gene.trans",all = T
)%>%mutate(treatment="Wounded")

result.ASDE.reac <- merge(
  x = merge(x = df.result_ASDE.reac.cis,y = df.result_ASDE.reac.total,by = "dmel_gene",all = T,suffixes = c(".cis",".total"),),
  y = df.result_ASDE.reac.trans%>%`colnames<-`(paste0(colnames(.),".trans")),
  by.x = "dmel_gene",by.y = "dmel_gene.trans",all = T
)%>%mutate(treatment="Response")


for (df in c("result.ASDE.ctrl","result.ASDE.reac","result.ASDE.trt")) {
  .GlobalEnv[[df]]$Category<-NA
  .GlobalEnv[[df]]$Category.factor<-NA
  
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Cis only"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Cis only"]<-1
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Ambiguous_cis.but.no.total"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_cis.but.no.total"]<-1.1
  
  # .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
  #                           &(.GlobalEnv[[df]]$FDR.cis>0.05)
  #                           &(abs(.GlobalEnv[[df]]$logFC.cis-.GlobalEnv[[df]]$logFC)-
  #                               abs(.GlobalEnv[[df]]$logFC.trans-.GlobalEnv[[df]]$logFC)<=-1)
  #                           &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Ambiguous_total.but.no.cis.no.trans.yet.larger.cis.contribution"
  # .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.cis.no.trans.yet.large.cisFC"]<-1.2 # large FC
  
  
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<0.05)]<-"Trans only"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Trans only"]<-2
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<0.05)]<-"Ambiguous_trans.but.no.total"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_trans.but.no.total"]<-2.1
  
  # .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
  #                           &(.GlobalEnv[[df]]$FDR.cis>0.05)
  #                           &(.GlobalEnv[[df]]$FDR.trans>0.05)
  #                           &(abs(.GlobalEnv[[df]]$logFC.trans-.GlobalEnv[[df]]$logFC.total)-
  #                               abs(.GlobalEnv[[df]]$logFC.cis-.GlobalEnv[[df]]$logFC.total)<=-1)]<-"Ambiguous_total.but.no.trans.no.cis.yet.larger.trans.contribution"
  # .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.trans.no.cis.yet.large.transFC"]<-2.2 # large FC
  
  
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<0.05)
                            &(.GlobalEnv[[df]]$logFC.cis*.GlobalEnv[[df]]$logFC.trans<0)]<-"CisxTrans"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="CisxTrans"]<-3
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<0.05)
                            &(.GlobalEnv[[df]]$logFC.cis*.GlobalEnv[[df]]$logFC.trans>0)]<-"Cis+Trans"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Cis+Trans"]<-4
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis<0.05)
                            &(.GlobalEnv[[df]]$FDR.trans<0.05)]<-"Compensatory"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Compensatory"]<-5
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total>0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)]<-"Conserved"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Conserved"]<-7
  
  .GlobalEnv[[df]]$Category[(.GlobalEnv[[df]]$FDR.total<0.05)
                            &(.GlobalEnv[[df]]$FDR.cis>0.05)
                            &(.GlobalEnv[[df]]$FDR.trans>0.05)
                            &(is.na(.GlobalEnv[[df]]$Category))]<-"Ambiguous_total.but.no.cis.or.trans"
  .GlobalEnv[[df]]$Category.factor[.GlobalEnv[[df]]$Category=="Ambiguous_total.but.no.cis.or.trans"]<-6 # large FC
  
  
 
  .GlobalEnv[[df]]$cat2.col <- NA
  .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor%in%c(1,1.1,1.2)] <- "*Cis* diverged"
  .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor%in%c(2,2.1,2.2)] <- "*Trans* diverged"
  .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor%in%c(3,4,5)] <- "Both *Cis* and *Trans* diverged"
  .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor==6] <- "Ambiguous"
  .GlobalEnv[[df]]$cat2.col[.GlobalEnv[[df]]$Category.factor==7] <- "Conserved"
  .GlobalEnv[[df]]$cat2 <- NA
  .GlobalEnv[[df]] <- bind_rows(.GlobalEnv[[df]],(.GlobalEnv[[df]])%>%filter(cat2.col=="Both *Cis* and *Trans* diverged")%>%mutate(cat2="*Trans*"))
  .GlobalEnv[[df]]$cat2[.GlobalEnv[[df]]$Category.factor%in%c(1,1.1,1.2)] <- "*Cis*"
  .GlobalEnv[[df]]$cat2[.GlobalEnv[[df]]$Category.factor%in%c(2,2.1,2.2)] <- "*Trans*"
  .GlobalEnv[[df]]$cat2[(.GlobalEnv[[df]]$Category.factor%in%c(3,4,5))&is.na(.GlobalEnv[[df]]$cat2)] <- "*Cis*"
  .GlobalEnv[[df]]$cat2[.GlobalEnv[[df]]$Category.factor==6] <- "Ambiguous"
  .GlobalEnv[[df]]$cat2[.GlobalEnv[[df]]$Category.factor==7] <- "Conserved"
  
}




plot.data <- bind_rows(result.ASDE.ctrl,result.ASDE.trt,result.ASDE.reac)%>%mutate(treatment=factor(treatment,levels = c("Control","Wounded","Response")))

################################
##### transcriptome ##
################################
plot.whole.1 <- ggplot()+
  scale_colour_manual(values = plot.colours)+
  geom_abline(slope = 1,intercept = 0,linetype = 3)+geom_abline(slope = 0,intercept = 0,linetype = 3)+geom_vline(xintercept = 0,linetype = 3)+
  annotate(geom = "text",x = 14.5,y = 0.5,label = "italic(trans)~only",parse=T,colour = "red",fontface = "bold", size = 5)+
  annotate(geom = "text",x = 14.5,y = 16,label = "italic(cis)~only",parse=T,colour = "black", fontface = "bold",size = 5, angle = 45)+
  geom_point(data = plot.data%>%filter(cat2.col=="Conserved",treatment=="Control"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="Ambiguous",treatment=="Control"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="Both *Cis* and *Trans* diverged",treatment=="Control"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="*Cis* diverged",treatment=="Control"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="*Trans* diverged",treatment=="Control"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  labs(
    x="Log<sub>2</sub>(Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)",
    y="Log<sub>2</sub>(Hybrid<sub>*mel*</sub>/Hybrid<sub>*sim*</sub>)")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size = 20),
        legend.position = "none",
        panel.grid = element_blank())+
  xlim(-17,17)+ylim(-17,17)

plot.whole.2<-ggplot(data = plot.data%>%filter(Category!='Conserved',treatment=="Control")%>%mutate(cat2=factor(cat2,levels = c("Conserved","Ambiguous","*Trans*","*Cis*"))),
                     mapping = aes())+
  scale_fill_manual(values = plot.colours)+
  geom_bar(aes(y = cat2, fill = cat2.col),position = position_stack(reverse = F))+
  labs(
    x="Number of genes",
    y="")+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text.y = element_markdown(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_markdown(size=20),
        # plot.margin=unit(c(2,0.5,0.5,0.5),"cm"),
        legend.position=c(-0.1, 1.2),
        legend.title = element_blank(),
        panel.grid = element_blank())

plot.data.hist <- bind_rows(
  plot.data%>%filter(cat2=="*Cis*")%>%select(dmel_gene,logFC.cis,cat2,treatment)%>%mutate(logFC=abs(logFC.cis))%>%select(-logFC.cis),
  plot.data%>%filter(cat2=="*Trans*")%>%select(dmel_gene,logFC.trans,cat2,treatment)%>%mutate(logFC=abs(logFC.trans))%>%select(-logFC.trans),
)
plot.whole.3<-ggplot( data=plot.data.hist%>%filter(treatment=="Control"), aes(x=logFC, fill=cat2)) +
  geom_histogram( color="black" ,position = position_dodge(),bins = 17) +
  scale_fill_manual(values=list("*Cis*" = "#000000","*Trans*" = "#C7050F")) +
  scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
                labels = trans_format(trans = "log10",format = math_format(.x)),
                limits=c(NA,10000),
                oob=squish_infinite)+
  labs(y='Log<sub>10</sub> Number of Genes',
       x='Magnitude of expression divergence (log<sub>2</sub>fold change)')+
  theme_bw()+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0,'lines'),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        legend.position=c(0.8, 0.8),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size=20),
        strip.background = element_rect(fill=F,colour='black'),
        strip.text.x = element_markdown(size=20))



pdf("cis.trans.plot_whole.pdf",height = 8,width = 15)
ggarrange(ggarrange(NULL,plot.whole.1,nrow = 2,heights = c(0.01,1)),
          NULL,
          ggarrange(plot.whole.2,plot.whole.3,nrow = 2,labels = c("B","C"),hjust = 1,vjust = 1,font.label = list(family="sans",face="plain",size=20)),
          ncol = 3,labels = c(NULL,"A"),vjust = 1,widths = c(1,0.1,1),font.label = list(family="sans",face="plain",size=20))+
  theme(plot.margin = margin(2,0.4,0.1,0.1, "cm"))
dev.off()



pdf("cis.trans.log.plot_cat2.pdf",width = 15,height = 6.3)
ggplot()+
  scale_colour_manual(values = plot.colours)+
  geom_abline(slope = 1,intercept = 0,linetype = 3)+geom_abline(slope = 0,intercept = 0,linetype = 3)+geom_vline(xintercept = 0,linetype = 3)+
  annotate(geom = "text",x = 14.5,y = 0.5,label = "italic(trans)~only",parse=T,colour = "red",fontface = "bold", size = 5)+
  annotate(geom = "text",x = 14.5,y = 16,label = "italic(cis)~only",parse=T,colour = "black", fontface = "bold",size = 5, angle = 45)+
  geom_point(data = plot.data%>%filter(cat2.col=="Conserved"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="Ambiguous"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="Both *Cis* and *Trans* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="*Cis* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data%>%filter(cat2.col=="*Trans* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  labs(
    x="Log<sub>2</sub>(Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)",
    y="Log<sub>2</sub>(Hybrid<sub>*mel*</sub>/Hybrid<sub>*sim*</sub>)")+
  theme_bw()+
  facet_grid(cols = vars(treatment))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size = 20),
        legend.position = "none",
        panel.grid = element_blank())+
  xlim(-17,17)+ylim(-17,17)
dev.off()

pdf("cis.trans.bar.plot_cat2.pdf",width = 15,height = 2.5)
ggplot(data = plot.data%>%filter(Category!='Conserved')%>%mutate(cat2=factor(cat2,levels = c("Conserved","Ambiguous","*Trans*","*Cis*"))),
       mapping = aes())+
  scale_fill_manual(values = plot.colours)+
  geom_bar(aes(y = cat2, fill = cat2.col),position = position_stack(reverse = F))+
  facet_grid(cols = vars(treatment))+
  labs(
    x="Number of genes",
    y="")+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text.y = element_markdown(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(1.2,0.5,0.5,0.5),"cm"),
        legend.position=c(0.55, 1.2),
        legend.title = element_blank(),
        # legend.position = "none",
        panel.grid = element_blank())
dev.off()


plot.data.hist <- bind_rows(
  plot.data%>%filter(cat2=="*Cis*")%>%select(dmel_gene,logFC.cis,cat2,treatment)%>%mutate(logFC=abs(logFC.cis))%>%select(-logFC.cis),
  plot.data%>%filter(cat2=="*Trans*")%>%select(dmel_gene,logFC.trans,cat2,treatment)%>%mutate(logFC=abs(logFC.trans))%>%select(-logFC.trans),
)
pdf('cis.trans.number.genes.hist.plot_cat2.pdf',width = 15,height = 5)
ggplot( data=plot.data.hist, aes(x=logFC, fill=cat2)) +
  geom_histogram( color="black" ,position = position_dodge(),bins = 17) +
  scale_fill_manual(values=list("*Cis*" = "#000000","*Trans*" = "#C7050F")) +
  facet_grid(cols = vars(treatment),scales = 'free')+
  scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
                labels = trans_format(trans = "log10",format = math_format(.x)),
                limits=c(NA,10000),
                oob=squish_infinite)+
  labs(y='Log<sub>10</sub> Number of Genes',
       x='Magnitude of expression divergence (log<sub>2</sub>fold change)')+
  theme_bw()+
  guides(fill = guide_legend(direction = 'horizontal'))+
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
# pdf('cis.trans.number.genes.hist.plot_cat2.pdf',width = 15,height = 5)
# ggplot(data = plot.data%>%mutate(logFC.total=abs(logFC.total))%>%filter(cat2%in%c("*Cis*","*Trans*")))+
#   geom_histogram(aes(x = logFC.total,fill = cat2),colour='black',position = position_dodge(),bins = 17)+
#   scale_fill_manual(values = list("*Cis*" = "#000000","*Trans*" = "#C7050F"))+
#   facet_grid(cols = vars(treatment))+
#   scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
#                 labels = trans_format(trans = "log10",format = math_format(.x)),
#                 limits=c(NA,10000),
#                 oob=squish_infinite)+
#   labs(y='Log<sub>10</sub> Number of Genes',
#        x='Magnitude of expression divergence (log<sub>2</sub>Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)')+
#   theme_bw()+
#   guides(fill = guide_legend(direction = 'horizontal'))+
#   theme(panel.grid = element_blank(),
#         panel.spacing = unit(0,'lines'),
#         legend.title = element_blank(),
#         legend.text = element_markdown(size=20),
#         plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
#         legend.position=c(0.92, 0.9),
#         axis.title.x = element_markdown(size=20),
#         axis.title.y = element_markdown(size=20),
#         axis.text = element_text(size=20),
#         strip.background = element_rect(fill=F,colour='black'),
#         strip.text.x = element_markdown(size=20))
# dev.off()







################################
##### immune responsive genes ##
################################
plot.data.immune <- plot.data%>%filter(dmel_gene%in%immune.genes$V1)


pdf("cis.trans.log.plot_cat2_immune.pdf",width = 15,height = 6.3)
ggplot()+
  scale_colour_manual(values = plot.colours)+
  geom_abline(slope = 1,intercept = 0,linetype = 3)+geom_abline(slope = 0,intercept = 0,linetype = 3)+geom_vline(xintercept = 0,linetype = 3)+
  annotate(geom = "text",x = 14.5,y = 0.5,label = "italic(trans)~only",parse=T,colour = "red",fontface = "bold", size = 5)+
  annotate(geom = "text",x = 14.5,y = 16,label = "italic(cis)~only",parse=T,colour = "black", fontface = "bold",size = 5, angle = 45)+
  geom_point(data = plot.data.immune%>%filter(cat2.col=="Conserved"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune%>%filter(cat2.col=="Ambiguous"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune%>%filter(cat2.col=="Both *Cis* and *Trans* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune%>%filter(cat2.col=="*Cis* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune%>%filter(cat2.col=="*Trans* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  labs(
    x="Log<sub>2</sub>(Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)",
    y="Log<sub>2</sub>(Hybrid<sub>*mel*</sub>/Hybrid<sub>*sim*</sub>)")+
  theme_bw()+
  facet_grid(cols = vars(treatment))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size = 20),
        legend.position = "none",
        panel.grid = element_blank())+
  xlim(-17,17)+ylim(-17,17)
dev.off()

pdf("cis.trans.bar.plot_cat2_immune.pdf",width = 15,height = 2.5)
ggplot(data = plot.data.immune%>%filter(Category!='Conserved')%>%mutate(cat2=factor(cat2,levels = c("Conserved","Ambiguous","*Trans*","*Cis*"))),
       mapping = aes())+
  scale_fill_manual(values = plot.colours)+
  geom_bar(aes(y = cat2, fill = cat2.col),position = position_stack(reverse = F))+
  facet_grid(cols = vars(treatment))+
  labs(
    x="Number of genes",
    y="")+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text.y = element_markdown(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(1.2,0.5,0.5,0.5),"cm"),
        legend.position=c(0.55, 1.2),
        legend.title = element_blank(),
        # legend.position = "none",
        panel.grid = element_blank())
dev.off()



plot.data.immune.hist <- bind_rows(
  plot.data.immune%>%filter(cat2=="*Cis*")%>%select(dmel_gene,logFC.cis,cat2,treatment)%>%mutate(logFC=abs(logFC.cis))%>%select(-logFC.cis),
  plot.data.immune%>%filter(cat2=="*Trans*")%>%select(dmel_gene,logFC.trans,cat2,treatment)%>%mutate(logFC=abs(logFC.trans))%>%select(-logFC.trans),
)
pdf('cis.trans.number.genes.hist.plot_cat2_immue.pdf',width = 15,height = 5)
ggplot( data=plot.data.immune.hist, aes(x=logFC, fill=cat2)) +
  geom_histogram( color="black" ,position = position_dodge(),bins = 17) +
  scale_fill_manual(values=list("*Cis*" = "#000000","*Trans*" = "#C7050F")) +
  facet_grid(cols = vars(treatment),scales = 'free')+
  scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
                labels = trans_format(trans = "log10",format = math_format(.x)),
                limits=c(NA,10000),
                oob=squish_infinite)+
  labs(y='Log<sub>10</sub> Number of Genes',
       x='Magnitude of expression divergence (log<sub>2</sub>fold change)')+
  theme_bw()+
  guides(fill = guide_legend(direction = 'horizontal'))+
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
# pdf('cis.trans.number.genes.hist.plot_cat2_immue.pdf',width = 15,height = 5)
# ggplot(data = plot.data.immune%>%mutate(logFC.total=abs(logFC.total))%>%filter(cat2%in%c("*Cis*","*Trans*")))+
#   geom_histogram(aes(x = logFC.total,fill = cat2),colour='black',position = position_dodge(),bins = 17)+
#   scale_fill_manual(values = list("*Cis*" = "#000000","*Trans*" = "#C7050F"))+
#   facet_grid(cols = vars(treatment))+
#   scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
#                 labels = trans_format(trans = "log10",format = math_format(.x)),
#                 limits=c(NA,10000),
#                 oob=squish_infinite)+
#   labs(y='Log<sub>10</sub> Number of Genes',
#        x='Magnitude of expression divergence (log<sub>2</sub>Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)')+
#   theme_bw()+
#   guides(fill = guide_legend(direction = 'horizontal'))+
#   theme(panel.grid = element_blank(),
#         panel.spacing = unit(0,'lines'),
#         legend.title = element_blank(),
#         legend.text = element_markdown(size=20),
#         plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
#         legend.position=c(0.92, 0.9),
#         axis.title.x = element_markdown(size=20),
#         axis.title.y = element_markdown(size=20),
#         axis.text = element_text(size=20),
#         strip.background = element_rect(fill=F,colour='black'),
#         strip.text.x = element_markdown(size=20))
# dev.off()


#######################################
##### immune responsive genes logfc1 ##
#######################################

immune.genes.logfc1 <- read.table("genes.sig.in.any.logfc1.txt")
immune.genes.logfc1 <- immune.genes.logfc1[!is.na(immune.genes.logfc1)]
plot.data.immune.logfc1 <- plot.data%>%filter(dmel_gene%in%immune.genes.logfc1)


# pdf("cis.trans.log.plot_cat2_immune_logfc1.pdf",width = 15,height = 6.3)
plot.cis.trans.1 <- ggplot()+
  scale_colour_manual(values = plot.colours)+
  geom_abline(slope = 1,intercept = 0,linetype = 3)+geom_abline(slope = 0,intercept = 0,linetype = 3)+geom_vline(xintercept = 0,linetype = 3)+
  annotate(geom = "text",x = 14.5,y = 0.5,label = "italic(trans)~only",parse=T,colour = "red",fontface = "bold", size = 5)+
  annotate(geom = "text",x = 14.5,y = 16,label = "italic(cis)~only",parse=T,colour = "black", fontface = "bold",size = 5, angle = 45)+
  geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="Conserved"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="Ambiguous"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="Both *Cis* and *Trans* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="*Cis* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="*Trans* diverged"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
  labs(
    x="Log<sub>2</sub>(Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)",
    y="Log<sub>2</sub>(Hybrid<sub>*mel*</sub>/Hybrid<sub>*sim*</sub>)")+
  theme_bw()+
  facet_grid(cols = vars(treatment))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size = 20),
        plot.margin=unit(c(2.5,0.5,0.5,0.5),"cm"),
        legend.position=c(0.5, 1.15),
        legend.text = element_markdown(size=20),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.direction = "horizontal")+
  guides(col = guide_legend(override.aes = list(shape = 15, size = 10)))+
  xlim(-17,17)+ylim(-17,17)
# dev.off()

# ggplot()+
#   scale_colour_manual(values = plot.colours)+
#   geom_abline(slope = 1,intercept = 0,linetype = 3)+geom_abline(slope = 0,intercept = 0,linetype = 3)+geom_vline(xintercept = 0,linetype = 3)+
#   annotate(geom = "text",x = 14.5,y = 0.5,label = "italic(trans)~only",parse=T,colour = "red",fontface = "bold", size = 5)+
#   annotate(geom = "text",x = 14.5,y = 16,label = "italic(cis)~only",parse=T,colour = "black", fontface = "bold",size = 5, angle = 45)+
#   geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="Conserved",treatment!="Response"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
#   geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="Ambiguous",treatment!="Response"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
#   geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="Both *Cis* and *Trans* diverged",treatment!="Response"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
#   geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="*Cis* diverged",treatment!="Response"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
#   geom_point(data = plot.data.immune.logfc1%>%filter(cat2.col=="*Trans* diverged",treatment!="Response"),aes(x = logFC.total,y = logFC.cis,colour = cat2.col))+
#   geom_text(data=plot.data.immune.logfc1%>%filter(treatment!="Response"),aes(x = logFC.total,y = logFC.cis,label=if_else(condition = abs(logFC.cis)>8,true = dmel_gene,false = "")))+
#   labs(
#     x="Log<sub>2</sub>(Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)",
#     y="Log<sub>2</sub>(Hybrid<sub>*mel*</sub>/Hybrid<sub>*sim*</sub>)")+
#   theme_bw()+
#   facet_grid(cols = vars(treatment))+
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 20),
#         axis.title.x = element_markdown(size=20),
#         axis.title.y = element_markdown(size=20),
#         axis.text = element_text(size = 20),
#         plot.margin=unit(c(0.5,6,0.5,0.5),"cm"),
#         legend.position=c(1.3, 0.5),
#         legend.text = element_markdown(size=20),
#         legend.title = element_blank(),
#         panel.grid = element_blank(),
#         legend.direction = "vertical")+
#   guides(col = guide_legend(override.aes = list(shape = 15, size = 10)))+
#   xlim(-17,17)+ylim(-17,17)

# pdf("cis.trans.bar.plot_cat2_immune_logfc1.pdf",width = 15,height = 2.5)
plot.cis.trans.2<-ggplot(data = plot.data.immune.logfc1%>%filter(Category!='Conserved',treatment!="Response")%>%mutate(cat2=factor(cat2,levels = c("Conserved","Ambiguous","*Trans*","*Cis*"))),
       mapping = aes())+
  scale_fill_manual(values = plot.colours)+
  geom_bar(aes(y = cat2, fill = cat2.col),position = position_stack(reverse = F))+
  facet_grid(cols = vars(treatment))+
  labs(
    x="Number of genes",
    y="")+
  # guides(fill = guide_legend(direction = 'horizontal'))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text.y = element_markdown(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_markdown(size=20),
        # plot.margin=unit(c(1.2,0.5,0.5,0.5),"cm"),
        # legend.position=c(0.55, 1.2),
        # legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
# dev.off()

# #### Contribution
# plot.data.immune.logfc1.contribution <- merge(plot.data.immune.logfc1%>%filter(treatment=="Wounded")%>%mutate(cat2=cat2.col)%>%select(-cat2.col),
#                                               plot.data.immune.logfc1%>%filter(treatment=="Control")%>%select(dmel_gene,cat2.col),by="dmel_gene")%>%distinct()
# 
# ggplot(data = plot.data.immune.logfc1.contribution%>%filter(Category!='Conserved'),
#        mapping = aes())+
#   scale_fill_manual(values = plot.colours)+
#   geom_bar(aes(y = cat2, fill = cat2.col),position = position_stack(reverse = F))+
#   labs(
#     x="Number of genes",
#     y="")+
#   guides(fill = guide_legend(direction = 'horizontal'))+
#   theme_bw()+
#   theme(strip.background = element_blank(),
#         strip.text = element_blank(),
#         axis.title.x = element_markdown(size=20),
#         axis.title.y = element_markdown(size=20),
#         axis.text.y = element_markdown(size = 20),
#         axis.text.x = element_text(size = 20),
#         legend.text = element_markdown(size=20),
#         plot.margin=unit(c(1.2,0.5,0.5,0.5),"cm"),
#         legend.position=c(0.55, 1.2),
#         legend.title = element_blank(),
#         # legend.position = "none",
#         panel.grid = element_blank())






plot.data.immune.logfc1.hist <- bind_rows(
  plot.data.immune.logfc1%>%filter(cat2=="*Cis*")%>%select(dmel_gene,logFC.cis,cat2,treatment)%>%mutate(logFC=abs(logFC.cis))%>%select(-logFC.cis),
  plot.data.immune.logfc1%>%filter(cat2=="*Trans*")%>%select(dmel_gene,logFC.trans,cat2,treatment)%>%mutate(logFC=abs(logFC.trans))%>%select(-logFC.trans),
)
# pdf('cis.trans.number.genes.hist.plot_2_immue_logfc1.pdf',width = 15,height = 5)
plot.cis.trans.3 <- ggplot( data=plot.data.immune.logfc1.hist%>%filter(treatment!="Response"), aes(x=logFC, fill=cat2)) +
  geom_histogram( color="black" ,position = position_dodge(),bins = 15) +
  scale_fill_manual(values=list("*Cis*" = "#000000","*Trans*" = "#C7050F")) +
  facet_grid(cols = vars(treatment),scales = "free")+
  scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
                labels = trans_format(trans = "log10",format = math_format(.x)),
                limits=c(NA,100),
                oob=squish_infinite)+
  labs(y='Log<sub>10</sub> Number of Genes',
       x='Magnitude of expression divergence (log<sub>2</sub>fold change)')+
  theme_bw()+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0,'lines'),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position=c(0.8, 0.8),
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size=20),
        strip.background = element_rect(fill=F,colour='black'),
        strip.text.x = element_markdown(size=20))
# dev.off()


# pdf('cis.trans.violin_cat2_immue_logfc1.pdf',width = 7,height = 7)
plot.cis.trans.4<-ggplot( data=plot.data.immune.logfc1.hist%>%filter(treatment!="Response"), aes(x=treatment, y=logFC,fill=treatment)) +
  geom_violin() +
  geom_boxplot(coef=10, width=0.1)+
  scale_fill_manual(values=list("*Cis*" = "#000000","*Trans*" = "#C7050F","Control"="white","Wounded"="gray40")) +
  facet_grid(cols = vars(cat2),scales = "free")+
  labs(x='',
       y='Magnitude of expression divergence (log<sub>2</sub>fold change)')+
  theme_bw()+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0,'lines'),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        legend.position="none",
        axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size=20),
        strip.background = element_rect(fill=F,colour='black'),
        strip.text.x = element_markdown(size=20))
# dev.off()


pdf("cid.trans_logfc1.pdf",width = 17,height = 13.3)
# ggarrange(
#   ggarrange(
#    ggarrange(NULL,plot.cis.trans.1,ncol = 2,widths = c(0.1,1)),
#     plot.cis.trans.2,
#     ggarrange(NULL,plot.cis.trans.3,ncol = 2,widths = c(0.1,1)),
#     nrow = 3,heights = c(2,0.8,1.5),labels = c("A","B","C"),vjust = c(4.5,1.5,1.5),font.label = list(face="plain",size=20)
#   ),
#   NULL,
#   ggarrange(NULL,
#             plot.cis.trans.4,
#             NULL,nrow = 3,heights = c(0.5,1,0.01),labels = "D",vjust = 5,font.label = list(face="plain",size=20)),
#   ncol = 3,widths = c(1.5,0.1,1))

ggarrange(
  ggarrange(ggarrange(NULL,plot.cis.trans.1,ncol = 2,widths = c(0.04,1)),NULL,ncol = 2,widths = c(1,0.08)),
  ggarrange(
    ggarrange(plot.cis.trans.2,
              ggarrange(NULL,plot.cis.trans.3,ncol = 2,widths = c(0.06,1)),nrow = 2,heights = c(0.6,1),labels = c("B","C"),font.label = list(face="plain",size=20)),
    NULL,
    ggarrange(plot.cis.trans.4,NULL,nrow = 2,heights =c(2,0.1),labels = "D",font.label = list(face="plain",size=20),hjust = 1.5,vjust = -1.5),ncol = 3,widths = c(2,0.1,1)
  ),nrow = 2,heights = c(1.1,1),labels = "A",font.label = list(face="plain",size=20)
)
dev.off()


# pdf('cis.trans.number.genes.hist.plot_2_immue_logfc1.pdf',width = 15,height = 5)
# ggplot(data = plot.data.immune.logfc1%>%mutate(logFC.total=abs(logFC.total))%>%filter(cat2%in%c("*Cis*","*Trans*")))+
#   geom_histogram(aes(x = logFC.total,fill = cat2),colour='black',position = position_dodge(),bins = 25)+
#   scale_fill_manual(values = list("*Cis*" = "#000000","*Trans*" = "#C7050F"))+
#   facet_grid(cols = vars(treatment))+
#   scale_y_log10(breaks = trans_breaks(trans = "log10",inv = function(x) 10^x,n = 5),
#                 labels = trans_format(trans = "log10",format = math_format(.x)),
#                 limits=c(NA,10000),
#                 oob=squish_infinite)+
#   labs(y='Log<sub>10</sub> Number of Genes',
#        x='Magnitude of expression divergence (log<sub>2</sub>Parent<sub>*mel*</sub>/Parent<sub>*sim*</sub>)')+
#   theme_bw()+
#   guides(fill = guide_legend(direction = 'horizontal'))+
#   theme(panel.grid = element_blank(),
#         panel.spacing = unit(0,'lines'),
#         legend.title = element_blank(),
#         legend.text = element_markdown(size=20),
#         plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
#         legend.position=c(0.92, 0.9),
#         axis.title.x = element_markdown(size=20),
#         axis.title.y = element_markdown(size=20),
#         axis.text = element_text(size=20),
#         strip.background = element_rect(fill=F,colour='black'),
#         strip.text.x = element_markdown(size=20))
# dev.off()



save(plot.data.immune.logfc1.hist,plot.data.immune.logfc1,file = "ASDE.RData")




##################################
#### divergence bar plot  ########
##################################
barplot.immune.genes.0 <- read.csv("immune.plot.genes.csv",header = T)

barplot.immune.genes.tr <- plot.data.immune.logfc1%>%
  select(dmel_gene,logFC.trans,treatment,Category)%>%
  filter(dmel_gene%in%barplot.immune.genes.0$dmel_gene)%>%
  mutate(divergence="*Trans*", row="*Cis* and *Trans*",
         Category = gsub("Cis only","",Category),
         Category = gsub("Trans only","*",Category),
         Category = gsub("Ambiguous_cis.but.no.total","",Category),
         Category = gsub("Ambiguous_total.but.no.cis.or.trans","",Category),
         Category = gsub("Ambiguous_trans.but.no.total","*",Category),
         Category = gsub("Conserved","",Category),
         Category = gsub("Compensatory","*",Category),
         Category = gsub("CisxTrans","*",Category),
         Category = gsub("Cis\\+Trans","*",Category))%>%
  `colnames<-`(.,gsub(pattern = "logFC.trans",replacement = "logFC",x = colnames(.)))%>%distinct()

barplot.immune.genes.ci <- plot.data.immune.logfc1%>%
  select(dmel_gene,logFC.cis,treatment,Category)%>%
  filter(dmel_gene%in%barplot.immune.genes.0$dmel_gene)%>%
  mutate(divergence="*Cis*", row="*Cis* and *Trans*",
         Category = gsub("Cis only","*",Category),
         Category = gsub("Trans only","",Category),
         Category = gsub("Ambiguous_cis.but.no.total","*",Category),
         Category = gsub("Ambiguous_total.but.no.cis.or.trans","",Category),
         Category = gsub("Ambiguous_trans.but.no.total","",Category),
         Category = gsub("Conserved","",Category),
         Category = gsub("Compensatory","*",Category),
         Category = gsub("CisxTrans","*",Category),
         Category = gsub("Cis\\+Trans","*",Category))%>%
  `colnames<-`(.,gsub(pattern = "logFC.cis",replacement = "logFC",x = colnames(.)))%>%distinct()

barplot.immune.genes.to <- plot.data.immune.logfc1%>%
  select(dmel_gene,logFC.total,treatment,Category)%>%
  filter(dmel_gene%in%barplot.immune.genes.0$dmel_gene)%>%
  mutate(divergence="Total", row="Total",
         Category = gsub("Cis only","*",Category),
         Category = gsub("Trans only","*",Category),
         Category = gsub("Ambiguous_cis.but.no.total","",Category),
         Category = gsub("Ambiguous_total.but.no.cis.or.trans","*",Category),
         Category = gsub("Ambiguous_trans.but.no.total","*",Category),
         Category = gsub("Conserved","",Category),
         Category = gsub("Compensatory","",Category),
         Category = gsub("CisxTrans","*",Category),
         Category = gsub("Cis\\+Trans","*",Category))%>%
  `colnames<-`(.,gsub(pattern = "logFC.total",replacement = "logFC",x = colnames(.)))%>%distinct()

barplot.immune.genes <- bind_rows(barplot.immune.genes.ci,barplot.immune.genes.tr,barplot.immune.genes.to)%>%
  mutate(dmel_gene=factor(dmel_gene, levels = c(unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Gram-negative bacteria"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Gram-positive bacteria"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Fungi"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Bomanin"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Edin"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="TEP"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Iron binding"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Carbohydrate binding"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Pathogen recognition"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Regulate Toll"]),
                                                unique(barplot.immune.genes.0$dmel_gene[barplot.immune.genes.0$cat2=="Melanisation"]))))
barplot.immune.genes$ysig <- barplot.immune.genes$logFC
for (i in 1:length(barplot.immune.genes$ysig)) {
  if (barplot.immune.genes$ysig[i]>=0) {
    barplot.immune.genes$ysig[i] <- barplot.immune.genes$ysig[i] + 0.3
  }
  else
    barplot.immune.genes$ysig[i] <- 0.3
}

pdf("divergence.barplot.control.pdf",width = 16,height = 10)
ggplot(data = barplot.immune.genes%>%filter(treatment=="Control"))+
  geom_col(aes(x=dmel_gene,y = logFC,fill=divergence,colour=divergence),position = "dodge2",width = 0.65)+
  geom_hline(yintercept = 0,colour="dark gray")+
  geom_text(aes(x=dmel_gene,y = ysig,group=divergence, colour=divergence, label =Category), 
            position = position_dodge(width = 0.7),size=5,show.legend = FALSE)+
  scale_fill_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="white"))+
  scale_colour_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="black"))+
  # facet_grid(rows = vars(row))+
  theme_bw()+
  ggtitle("Control")+
  ylab("Log<sub>2</sub>(*D. melanogaster* / *D. simulans*)")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(2,0.5,0.5,0.5),"cm"),
        legend.position=c(0.85, 1),
        legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=20),
        axis.text.y  = element_text(size=20),
        axis.text.x  = element_text(size=20,angle = 45,hjust = 1),
        strip.background = element_rect(fill=F,colour='black'),strip.text.y = element_markdown(size=20),
        plot.title = element_text(color="black", size=20, face="plain",hjust = 0.5))
dev.off()






pdf("divergence.barplot.wounded.pdf",width = 16,height = 10)
ggplot(data = barplot.immune.genes%>%filter(treatment=="Wounded"))+
  geom_col(aes(x=dmel_gene,y = logFC,fill=divergence,colour=divergence),position = "dodge2",width = 0.65)+
  geom_hline(yintercept = 0,colour="dark gray")+
  geom_text(aes(x=dmel_gene,y = ysig,group=divergence, colour=divergence, label =Category), 
            position = position_dodge(width = 0.7),size=5,show.legend = FALSE)+
  scale_fill_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="white"))+
  scale_colour_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="black"))+
  # facet_grid(rows = vars(row))+
  theme_bw()+
  ggtitle("Wounded")+
  ylab("Log<sub>2</sub>(*D. melanogaster* / *D. simulans*)")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=20),
        axis.text.y  = element_text(size=20),
        axis.text.x  = element_text(size=20,angle = 45,hjust = 1),
        strip.background = element_rect(fill=F,colour='black'),strip.text.y = element_markdown(size=20),
        plot.title = element_text(color="black", size=20, face="plain",hjust = 0.5))
dev.off()





pdf("divergence.barplot.control_wounded.pdf",width = 18,height = 20)
ggplot(data = barplot.immune.genes%>%filter(treatment!="Response"))+
  geom_col(aes(x=dmel_gene,y = logFC,fill=divergence,colour=divergence),position = "dodge2",width = 0.65)+
  geom_hline(yintercept = 0,colour="dark gray")+
  geom_text(aes(x=dmel_gene,y = ysig,group=divergence, colour=divergence, label =Category), 
            position = position_dodge(width = 0.7),size=5,show.legend = FALSE)+
  scale_fill_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="white"))+
  scale_colour_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="black"))+
  facet_grid(rows = vars(treatment),scales = "free_x")+
  theme_bw()+
  ylab("Log<sub>2</sub>(*D. melanogaster* / *D. simulans*)")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(2,0.5,0.5,0.5),"cm"),
        legend.position=c(0.85, 1),
        legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=20),
        axis.text.y  = element_text(size=20),
        axis.text.x  = element_text(size=20,angle = 45,hjust = 1),
        strip.background = element_rect(fill=F,colour='black'),strip.text.y = element_markdown(size=20),
        plot.title = element_text(color="black", size=20, face="plain",hjust = 0.5))
dev.off()












pdf("divergence.barplot.response.pdf",width = 16,height = 10)
ggplot(data = barplot.immune.genes%>%filter(treatment=="Response"))+
  geom_col(aes(x=dmel_gene,y = logFC,fill=divergence,colour=divergence),position = "dodge2",width = 0.7)+
  geom_hline(yintercept = 0,colour="dark gray")+
  geom_text(aes(x=dmel_gene,y = ysig,group=divergence, colour=divergence, label =Category), 
            position = position_dodge(width = 0.8),size=5,show.legend = FALSE)+
  scale_fill_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="white"))+
  scale_colour_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="black"))+
  facet_grid(rows = vars(row))+
  theme_bw()+
  ggtitle("Response")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=20),
        axis.text.y  = element_text(size=20),
        axis.text.x  = element_text(size=20,angle = 45,hjust = 1),
        strip.background = element_rect(fill=F,colour='black'),strip.text.y = element_markdown(size=20),
        plot.title = element_text(color="black", size=20, face="plain",hjust = 0.5))+
  ylim(-10,10)
dev.off()



ggplot(data = barplot.immune.genes%>%filter(treatment!="Response"),
       aes(x= treatment,y = logFC,colour = divergence,group = divergence))+
  scale_colour_manual(values = c("*Cis*"="#000000","*Trans*"="#C7050F","Total"="blue"))+
  geom_line()+
  geom_point()+
  facet_wrap(facets = ~dmel_gene,ncol = 4,scales = "free")
  










