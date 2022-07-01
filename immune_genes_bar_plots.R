





plot.immune.genes <- read.csv("immune.plot.genes.csv",header = T)
plot.immune.genes <- merge.data.frame(
  x = plot.immune.genes,
  y = result.within.DE%>%select(dmel_gene,logFC,species),
  by = "dmel_gene",all.x = T
 )

pdf("immune.genes.barplot.1.pdf",width = 14,height = 4)
ggplot(data = plot.immune.genes%>%filter(cat1=="AMP")%>%
         mutate(dmel_gene=factor(dmel_gene, levels = c(unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Gram-negative bacteria"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Gram-positive bacteria"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Fungi"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Bomanin"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Edin"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="TEP"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Iron binding"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Carbohydrate binding"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Pathogen recognition"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Regulate Toll"]),
                                                  unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Melanisation"])))))+
  geom_col(aes(x=dmel_gene,y=logFC,fill=species),position = position_dodge())+
  scale_fill_manual(values = list("*D. melanogaster*" = "#BA4C5B","*D. simulans*"="#2970C2","Hybrid"="#6D9C76"))+
  geom_linerange(aes(x=NULL,xmin=1,xmax=8,y=5.5))+annotate(geom = "text",size=5,x = 5,y = 6,label="Gram-negative bacteria")+
  geom_linerange(aes(x=NULL,xmin=9,xmax=12,y=5.5))+annotate(geom = "text",size=5,x = 10.5,y = 6,label="Gram-positive bacteria")+
  geom_linerange(aes(x=NULL,xmin=13,xmax=15,y=5.5))+annotate(geom = "text",size=5,x = 14,y = 6,label="Fungi")+
  geom_linerange(aes(x=NULL,xmin=16,xmax=20,y=5.5))+annotate(geom = "text",size=5,x = 18,y = 6,label="Bomanins")+
  geom_linerange(aes(x=NULL,xmin=20.5,xmax=21.5,y=5.5))+annotate(geom = "text",size=5,x = 21,y = 6,label="Edin")+
  theme_bw()+
  guides(fill = guide_legend(direction = 'horizontal'))+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(2,0.5,0.5,0.5),"cm"),
        legend.position=c(0.8, 1.3),
        # legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=20),
        axis.text.y  = element_text(size=20),
        axis.text.x  = element_text(size=20,angle = 45,hjust = 1),
        strip.background = element_rect(fill=F,colour='black'))
dev.off()  



pdf("immune.genes.barplot.2.pdf",width = 8.5,height = 4)
ggplot(data = plot.immune.genes%>%filter(cat1==2)%>%
         mutate(dmel_gene=factor(dmel_gene, levels = c(unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Iron binding"]),
                                                       unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="TEP"]),
                                                       unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Carbohydrate binding"]),
                                                       unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Pathogen recognition"])))))+
  geom_col(aes(x=dmel_gene,y=logFC,fill=species),position = position_dodge())+
  scale_fill_manual(values = list("*D. melanogaster*" = "#BA4C5B","*D. simulans*"="#2970C2","Hybrid"="#6D9C76"))+
  geom_linerange(aes(x=NULL,xmin=0.5,xmax=1.5,y=5.5))+annotate(geom = "text",size=5,x = 1.2,y = 6,label="Iron binding")+
  geom_linerange(aes(x=NULL,xmin=1.6,xmax=3.5,y=5.5))+annotate(geom = "text",size=5,x = 2.5,y = 6,label="TEP")+
  geom_linerange(aes(x=NULL,xmin=3.6,xmax=5.5,y=5.5))+annotate(geom = "text",size=5,x = 4.6,y = 6,label="Carbohydrate binding")+
  geom_linerange(aes(x=NULL,xmin=5.6,xmax=8.5,y=5.5))+annotate(geom = "text",size=5,x = 7,y = 6,label="Pathogen recognition")+
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=20),
        axis.text.y  = element_text(size=20),
        axis.text.x  = element_text(size=20,angle = 45,hjust = 1),
        strip.background = element_rect(fill=F,colour='black'))
dev.off()  









pdf("immune.genes.barplot.3.pdf",width = 5.5,height = 4)
ggplot(data = plot.immune.genes%>%filter(cat1=="other")%>%
         mutate(dmel_gene=factor(dmel_gene, levels = c(unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Regulate Toll"]),
                                                       unique(plot.immune.genes$dmel_gene[plot.immune.genes$cat2=="Melanisation"])))))+
  geom_col(aes(x=dmel_gene,y=logFC,fill=species),position = position_dodge())+
  scale_fill_manual(values = list("*D. melanogaster*" = "#BA4C5B","*D. simulans*"="#2970C2","Hybrid"="#6D9C76"))+
  geom_linerange(aes(x=NULL,xmin=0.5,xmax=2.4,y=2.5))+annotate(geom = "text",size=5,x = 1.5,y = 2.7,label="Toll")+
  geom_linerange(aes(x=NULL,xmin=2.6,xmax=4.5,y=2.5))+annotate(geom = "text",size=5,x = 3.5,y =2.7,label="Melanisation")+
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(),
        legend.title = element_blank(),
        legend.text = element_markdown(size=20),
        # plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        # legend.position = "right",
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=20),
        axis.text.y  = element_text(size=20),
        axis.text.x  = element_text(size=20,angle = 45,hjust = 1),
        strip.background = element_rect(fill=F,colour='black'))
dev.off()  





  
  
  
  
  
  
  
  
  
  
