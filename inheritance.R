

plot.colours<-c(
  '*D. melanogaster* dominant' = "#F8766D",
  '*D. simulans* dominant' = "#00BFC4",
  "Additive" = "#FFC300",
  "Overdominant" = "#089911",
  "Underdominant" = "#9B5DE5",
  "Conserved" = "#778DA9"
  
)

inherit.data <- bind_rows(
  merge.data.frame(
    x = df.result_DE.hyb.mel.0%>%mutate(logFC_f1.mel=logFC)%>%select(logFC_f1.mel,dmel_gene),
    y = df.result_DE.hyb.sim.0%>%mutate(logFC_f1.sim=logFC)%>%select(logFC_f1.sim,dmel_gene), by="dmel_gene"
  )%>%mutate(treatment="Control"),
  merge.data.frame(
    x = df.result_DE.hyb.mel.1%>%mutate(logFC_f1.mel=logFC)%>%select(logFC_f1.mel,dmel_gene),
    y = df.result_DE.hyb.sim.1%>%mutate(logFC_f1.sim=logFC)%>%select(logFC_f1.sim,dmel_gene), by="dmel_gene"
  )%>%mutate(treatment="Wounded"),
  merge.data.frame(
    x = df.result_DE.hyb.mel.r%>%mutate(logFC_f1.mel=logFC)%>%select(logFC_f1.mel,dmel_gene),
    y = df.result_DE.hyb.sim.r%>%mutate(logFC_f1.sim=logFC)%>%select(logFC_f1.sim,dmel_gene), by="dmel_gene"
  )%>%mutate(treatment="Response")
)
for (i in 1:length(inherit.data$dmel_gene)) {
  if (abs(inherit.data$`logFC_f1.sim`[i])>=1.25&abs(inherit.data$`logFC_f1.mel`[i])<1) {
    inherit.data$inher[i] <- "*D. melanogaster* dominant"
  }
  if (abs(inherit.data$`logFC_f1.sim`[i])<1.25&abs(inherit.data$`logFC_f1.mel`[i])>=1) {
    inherit.data$inher[i] <- "*D. simulans* dominant"
  }
  if (abs(inherit.data$`logFC_f1.sim`[i])<1.25&abs(inherit.data$`logFC_f1.mel`[i])<1) {
    inherit.data$inher[i] <- "Conserved"
  }
  if (inherit.data$`logFC_f1.sim`[i]>=1.25&inherit.data$`logFC_f1.mel`[i]>=1) {
    inherit.data$inher[i] <- "Overdominant"
  }
  if (inherit.data$`logFC_f1.sim`[i]< -1.25&inherit.data$`logFC_f1.mel`[i]< -1) {
    inherit.data$inher[i] <- "Underdominant"
  }
  if (inherit.data$`logFC_f1.sim`[i]< -1.25&inherit.data$`logFC_f1.mel`[i] >=1) {
    inherit.data$inher[i] <- "Additive"
  }
  if (inherit.data$`logFC_f1.sim`[i]>=1.25&inherit.data$`logFC_f1.mel`[i] < -1) {
    inherit.data$inher[i] <- "Additive"
  }
}

# save(inherit.data,file = "inherit.data.RData")



pdf(file = "inheritance.scatter.plot.pdf",width = 7,height = 3)
ggplot(data = inherit.data%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. melanogaster* dominant",
                                              "*D. simulans* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wounded','Response')))%>%arrange(inher,treatment),                  
       mapping = aes(x = `logFC_f1.mel`,y = `logFC_f1.sim`,colour = inher))+
  scale_colour_manual(values = plot.colours)+
  geom_point(size = 1,shape=1)+
  geom_vline(xintercept = 0,colour='light grey',linetype=3)+geom_hline(yintercept = 0,colour='light grey',linetype=3)+
  # xlim(-12,19)+ylim(-12,19)+
  guides(colour = guide_legend(override.aes =  list(size = 3.5,shape = 15)))+
  facet_grid(cols = vars(treatment))+
  labs(x = "Log<sub>2</sub>( Hybrid / *D. melanogaster* )",
       y = "Log<sub>2</sub>( Hybrid / *D. simulans* )")+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())
dev.off()

pdf(file = "inheritance.bar.plot.pdf",width = 7,height = 1.5)
ggplot(data = inherit.data%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. melanogaster* dominant",
                                              "*D. simulans* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wounded','Response')))%>%arrange(inher,treatment))+
  scale_fill_manual(values = plot.colours)+
  geom_bar(mapping = aes(y=inher,fill = inher))+
  facet_grid(cols = vars(treatment))+
  labs(x = 'Number of genes',y='')+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.text.y = element_markdown())
dev.off()











pdf(file = "inheritance.scatter.plot_immue.pdf",width = 7,height = 3)
ggplot(data = inherit.data%>%filter(dmel_gene%in%genes.sig.in.any)%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. melanogaster* dominant",
                                              "*D. simulans* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wounded','Response')))%>%arrange(inher,treatment),                  
       mapping = aes(x = `logFC_f1.mel`,y = `logFC_f1.sim`,colour = inher))+
  scale_colour_manual(values = plot.colours)+
  geom_point(size = 1,shape=1)+
  geom_vline(xintercept = 0,colour='light grey',linetype=3)+geom_hline(yintercept = 0,colour='light grey',linetype=3)+
  # xlim(-12,19)+ylim(-12,19)+
  guides(colour = guide_legend(override.aes =  list(size = 3.5,shape = 15)))+
  facet_grid(cols = vars(treatment))+
  labs(x = "Log<sub>2</sub>( Hybrid / *D. melanogaster* )",
       y = "Log<sub>2</sub>( Hybrid / *D. simulans* )")+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())
dev.off()

pdf(file = "inheritance.bar.plot_immue.pdf",width = 7,height = 1.5)
ggplot(data = inherit.data%>%filter(dmel_gene%in%genes.sig.in.any)%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. melanogaster* dominant",
                                              "*D. simulans* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wounded','Response')))%>%arrange(inher,treatment))+
  scale_fill_manual(values = plot.colours)+
  geom_bar(mapping = aes(y=inher,fill = inher))+
  facet_grid(cols = vars(treatment))+
  labs(x = 'Number of genes',y='')+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.text.y = element_markdown())
dev.off()








############################
#### logFC1
############################

pdf(file = "inheritance.scatter.plot_immue.logfc1.pdf",width = 5,height = 3)
ggplot(data = inherit.data%>%filter(dmel_gene%in%genes.sig.in.any.logfc1,treatment!="Response")%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. melanogaster* dominant",
                                              "*D. simulans* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wounded','Response')))%>%arrange(inher,treatment),                  
       mapping = aes(x = `logFC_f1.mel`,y = `logFC_f1.sim`,colour = inher))+
  scale_colour_manual(values = plot.colours)+
  geom_point(size = 1,shape=1)+
  geom_vline(xintercept = 0,colour='light grey',linetype=3)+geom_hline(yintercept = 0,colour='light grey',linetype=3)+
  # xlim(-12,19)+ylim(-12,19)+
  guides(colour = guide_legend(override.aes =  list(size = 3.5,shape = 15)))+
  facet_grid(cols = vars(treatment))+
  labs(x = "Log<sub>2</sub>( Hybrid / *D. melanogaster* )",
       y = "Log<sub>2</sub>( Hybrid / *D. simulans* )")+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())
dev.off()

pdf(file = "inheritance.scatter.plot_immue.logfc1_scratch.pdf",width = 7,height = 3)
ggplot(data = inherit.data%>%filter(dmel_gene%in%genes.sig.in.any.logfc1)%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. melanogaster* dominant",
                                              "*D. simulans* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wounded','Response')))%>%arrange(inher,treatment),                  
       mapping = aes(x = `logFC_f1.mel`,y = `logFC_f1.sim`,colour = inher))+
  scale_colour_manual(values = plot.colours)+
  geom_point(size = 1,shape=1)+
  geom_text_repel(aes(label = if_else(condition = abs(logFC_f1.mel)>=5|abs(logFC_f1.sim)>=5|inher=="Additive",true = dmel_gene,false = "" )))+
  geom_vline(xintercept = 0,colour='light grey',linetype=3)+geom_hline(yintercept = 0,colour='light grey',linetype=3)+
  # xlim(-12,19)+ylim(-12,19)+
  guides(colour = guide_legend(override.aes =  list(size = 3.5,shape = 15)))+
  facet_grid(cols = vars(treatment))+
  labs(x = "Log<sub>2</sub>( Hybrid / *D. melanogaster* )",
       y = "Log<sub>2</sub>( Hybrid / *D. simulans* )")+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())
dev.off()

pdf(file = "inheritance.bar.plot_immue.logfc1.pdf",width = 7.7,height = 2)
ggplot(data = inherit.data%>%filter(dmel_gene%in%genes.sig.in.any.logfc1,treatment!="Response")%>%
         mutate(inher=factor(inher,levels = c("Conserved",
                                              "*D. melanogaster* dominant",
                                              "*D. simulans* dominant",
                                              "Additive",
                                              "Overdominant",
                                              "Underdominant")),
                treatment=factor(treatment,levels = c('Control','Wounded','Response')))%>%arrange(inher,treatment))+
  scale_fill_manual(values = plot.colours)+
  geom_bar(mapping = aes(y=inher,fill = inher))+
  facet_grid(cols = vars(treatment))+
  labs(x = 'Number of genes',y='')+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position = 'none',
        axis.text.y = element_markdown(),text = element_text(size = 15))
dev.off()

a<-(inherit.data%>%filter(dmel_gene%in%genes.sig.in.any.logfc1,inher=="Overdominant",treatment=="Control"))$dmel_gene
b<-(inherit.data%>%filter(dmel_gene%in%genes.sig.in.any.logfc1,inher=="Overdominant",treatment=="Wounded"))$dmel_gene
c<-(inherit.data%>%filter(dmel_gene%in%genes.sig.in.any.logfc1,inher=="Overdominant",treatment=="Response"))$dmel_gene
intersect(a,b)
write.table(x = a,file = "overdom.inf1_control.txt",quote = F,row.names = F,col.names = F)
write.table(x = b,file = "overdom.inf1_treatment.txt",quote = F,row.names = F,col.names = F)


# load("ASDE.RData")
# 
# plot.data.immune.logfc1.hist <- merge(plot.data.immune.logfc1.hist,inherit.data%>%select(dmel_gene,treatment,inher),
#                                       by = c("dmel_gene","treatment"))
# pdf('.pdf',width = 7,height = 7)
# ggplot( data=plot.data.immune.logfc1.hist, aes(x=cat2, y=logFC)) +
#   geom_violin() +
#   scale_fill_manual(values=list("*Cis*" = "#000000","*Trans*" = "#C7050F","Control"="white","Treated"="gray40")) +
#   facet_grid(cols = vars(treatment),scales = "free")+
#   labs(x='',
#        y='Magnitude of expression divergence (log<sub>2</sub>fold change)')+
#   theme_bw()+
#   guides(fill = guide_legend(direction = 'horizontal'))+
#   theme(panel.grid = element_blank(),
#         panel.spacing = unit(0,'lines'),
#         legend.title = element_blank(),
#         legend.text = element_markdown(size=20),
#         legend.position="none",
#         axis.title.x = element_markdown(size=20),
#         axis.title.y = element_markdown(size=20),
#         axis.text = element_text(size=20),
#         strip.background = element_rect(fill=F,colour='black'),
#         strip.text.x = element_markdown(size=20))
# dev.off()
# 
# 
