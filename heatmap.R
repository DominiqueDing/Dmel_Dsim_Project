


data.heatmap <- merge.data.frame(
  result.within.DE%>%filter(species=="*D. melanogaster*",dmel_gene%in%genes.sig.in.any)%>%select(dmel_gene,logFC)%>%`colnames<-`(.,c("dmel_gene","*D. melanogaster*")),
  result.within.DE%>%filter(species=="*D. simulans*",dmel_gene%in%genes.sig.in.any)%>%select(dmel_gene,logFC)%>%`colnames<-`(.,c("dmel_gene","*D. simulans*"))
)
data.heatmap <- merge.data.frame(
  data.heatmap,
  result.within.DE%>%filter(species=="Hybrid",dmel_gene%in%genes.sig.in.any)%>%select(dmel_gene,logFC)%>%`colnames<-`(.,c("dmel_gene","Hybrid"))
)


annotdf.reac <- merge.data.frame(
  data.heatmap%>%select(dmel_gene),
  inherit.data%>%filter(treatment=="Response")%>%select(dmel_gene,inher),by="dmel_gene",all.x = T
)%>%filter(dmel_gene%in%genes.sig.in.any.logfc1)
rownames(annotdf.reac)<-annotdf.reac$dmel_gene
annotdf.reac<-annotdf.reac%>%select(-dmel_gene)

rownames(data.heatmap)<-data.heatmap$dmel_gene
data.heatmap<-data.heatmap%>%select(-dmel_gene)
data.heatmap<-t(as.matrix(data.heatmap))

library("pheatmap")
pheatmap(data.heatmap,cutree_rows = 3,cutree_cols = 1,
         border_color = NA,drop_levels = T,show_colnames = F,
         # treeheight_row = 0,treeheight_col = 0,,height = 2.5,width = 8,
         height = 3,width = 8,
         labels_row = c(expression(italic(D.~melanogaster)),expression(italic(D.~simulans)),"Hybrid"),
         filename = "heatmap.pdf")

heat <- pheatmap(data.heatmap,cutree_rows = 3,cutree_cols = 4,
                 border_color = NA,drop_levels = T,show_colnames = F,
                 # treeheight_row = 0,treeheight_col = 0,,height = 2.5,width = 8,
                 height = 3,width = 8,
                 labels_row = c(expression(italic(D.~melanogaster)),expression(italic(D.~simulans)),"Hybrid"))
heat$tree_col
col <- cutree(heat$tree_col,k=3)
col[col[]!=1]






data.heatmap.logfc1 <- data.heatmap[,colnames(data.heatmap)%in%genes.sig.in.any.logfc1]
data.heatmap.logfc1 <- data.heatmap.logfc1[,c(rownames(annotdf.reac%>%filter(inher=="*D. melanogaster* dominant")),
                                              rownames(annotdf.reac%>%filter(inher=="*D. simulans* dominant")),
                                              rownames(annotdf.reac%>%filter(inher=="Additive")),
                                              rownames(annotdf.reac%>%filter(inher=="Overdominant")),
                                              rownames(annotdf.reac%>%filter(inher=="Underdominant")),
                                              rownames(annotdf.reac%>%filter(inher=="Conserved")))]
pheatmap(data.heatmap.logfc1,cutree_rows = 3,cutree_cols = 1,
         border_color = NA,drop_levels = T,show_colnames = F,
         # treeheight_row = 0,treeheight_col = 0,,height = 2.5,width = 8,
         height = 3,width = 8,
         labels_row = c(expression(italic(D.~melanogaster)),expression(italic(D.~simulans)),"Hybrid"),
         filename = "heatmap.logfc1.pdf")

pheatmap(data.heatmap.logfc1,cutree_rows = 3,cutree_cols = 1,
         border_color = NA,drop_levels = T,show_colnames = T,
         # treeheight_row = 0,treeheight_col = 0,,height = 2.5,width = 8,
         height = 3,width = 25,
         labels_row = c(expression(italic(D.~melanogaster)),expression(italic(D.~simulans)),"Hybrid"),
         filename = "heatmap.logfc1.scratch.pdf")





pheatmap(data.heatmap.logfc1,cutree_rows = 3,cutree_cols = 1,
         border_color = NA,drop_levels = T,show_colnames = F,
         height = 3,width = 8,
         labels_row = c(expression(italic(D.~melanogaster)),expression(italic(D.~simulans)),"Hybrid"),
         filename = "heatmap.logfc1_response.pdf",
         cluster_cols = FALSE,
         annotation_col = annotdf.reac,
         annotation_colors = list(
           inher=c('*D. melanogaster* dominant' = "#F8766D",
           '*D. simulans* dominant' = "#00BFC4",
           "Additive" = "#FFC300",
           "Overdominant" = "#089911",
           "Underdominant" = "#9B5DE5",
           "Conserved" = "#778DA9"))
)


