
data.misexp.f1 <- plot.data.immune.logfc1%>%select(dmel_gene,logFC.cis,logFC.trans,treatment)
load("inherit.data.RData")
data.misexp.f1<-merge(x = data.misexp.f1,
                      y = inherit.data%>%select(dmel_gene,treatment,inher),
                      by = c("dmel_gene","treatment"))%>%mutate(direction=logFC.cis*logFC.trans)
data.misexp.f1$direction[data.misexp.f1$direction<0]<-"oppo"
data.misexp.f1$direction[data.misexp.f1$direction!="oppo"]<-"same"
data.misexp.f1$inher<-gsub(pattern = "Overdominant",replacement = "1",x = data.misexp.f1$inher)
data.misexp.f1$inher<-gsub(pattern = "Underdominant",replacement = "0",x = data.misexp.f1$inher)
data.misexp.f1$inher[data.misexp.f1$inher!="1"]<-"0"
data.misexp.f1<-data.misexp.f1%>%mutate(direction=factor(direction,levels = c("same","oppo")))%>%distinct()

# misexp.fit <- glm(formula = direction~inher,family = "binomial",
#                   data = data.misexp.f1%>%filter(treatment=="Control"))
# 
# summary(misexp.fit)



data.misexp.f1.count.data.control <- data.frame(oppo=c(length((data.misexp.f1%>%filter(direction=="oppo",inher=="1",treatment=="Control"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="oppo",inher=="0",treatment=="Control"))$dmel_gene)),
                                                same=c(length((data.misexp.f1%>%filter(direction=="same",inher=="1",treatment=="Control"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="same",inher=="0",treatment=="Control"))$dmel_gene)))%>%
  `row.names<-`(.,c("Misexpression","Non-misexpression"))
data.misexp.f1.count.data.control <- t(data.misexp.f1.count.data.control)

# text.chi.control <- chisq.test(data.misexp.f1.count.data.control)
# text.chi.control
text.fish.control <- fisher.test(data.misexp.f1.count.data.control)
text.fish.control


data.misexp.f1.count.data.wounded <- data.frame(oppo=c(length((data.misexp.f1%>%filter(direction=="oppo",inher=="1",treatment=="Wounded"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="oppo",inher=="0",treatment=="Wounded"))$dmel_gene)),
                                                same=c(length((data.misexp.f1%>%filter(direction=="same",inher=="1",treatment=="Wounded"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="same",inher=="0",treatment=="Wounded"))$dmel_gene)))%>%
  `row.names<-`(.,c("Misexpression","Non-misexpression"))
data.misexp.f1.count.data.wounded <- t(data.misexp.f1.count.data.wounded)

# text.chi.wounded <- chisq.test(data.misexp.f1.count.data.wounded)
# text.chi.wounded
text.fish.wounded <- fisher.test(data.misexp.f1.count.data.wounded)
text.fish.wounded






















data.misexp.f1 <- plot.data.immune.logfc1%>%select(dmel_gene,logFC.cis,logFC.trans,treatment)%>%filter(dmel_gene%in%plot.immune.genes$dmel_gene)%>%distinct()

data.misexp.f1<-merge(x = data.misexp.f1,
                      y = inherit.data%>%select(dmel_gene,treatment,inher),
                      by = c("dmel_gene","treatment"))%>%mutate(direction=logFC.cis*logFC.trans)
data.misexp.f1$direction[data.misexp.f1$direction<0]<-"oppo"
data.misexp.f1$direction[data.misexp.f1$direction!="oppo"]<-"same"
data.misexp.f1$inher<-gsub(pattern = "Overdominant",replacement = "1",x = data.misexp.f1$inher)
data.misexp.f1$inher<-gsub(pattern = "Underdominant",replacement = "0",x = data.misexp.f1$inher)
data.misexp.f1$inher[data.misexp.f1$inher!="1"]<-"0"
data.misexp.f1<-data.misexp.f1%>%mutate(direction=factor(direction,levels = c("same","oppo")))%>%distinct()

misexp.fit <- glm(formula = direction~inher,family = "binomial",
                  data = data.misexp.f1%>%filter(treatment=="Control"))

summary(misexp.fit)



data.misexp.f1.count.data.control <- data.frame(oppo=c(length((data.misexp.f1%>%filter(direction=="oppo",inher=="1",treatment=="Control"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="oppo",inher=="0",treatment=="Control"))$dmel_gene)),
                                                same=c(length((data.misexp.f1%>%filter(direction=="same",inher=="1",treatment=="Control"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="same",inher=="0",treatment=="Control"))$dmel_gene)))%>%
  `row.names<-`(.,c("Misexpression","Non-misexpression"))

text.fish.control <- fisher.test(data.misexp.f1.count.data.control)
text.fish.control


data.misexp.f1.count.data.wounded <- data.frame(oppo=c(length((data.misexp.f1%>%filter(direction=="oppo",inher=="1",treatment=="Wounded"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="oppo",inher=="0",treatment=="Wounded"))$dmel_gene)),
                                                same=c(length((data.misexp.f1%>%filter(direction=="same",inher=="1",treatment=="Wounded"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="same",inher=="0",treatment=="Wounded"))$dmel_gene)))%>%
  `row.names<-`(.,c("Misexpression","Non-misexpression"))

text.fish.wounded <- fisher.test(data.misexp.f1.count.data.wounded)
text.fish.wounded














data.misexp.f1 <- plot.data.immune.logfc1%>%select(dmel_gene,logFC.cis,logFC.trans,treatment)%>%
  filter(dmel_gene%in%c("AttA","AttB","AttC","AttD","CecA1",""))%>%distinct()

data.misexp.f1<-merge(x = data.misexp.f1,
                      y = inherit.data%>%select(dmel_gene,treatment,inher),
                      by = c("dmel_gene","treatment"))%>%mutate(direction=logFC.cis*logFC.trans)
data.misexp.f1$direction[data.misexp.f1$direction<0]<-"oppo"
data.misexp.f1$direction[data.misexp.f1$direction!="oppo"]<-"same"
data.misexp.f1$inher<-gsub(pattern = "Overdominant",replacement = "1",x = data.misexp.f1$inher)
data.misexp.f1$inher<-gsub(pattern = "Underdominant",replacement = "0",x = data.misexp.f1$inher)
data.misexp.f1$inher[data.misexp.f1$inher!="1"]<-"0"
data.misexp.f1<-data.misexp.f1%>%mutate(direction=factor(direction,levels = c("same","oppo")))%>%distinct()

misexp.fit <- glm(formula = direction~inher,family = "binomial",
                  data = data.misexp.f1%>%filter(treatment=="Control"))

summary(misexp.fit)



data.misexp.f1.count.data.control <- data.frame(oppo=c(length((data.misexp.f1%>%filter(direction=="oppo",inher=="1",treatment=="Control"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="oppo",inher=="0",treatment=="Control"))$dmel_gene)),
                                                same=c(length((data.misexp.f1%>%filter(direction=="same",inher=="1",treatment=="Control"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="same",inher=="0",treatment=="Control"))$dmel_gene)))%>%
  `row.names<-`(.,c("Misexpression","Non-misexpression"))

text.fish.control <- fisher.test(data.misexp.f1.count.data.control)
text.fish.control


data.misexp.f1.count.data.wounded <- data.frame(oppo=c(length((data.misexp.f1%>%filter(direction=="oppo",inher=="1",treatment=="Wounded"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="oppo",inher=="0",treatment=="Wounded"))$dmel_gene)),
                                                same=c(length((data.misexp.f1%>%filter(direction=="same",inher=="1",treatment=="Wounded"))$dmel_gene),
                                                       length((data.misexp.f1%>%filter(direction=="same",inher=="0",treatment=="Wounded"))$dmel_gene)))%>%
  `row.names<-`(.,c("Misexpression","Non-misexpression"))

text.fish.wounded <- fisher.test(data.misexp.f1.count.data.wounded)
text.fish.wounded
















