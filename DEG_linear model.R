##extract gene expression from GSE
library(data.table)
library(ggplot2)
path<-setwd("/scratch/users/s/h/shifang/HEM/RNA_seq/seq_list")
fileNames = list.files(path=path,pattern="*.tsv", full.names = TRUE)
tbl3<-read.csv("meta.txt",header=T)
alphabeita <- matrix(nrow =198284 ,
                     ncol = 38)
alphabeita<-data.frame(alphabeita)
for(i in 1:length(fileNames)){
  gumlm<-fread(fileNames[i],header=T)
  gumlm1<-gumlm[,c(1,3)]
  data<-gumlm1[match(tbl3$gene, gumlm1$gene),]
  alphabeita[,i]<-data[, 2]
}
write.csv(alphabeita,"HEM_all_gene.csv",quote = FALSE)

##different Gene expression using linear model correcting the effect of sex and BMI
data1<-read.csv("/scratch/users/s/h/shifang/HEM/RNA_seq/HEM_all_gene.csv",row.names = 1)
data2<-read.csv("/scratch/users/s/h/shifang/HEM/RNA_seq/meta.csv",header=F)
data3<-data.frame(colnames(data1))
colnames(data1)<-data2$V1
data0<-read.csv("/scratch/users/s/h/shifang/HEM/RNA_seq/HO_meta.csv")
head (data1)
data12<-data.frame(t(data1[c("ENST00000185150.8:protein_coding:ERLEC1-001","ENST00000375544.7:protein_coding:ASPN-001"),]))
data12<-data12[data0$ID,]
data12<-cbind(data12,data0)

alphabeita <- matrix(nrow =2,
                     ncol =3)
colnames(alphabeita)<-c("pvalue_sex","pvalue_Group","pvalue_BMI")
rownames(alphabeita)<-c("ERLEC1","ASPN")
for (j in c(1:2)){
  fit<-lm(log2(data12[,j])~data12$Sex+data12$Group+data12$BMI)
  alphabeita[j,1]<- coef(summary(fit))[2,c(4)]
  alphabeita[j,2]<- coef(summary(fit))[3,c(4)]
  alphabeita[j,3]<- coef(summary(fit))[4,c(4)]
}

#alphabeita
#        pvalue_sex pvalue_Group pvalue_BMI
#ERLEC1  0.3619023  0.007120946  0.0922203
#ASPN    0.2122636  0.190137011  0.5865649

data0$Group<-factor(data0$Group,levels=c("Control","HEM"))
p1<-ggplot(data=data0, aes(x = Group, y =gene, fill = Group))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 3), 
                                                                                                   legend.text = element_text(size = 3))+labs(x="",y="")+geom_jitter(shape=16,position=position_jitter(0.2))+theme(legend.position="none")+labs(x="",y="log2(RPM of ERLEC1)")+scale_fill_manual(values = pal <- c("#01B7FF","#D9AB42"))
p2<-ggplot(data=data0, aes(x = Group, y =gene, fill = Group))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 3), 
                                                                                                   legend.text = element_text(size = 3))+labs(x="",y="")+geom_jitter(shape=16,position=position_jitter(0.2))+theme(legend.position="none")+labs(x="",y="log2(RPM of ASPN)")+scale_fill_manual(values = pal <- c("#01B7FF","#D9AB42"))
p1+p2