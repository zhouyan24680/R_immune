install.packages("ggplot2")

source("http://bioconductor.org/biocLite.R")
biocLite("limma")

source("http://bioconductor.org/biocLite.R")
biocLite("impute")

source("http://bioconductor.org/biocLite.R")
biocLite("pheatmap")

setwd("LOCATION") 

library(limma)
library("impute")


geo_data<-read.table("data.txt",sep="\t",header=T)
geo_data<-as.matrix(geo_data)
rownames(geo_data)=geo_data[,1]
geo_exp<-geo_data[,2:ncol(geo_data)]
dimnames<-list(rownames(geo_exp),colnames(geo_exp))
geo_exp<-matrix(as.numeric(as.matrix(geo_exp)),nrow=nrow(geo_exp),dimnames=dimnames)


mat=impute.knn(geo_exp)
geo_data=mat$data
geo_data=avereps(geo_data)     

pdf(file="raw_box.pdf")
boxplot(geo_data,col = "green",xaxt = "n",outline = F)
dev.off()
geo_data=normalizeBetweenArrays(as.matrix(geo_data))
pdf(file="normal_box.pdf")
boxplot(geo_data,col = "red",xaxt = "n",outline = F)
dev.off()

#geo_data=log2(geo_data)            
class <- c(rep("normal",27),rep("treatment",88))   
design <- model.matrix(~0+factor(class))
colnames(design) <- c("normal","treatment"
fit <- lmFit(geo_data,design)
cont.matrix<-makeContrasts(treatment-normal,levels=design
fit1 <- contrasts.fit(fit, cont.matrix
fit1 <- eBayes(fit1)

allgene<-topTable(fit1,adjust='fdr',number=100000)
write.table(allgene,"allgene.xls",sep="\t",quote=F)
normaldata<-allgene[order(allgene$logFC),]
normaldata1<-rbind(Gene=colnames(normaldata),normaldata)
write.table(normaldata,"normaldata.txt",sep="\t",quote=F,col.names=F)



diffgene <- allgene[with(allgene, (abs(logFC)>=1 & adj.P.Val < 0.926 )), ]
write.table(diffgene,"diffgene.xls",sep="\t",quote=F)

Upgene <- allgene[with(allgene, (logFC>=1 & adj.P.Val < 0.926 )), ]
write.table(Upgene,"upgene.xls",sep="\t",quote=F)

Downgene <- allgene[with(allgene, (logFC<=(-1) & adj.P.Val < 0.926 )), ]
write.table(Downgene,"down.xls",sep="\t",quote=F)

diffexp=geo_data[rownames(diffgene),]
diffexp1=rbind(id=colnames(diffexp),diffexp)
write.table(diffexp1,"diffexp.txt",sep="\t",quote=F,col.names=F)


inputfile_vol="vol.txt"
mydata<-read.table(inputfile_vol,header=T,row.names=1,check.names=F)
down <- mydata[mydata$logFC <= -1 & mydata$adj.P.Val<0.23,] 
up <- mydata[mydata$logFC >= 1 & mydata$adj.P.Val<0.23,]
no <- mydata[(mydata$logFC > -1 & mydata$logFC <1),]
down<- transform(down,adj.P.Val=-log10(down$adj.P.Val))
up<- transform(up,adj.P.Val=-log10(up$adj.P.Val))
no<- transform(no,adj.P.Val=-log10(no$adj.P.Val))
pdf("vol1.pdf")
xm=max(abs(mydata$logFC))
ym=max(-log10(mydata$adj.P.Val))
plot(no$logFC,no$adj.P.Val,xlim = c(-xm,xm),ylim=c(0,ym),col="gray",pch=16,cex=0.9,main = "Volcano",xlab = "logFC",ylab="-log10(adj.P.Val)")
points(up$logFC,up$adj.P.Val,col="red",pch=16,cex=0.9)
points(down$logFC,down$adj.P.Val,col="blue",pch=16,cex=0.9)
abline(v=0,lwd=2)
dev.off()


inputheatmap="diffexp.txt"
library(pheatmap)
data<-read.table(inputheatmap,sep="\t",header=T,row.names=1,check.names=F)
data<-data[1:nrow(data),] 
data<-data[1:20,]    
pdf("heatmap.pdf")
pheatmap(data,display_numbers = F,fontsize_row=7,fontsize_col=5,cluster_cols = T,cluster_rows = T,color = colorRampPalette(c("green", "black", "red"))(50))
dev.off()




