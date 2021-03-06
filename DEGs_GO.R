#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      
qvalueFilter=0.926      

setwd("LOCATION")     
rt=read.table("id.txt",sep="\t",header=T,check.names=F)     
rt=rt[is.na(rt[,"entrezID"])==F,]                           
gene=rt$entrezID
geneFC=rt$logFC
#geneFC=2^geneFC
names(geneFC)=gene
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}
pdf(file="barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		

pdf(file="bubble.pdf",width = 9,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()
		

pdf(file="circos.pdf",width = 9,height = 6.5)
cnet=cnetplot(kk, foldChange=geneFC, showCategory = 5, circular = TRUE, colorEdge = TRUE)
print(cnet)
dev.off()
