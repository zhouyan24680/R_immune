#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot"
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05     
qvalueFilter=0.926    

setwd("LOCATION")     
rt=read.table("id.txt",sep="\t",header=T,check.names=F) 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
colnames(rt)[1]="Gene"
gene=rt$entrezID
geneFC=rt$logFC
#geneFC=2^gen
names(geneFC)=gene
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)


showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}


pdf(file="barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()


pdf(file="bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

pdf(file="circos.pdf",width = 11,height = 7)
kkx=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx, foldChange=geneFC,showCategory = 5, circular = TRUE, colorEdge = TRUE,node_label="all")
dev.off()
