

#install.packages("pheatmap")

setwd("LOCATION")    
rt=read.table("CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=t(rt)

library(pheatmap)
Type=c(rep("con",25),rep("treat1",76))    #修改对照和处理组样品数目
names(Type)=colnames(rt)
Type=as.data.frame(Type)

pdf("heatmap.pdf",height=8,width=20)
pheatmap(rt, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         fontsize = 5,
         fontsize_row=13,
         fontsize_col=10)
dev.off()
