

#install.packages("corrplot")

setwd("LOCATION")     
rt=read.table("CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F)

library(corrplot)
pdf("corHeatmap.pdf",height=13,width=13)             
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         )
dev.off()

