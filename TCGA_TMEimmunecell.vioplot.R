#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("vioplot")


library(limma)
library(vioplot)

immuneFile="CIBERSORT-Results.txt"       
expFile="symbol.txt"                     
gene="ORM1"                              
pFilter=0.05                              
setwd("LOCATION")    


immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])


rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]


lowName=colnames(data)[data[gene,]<=median(data[gene,])]
highName=colnames(data)[data[gene,]>median(data[gene,])]

lowImm=intersect(row.names(immune),lowName)
highImm=intersect(row.names(immune),highName)
rt=rbind(immune[lowImm,],immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)


outTab=data.frame()
pdf("vioplot.pdf",height=8,width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")


for(i in 1:ncol(rt)){
	  if(sd(rt[1:lowNum,i])==0){
	    rt[1,i]=0.001
	  }
	  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
	    rt[(lowNum+1),i]=0.001
	  }
	  lowData=rt[1:lowNum,i]
	  highData=rt[(lowNum+1):(lowNum+highNum),i]
	  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = 'green')
	  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(lowData,highData)
	  p=wilcoxTest$p.value
	  if(p<pFilter){
	      cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
		  outTab=rbind(outTab,cellPvalue)
	  }
	  mx=max(c(lowData,highData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()


write.table(outTab,file="diff.result.txt",sep="\t",row.names=F,quote=F)

