

#install.packages("vioplot")

library(vioplot)                                                   
setwd("LOCATION")         
pcr=25                                                            
npcr=76                                                            

rt=read.table("CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F)  
pdf("vioplot.pdf",height=8,width=15)              
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
	  pcrData=rt[1:pcr,i]
    npcrData=rt[(pcr+1):(pcr+npcr),i]
	  vioplot(pcrData,at=3*(i-1),lty=1,add = T,col = 'green')
	  vioplot(npcrData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(pcrData,npcrData)
	  p=round(wilcoxTest$p.value,3)
	  mx=max(c(pcrData,npcrData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

