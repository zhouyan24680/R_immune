

setwd("LOCATION")
input="CIBERSORT.filter.txt"
outpdf="barplot.pdf"

data <- read.table(input,header=T,sep="\t",check.names=F,row.names=1)
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)

pdf(outpdf,height=6,width=10)
par(las=1,mar=c(8,5,4,13))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1,cex=1);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.97,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=0.8)
dev.off()
