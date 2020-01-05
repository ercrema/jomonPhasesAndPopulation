### Figure 1 ####
library(rcarbon)
# Load MCMC Posterior Samples 
load("../../R_images/posteriorSamples.RData")
phases =c("S0","S1.1","S1.2","S2.1","S2.2",paste0("S",3:8),paste0("Z",1:7),"C1","C234","C56","C78",paste0("C",9:14),paste0("K",1:8),paste0("B",1:6))
# Convert any post-bomb age to postbomb
postTrapezoid$posterior[which(postTrapezoid$posterior>0)]=-1
pdf(file = "./figure1.pdf",width = 7,height = 9)
plot(0,0,xlim=c(18000,2000),ylim=c(0,length(phases)+1),type='n', ylab="",xlab="",axes=F)
axis(1,at=seq(18000,2000,-1000),labels=seq(18000,2000,-1000)/1000,tck=-0.01,padj=-1,cex=0.7)
axis(3,at=seq(18000,2000,-1000),labels=seq(18000,2000,-1000)/1000,tck=-0.01,padj=1,cex=0.7)
axis(1,at=seq(18000,2000,-500),labels=NA,tck=-0.005)
axis(3,at=seq(18000,2000,-500),labels=NA,tck=-0.005)
mtext(1,line=2,text = 'k cal BP')
mtext(3,line=2,text = 'k cal BP')

n=length(18000:2000)
#axis(2,at=1:length(phases)+0.5,labels=phases,las=2,cex.axis=0.8,tck=-0.01,hadj=0.55)
abline(v=seq(18000,2000,-1000),lty=3,lwd=1.5,col="lightgrey")
abline(v=seq(18000,2000,-500),lty=3,lwd=0.5,col="lightgrey")

#abline(h=1:length(phases),lty=4,col="lightgrey")

for (k in 1:length(phases))
{
  
  dens1 = density(BCADtoBP(postTrapezoid$posterior[k,1,]))
  dens1.x = dens1$x
  dens1.y = reScale(dens1$y)+k
  
  dens2 = density(BCADtoBP(postTrapezoid$posterior[k,2,]))
  dens2.x = dens2$x
  dens2.y = reScale(dens2$y)+k
  
  dens3 = density(BCADtoBP(postTrapezoid$posterior[k,3,]))
  dens3.x = dens3$x
  dens3.y = reScale(dens3$y)+k
  
  dens4 = density(BCADtoBP(postTrapezoid$posterior[k,4,]))
  dens4.x = dens4$x
  dens4.y = reScale(dens4$y)+k
  
  polygon(c(dens1.x,rev(dens1.x)),c(dens1.y,rep(k,length(dens1.y))),col=rgb(0.4,0.76,0.65,0.8),border=NA)
  polygon(c(dens2.x,rev(dens2.x)),c(dens2.y,rep(k,length(dens2.y))),col=rgb(0.91,0.54,0.76,0.8),border=NA)
  polygon(c(dens3.x,rev(dens3.x)),c(dens3.y,rep(k,length(dens3.y))),col=rgb(0.55,0.63,0.79,0.8),border=NA)
  polygon(c(dens4.x,rev(dens4.x)),c(dens4.y,rep(k,length(dens4.y))),col=rgb(0.99,0.55,0.38,0.8),border=NA)
  
  text(x=median(BCADtoBP(postTrapezoid$posterior[k,1,]))+700,y=0.5+k,phases[k],cex=0.8)
}

rect(xleft=18400,xright=10000,ybottom=36,ytop=44,border=NA,col="white")
legend("topleft",legend=c("a","b","c","d"),fill=c(rgb(0.4,0.76,0.65,0.8),rgb(0.91,0.54,0.76,0.8),rgb(0.55,0.63,0.79,0.8),rgb(0.99,0.55,0.38,0.8)),bty="n")
lines(x=c(15000,14500,12200,11000,15000),y=c(38,43,43,38,38))
text(x=15300,y=37.5,"a")
text(x=14750,y=43.5,"b")
text(x=11850,y=43.5,"c")
text(x=10800,y=37.5,"d")
dev.off()

### Figure 2 ####
load("../../R_images/simdatesPithouses.RData")

pdf(file = "./figure2.pdf",width = 5,height = 8)
par(mfrow=c(2,1),mar=c(3,4,1,1))
plot(0,0,xlim=c(8000,3000),ylim=range(densSummary.trap.swkanto[,-1],na.rm=TRUE),type="n",xlab="",ylab="",axes=FALSE)
polygon(x=c(densSummary.trap.swkanto$CalBP,rev(densSummary.trap.swkanto$CalBP)),y=c(densSummary.trap.swkanto$lo,rev(densSummary.trap.swkanto$hi)),col=rgb(0,0,0,0.2),border=NA)
lines(densSummary.trap.swkanto$CalBP,densSummary.trap.swkanto$mean)
axis(side=1,at=seq(8000,3000,-1000),labels=seq(8000,3000,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=seq(8000,3000,-500),labels=NA,tck=-0.02)
axis(side=1,at=seq(8000,3000,-100),labels=NA,tck=-0.01)
mtext(2,line=3,text = 'Composite Kernel Density Estimate')
mtext(1,line=2,text = 'cal BP')
axis(side=2)
box()
legend("topright",legend="a",cex=1,bty='n')
legend("topleft",legend=c("95 Percentile Interval","Mean KDE"),lwd=c(7,1),col=c("lightgrey","black"),bty='n',cex=0.9)
plot(0,0,xlim=c(8000,3000),ylim=range(densSummary.trap.chubukochi[,-1],na.rm=TRUE),type="n",xlab="",ylab="",axes=FALSE)
polygon(x=c(densSummary.trap.chubukochi$CalBP,rev(densSummary.trap.chubukochi$CalBP)),y=c(densSummary.trap.chubukochi$lo,rev(densSummary.trap.chubukochi$hi)),col=rgb(0,0,0,0.2),border=NA)
lines(densSummary.trap.chubukochi$CalBP,densSummary.trap.chubukochi$mean)
axis(side=1,at=seq(8000,2500,-1000),labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=seq(8000,2500,-500),labels=NA,tck=-0.02)
axis(side=1,at=seq(8000,2500,-100),labels=NA,tck=-0.01)
mtext(2,line=3,text = 'Composite Kernel Density Estimate')
mtext(1,line=2,text = 'cal BP')
axis(side=2)
box()
legend("topright",legend="b",cex=1,bty='n')
legend("topleft",legend=c("95 Percentile Interval","Mean KDE"),lwd=c(7,1),col=c("lightgrey","black"),bty='n',cex=0.9)
dev.off()



### Figure 3 ####
load("../../R_images/simdatesPithouses.RData")
load("../../R_images/spdRes.RData")
load("../../R_images/comp.RData")
source("../../R/utilities.R")

pdf(file = "./figure3.pdf",width = 5.5,height = 7)
par(mfrow=c(3,1),mar=c(0,4,1,1)+0.1)
b=barplot(apply(tblocks.trap,1,mean),border=NA,col="darkorange",ylim=c(0,1800),space = 0)
error.bar(b,upper=apply(tblocks.trap,1,quantile,0.975),lower=apply(tblocks.trap,1,quantile,0.025),length=0.02)
abline(v=which(tbs2%in%seq(8000,2500,-500))-1,col="white",lty=3)
#mtext(1,line=2,text = 'cal BP')
axis(side=1,at=which(tbs2%in%seq(8000,2500,-1000))-1,labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-500))-1,labels=NA,tck=-0.02)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-100))-1,labels=NA,tck=-0.01)
mtext(2,line=3,text = 'Number of Pit-dwellings',cex=0.8)
legend("topleft",legend="a",bty = 'n',cex=2)

par(mar=c(5,4,2,1)+0.1)
b2=barplot(spdDataSPD_blocks$sumblock,border=NA,col="royalblue",space = 0)
abline(v=which(tbs2%in%seq(8000,2500,-500))-1,col="white",lty=3)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-1000))-1,labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-500))-1,labels=NA,tck=-0.02)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-100))-1,labels=NA,tck=-0.01)
#mtext(1,line=2,text = 'cal BP')
mtext(2,line=3,text = 'Summed Probability of 14C Dates',cex=0.8)
legend("topleft",legend="b",bty = 'n',cex=2)

par(mar=c(4,4,0,1)+0.1)
plot(b2,apply(tblockRoll10.trap,1,mean,na.rm=TRUE),ylim=c(-1,1),type='n',axes=FALSE,xlab="",ylab="")
polygon(c(b2,rev(b2)),c(apply(tblockRoll10.trap,1,quantile,0.025,na.rm=TRUE),rev(apply(tblockRoll10.trap,1,quantile,0.975,na.rm=TRUE))),border=NA,col='lightgrey')
abline(v=which(tbs2%in%seq(8000,2500,-500))-1,col="white",lty=3)
lines(b2,apply(tblockRoll10.trap,1,mean,na.rm=TRUE),lwd=2)
abline(h=0,lty=2,col='red')
mtext(2,line=3,text = '1000 yrs Rolling Correlation',cex=0.8)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-1000))-1,labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-500))-1,labels=NA,tck=-0.02)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-100))-1,labels=NA,tck=-0.01)
axis(side=2,at=c(-1,-0.5,0,0.5,1))
mtext(1,line=2,text = 'cal BP')
legend("topleft",legend="c",bty = 'n',cex=2)
dev.off()


### Figure 4 ####

load("../../R_images/comp.RData")
library(rcarbon)
pdf(file = "./figure4.pdf",width = 6,height = 5)

plot(res.compare,type='roc',xlim=c(7500,3500),ylim=c(-1,2),drawaxes = FALSE)
axis(side=1,at=seq(8000,2500,-1000),labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=seq(8000,2500,-500),labels=NA,tck=-0.02)
axis(side=1,at=seq(8000,2500,-100),labels=NA,tck=-0.01)
axis(side=2)
mtext(2,line=3,text = 'Annual Growth Rate (%)')
mtext(1,line=2,text = 'cal BP')
box()
legend("topleft",legend=c("Expectation from Pithouse Data","Positive Deviation","Negative Deviation"),bg="white",cex=0.8,fill=c("lightgrey",rgb(0.7,0,0.0,0.3),rgb(0,0,0.7,0.3)))
dev.off()
