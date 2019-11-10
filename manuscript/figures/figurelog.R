### Figure 1 ####

# Load MCMC Posterior Samples 
load("../../R_images/posteriorSamples.RData")

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
load("../../R_images/pithouseData.RData")
load("../../R_images/spdRes.RData")

nsim = ncol(simTrapezoid)
#tbs = seq(14950,650,-100)
tbs = seq(7950,2550,-100)
tbs2 = seq(8000,2500,-100)

tblocks = tblocksCal =matrix(NA,nrow=length(tbs),ncol=nsim)

for (s in 1:nsim)
{
  #tblocks[,s]=as.numeric(rev(table(cut(simTrapezoid[,s],breaks=seq(15000,600,-100)))))
  tblocks[,s]=as.numeric(rev(table(cut(simTrapezoid[,s],breaks=tbs2))))
  tblocksCal[,s]=as.numeric(rev(table(cut(t(westKantoSPD_sampled$sdates)[,s],breaks=tbs2))))
}

tblockRoll10 = rollCor(tblocks,tblocksCal,rollsize = 10)
overallCorr = numeric(length=nsim)
for (s in 1:nsim)
{
  overallCorr[s] = cor(tblocks[,s],tblocksCal[,s])
}
mean(overallCorr)
quantile(overallCorr,c(0.025,0.975))




pdf(file = "./figure2.pdf",width = 6,height = 7)
par(mfrow=c(3,1),mar=c(0,4,0,1)+0.1)
b=barplot(apply(tblocks,1,mean),border=NA,col="darkorange",ylim=range(tblocks),space = 0)
error.bar(b,upper=apply(tblocks,1,quantile,0.975),lower=apply(tblocks,1,quantile,0.025),length=0.02)
abline(v=which(tbs2%in%seq(8000,2500,-500))-1,col="white",lty=3)
#mtext(1,line=2,text = 'cal BP')
axis(side=1,at=which(tbs2%in%seq(8000,2500,-1000))-1,labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-500))-1,labels=NA,tck=-0.02)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-100))-1,labels=NA,tck=-0.01)
mtext(2,line=3,text = 'Number of Pit-dwellings',cex=0.7)
legend("topleft",legend="a",bty = 'n',cex=2)

par(mar=c(5,4,2,1)+0.1)
b2=barplot(westKantoSPD2$sumblock,border=NA,col="royalblue",space = 0)
abline(v=which(tbs2%in%seq(8000,2500,-500))-1,col="white",lty=3)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-1000))-1,labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-500))-1,labels=NA,tck=-0.02)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-100))-1,labels=NA,tck=-0.01)
#mtext(1,line=2,text = 'cal BP')
mtext(2,line=3,text = 'Summed Probability of 14C Dates',cex=0.7)
legend("topleft",legend="b",bty = 'n',cex=2)

par(mar=c(4,4,0,1)+0.1)
plot(b2,apply(tblockRoll10,1,mean,na.rm=TRUE),ylim=c(-1,1),type='n',axes=FALSE,xlab="",ylab="")
polygon(c(b2,rev(b2)),c(apply(tblockRoll10,1,quantile,0.025,na.rm=TRUE),rev(apply(tblockRoll10,1,quantile,0.975,na.rm=TRUE))),border=NA,col='lightgrey')
abline(v=which(tbs2%in%seq(8000,2500,-500))-1,col="white",lty=3)
lines(b2,apply(tblockRoll10,1,mean,na.rm=TRUE),lwd=2)
abline(h=0,lty=2,col='red')
mtext(2,line=3,text = '1000 yrs Rolling Correlation',cex=0.7)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-1000))-1,labels=seq(8000,2500,-1000),tck=-0.03,padj=-0.6)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-500))-1,labels=NA,tck=-0.02)
axis(side=1,at=which(tbs2%in%seq(8000,2500,-100))-1,labels=NA,tck=-0.01)
axis(side=2,at=c(-1,-0.5,0,0.5,1))
mtext(1,line=2,text = 'cal BP')
legend("topleft",legend="c",bty = 'n',cex=2)
dev.off()


