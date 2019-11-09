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



### Figure 3 ####



