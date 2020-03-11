# Load R packages ####
library(magrittr)
library(dplyr)
library(oxcAAR)
library(rcarbon)
library(TTR)
library(trapezoid)
library(magrittr)
library(readr)

# Source R functions ####
source("./R/utilities.R")
source("./R/oxcalReadjs.R")
source("./R/oxcalScriptCreator.R")
source("./R/outlierAnalysis.R")
source("./R/mcsim.R")

# General Settings ####
nsim = 5000

## Read and Prepare Data for Ceramic Phase Models####
c14data = read.csv("./data/c14dates.csv")

## Outlier Analysis ###
grps = unique(c14data$CombineGroup)
grps = grps[which(!is.na(grps))]
c14data$outlier=FALSE

for (i in 1:length(grps))
{
  print(paste0(i," out of ",length(grps)))
  tmp=subset(c14data,CombineGroup==grps[i])
  tmp2=outlierExcluder(c14=tmp$CRA,errors=tmp$Error,id=tmp$LabCode)
  c14data$outlier[match(tmp2$df$id,c14data$LabCode)]=tmp2$df$exclude
}
#remove outliers
c14data=subset(c14data,!outlier)
#create R image file
save(c14data,file="./R_images/c14data.RData")

## Generate Summary Table ###
phases =c("S0","S1.1","S1.2","S2.1","S2.2",paste0("S",3:8),paste0("Z",1:7),"C1","C234","C56","C78",paste0("C",9:14),paste0("K",1:8),paste0("B",1:6))

table1 = data.frame(phases,samples=NA,effsamples=NA,nsites=NA)


for (i in 1:length(phases))
{
  tmp = subset(c14data,PhaseAnalysis==phases[i])
  table1$samples[i] = nrow(tmp)
  table1$effsamples[i] = sum(is.na(tmp$CombineGroup)) + length(unique(tmp$CombineGroup[which(!is.na(tmp$CombineGroup))]))
  table1$nsites[i] = length(unique(tmp$SiteID))
}

# store base table (to be manually edited)
write.csv(table1,file="./manuscript/tables/table1_base.csv",row.names=FALSE)







## Ceramic Phase Modelling (via OxCal) ####

## Generate OxCal Script ##
oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/gaussian.oxcal",mcname="mcmcGaussian",model="gaussian",mcnsim=nsim)
oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/uniform.oxcal",mcname="mcmcUniform",model="uniform",mcnsim=nsim)
oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/trapezoid.oxcal",mcname="mcmcTrapezoid",model="trapezoid",mcnsim=nsim)

## Retrieve Oxcal Output (from Oxcal Online - notice this requires about 70-90 hours of analysis) ##
df = data.frame(id=as.character(c14data$LabCode),grp=c14data$CombineGroup,stringsAsFactors = FALSE)
# Read js files (this takes some time) 
gaussian.agreement = oxcalReadjs(x=df, model='gaussian',path='./oxcal/results/')
uniform.agreement = oxcalReadjs(x=df, model='uniform',path='./oxcal/results/')
trapezoid.agreement = oxcalReadjs(x=df, model='trapezoid',path='./oxcal/results/')

## Generate subsets of dates with agreement indices above 60 ##
c14data.gaussian.rerun = left_join(c14data,gaussian.agreement$df,by=c("LabCode"="id")) %>%
  subset(agreement>60|combine.agreement>60)
c14data.uniform.rerun = left_join(c14data,uniform.agreement$df,by=c("LabCode"="id")) %>%
  subset(agreement>60|combine.agreement>60)
c14data.trapezoid.rerun = left_join(c14data,trapezoid.agreement$df,by=c("LabCode"="id")) %>%
  subset(agreement>60|combine.agreement>60)

## Regenerate OxCal Scripts ##
oxcalScriptGen(id=c14data.gaussian.rerun$LabCode,c14age=c14data.gaussian.rerun$CRA,errors=c14data.gaussian.rerun$Error,group=c14data.gaussian.rerun$CombineGroup,phases=c14data.gaussian.rerun$PhaseAnalysis,fn="./oxcal/oxcalscripts/gaussianR.oxcal",mcname="mcmcGaussianR",model="gaussian",mcnsim=nsim)
oxcalScriptGen(id=c14data.uniform.rerun$LabCode,c14age=c14data.uniform.rerun$CRA,errors=c14data.uniform.rerun$Error,group=c14data.uniform.rerun$CombineGroup,phases=c14data.uniform.rerun$PhaseAnalysis,fn="./oxcal/oxcalscripts/uniformR.oxcal",mcname="mcmcUniformR",model="uniform",mcnsim=nsim)
oxcalScriptGen(id=c14data.trapezoid.rerun$LabCode,c14age=c14data.trapezoid.rerun$CRA,errors=c14data.trapezoid.rerun$Error,group=c14data.trapezoid.rerun$CombineGroup,phases=c14data.trapezoid.rerun$PhaseAnalysis,fn="./oxcal/oxcalscripts/trapezoidR.oxcal",mcname="mcmcTrapezoidR",model="trapezoid",mcnsim=nsim)

## Retrieve Oxcal Output (from Oxcal Online - notice this requires about 70-90 hours of analysis) ##
df.gaussian.rerun = data.frame(id=as.character(c14data.gaussian.rerun$LabCode),grp=c14data.gaussian.rerun$CombineGroup,stringsAsFactors = FALSE)
df.uniform.rerun = data.frame(id=as.character(c14data.uniform.rerun$LabCode),grp=c14data.uniform.rerun$CombineGroup,stringsAsFactors = FALSE)
df.trapezoid.rerun = data.frame(id=as.character(c14data.trapezoid.rerun$LabCode),grp=c14data.trapezoid.rerun$CombineGroup,stringsAsFactors = FALSE)

gaussian.agreementR = oxcalReadjs(x=df.gaussian.rerun, model='gaussianR',path='./oxcal/results/')
uniform.agreementR = oxcalReadjs(x=df.uniform.rerun, model='uniformR',path='./oxcal/results/')
trapezoid.agreementR = oxcalReadjs(x=df.trapezoid.rerun, model='trapezoidR',path='./oxcal/results/')

# Read MCMC samples (excluding first (pass number) and last (empty) column)
gaussian.samples = read.csv("./oxcal/results/mcmcGaussianR.csv")[,-c(1,86)] 
uniform.samples = read.csv("./oxcal/results/mcmcUniformR.csv")[,-c(1,86)]
trapezoid.samples = read.csv("./oxcal/results/mcmcTrapezoidR.csv")[,-c(1,170)]

# Convert into an Array
phases =c("S0","S1.1","S1.2","S2.1","S2.2",paste0("S",3:8),paste0("Z",1:7),"C1","C234","C56","C78",paste0("C",9:14),paste0("K",1:8),paste0("B",1:6))
postGaussian=convertToArray(gaussian.samples,type="gaussian",phases)
postUniform=convertToArray(uniform.samples,type="uniform",phases)
postTrapezoid=convertToArray(trapezoid.samples,type="trapezium",phases)

# Save output in an R image file
save(df.uniform.rerun,df.uniform.rerun,df.trapezoid.rerun,gaussian.agreementR,uniform.agreementR,trapezoid.agreementR,gaussian.agreement,uniform.agreement,trapezoid.agreement,postTrapezoid,postUniform,postGaussian,file="./R_images/posteriorSamples.RData")


## Simulate Pithouse Dates ####
load("./R_images/pithouseData.RData")
simGaussian=mcsim(pithouseData[,-c(3,4)],nsim=nsim,posterior=postGaussian,weights="variance")
simUniform=mcsim(pithouseData[,-c(3,4)],nsim=nsim,posterior=postUniform,weights="variance")
# Convert any post-bomb age to postbomb
postTrapezoid$posterior[which(postTrapezoid$posterior>0)]=-1
simTrapezoid=mcsim(pithouseData[,-c(3,4)],nsim=nsim,posterior=postTrapezoid,weights="variance")

## Group dates in 100-yrs bins between 8000 and 3000 cal BP
tbs = seq(7950,3050,-100)
tbs2 = seq(8000,3000,-100)
tblocks.gauss = tblocks.unif = tblocks.trap =matrix(NA,nrow=length(tbs),ncol=nsim)

for (s in 1:nsim)
{
  tblocks.trap[,s]=as.numeric(rev(table(cut(simTrapezoid[,s],breaks=tbs2))))
  tblocks.gauss[,s]=as.numeric(rev(table(cut(simGaussian[,s],breaks=tbs2))))
  tblocks.unif[,s]=as.numeric(rev(table(cut(simUniform[,s],breaks=tbs2))))
}

## Composite Kernel Density Estimate for Each subregion
swkanto = which(pithouseData$Region == "SWKanto")
chubukochi = which(pithouseData$Region == "ChubuHighlands")

bws.trap = sapply(1:nsim,function(x,y){density(y[,x])$bw},y=simTrapezoid)
bws.gauss = sapply(1:nsim,function(x,y){density(y[,x])$bw},y=simGaussian)
bws.unif = sapply(1:nsim,function(x,y){density(y[,x])$bw},y=simUniform)

densMat.trap.swkanto = densMat.gauss.swkanto= densMat.unif.swkanto = 
densMat.trap.chubukochi = densMat.gauss.chubukochi= densMat.unif.chubukochi = 
densMat.trap = densMat.gauss= densMat.unif = matrix(NA,nrow=length(8000:3000),ncol=nsim)

for (s in 1:nsim)
{
  tmp.trap.swkanto=density(simTrapezoid[swkanto,s],bw=mean(bws.trap))
  tmp.gauss.swkanto=density(simGaussian[swkanto,s],bw=mean(bws.gauss))
  tmp.unif.swkanto=density(simUniform[swkanto,s],bw=mean(bws.unif))
  tmp.trap.chubukochi=density(simTrapezoid[chubukochi,s],bw=mean(bws.trap))
  tmp.gauss.chubukochi=density(simGaussian[chubukochi,s],bw=mean(bws.gauss))
  tmp.unif.chubukochi=density(simUniform[chubukochi,s],bw=mean(bws.unif))
  tmp.trap=density(simTrapezoid[,s],bw=mean(bws.trap))
  tmp.gauss=density(simGaussian[,s],bw=mean(bws.gauss))
  tmp.unif=density(simUniform[,s],bw=mean(bws.unif))
  
  
  densMat.trap.swkanto[,s]=approx(x = tmp.trap.swkanto$x,y=tmp.trap.swkanto$y,xout=8000:3000)$y
  densMat.gauss.swkanto[,s]=approx(x = tmp.gauss.swkanto$x,y=tmp.gauss.swkanto$y,xout=8000:3000)$y
  densMat.unif.swkanto[,s]=approx(x = tmp.unif.swkanto$x,y=tmp.unif.swkanto$y,xout=8000:3000)$y  
  densMat.trap.chubukochi[,s]=approx(x = tmp.trap.chubukochi$x,y=tmp.trap.chubukochi$y,xout=8000:3000)$y
  densMat.gauss.chubukochi[,s]=approx(x = tmp.gauss.chubukochi$x,y=tmp.gauss.chubukochi$y,xout=8000:3000)$y
  densMat.unif.chubukochi[,s]=approx(x = tmp.unif.chubukochi$x,y=tmp.unif.chubukochi$y,xout=8000:3000)$y
  densMat.trap[,s]=approx(x = tmp.trap$x,y=tmp.trap$y,xout=8000:3000)$y
  densMat.gauss[,s]=approx(x = tmp.gauss$x,y=tmp.gauss$y,xout=8000:3000)$y
  densMat.unif[,s]=approx(x = tmp.unif$x,y=tmp.unif$y,xout=8000:3000)$y
}

# Extract only lo, hi, and median values for each CKDE
densSummary.trap.swkanto = data.frame(CalBP=8000:3000,
                                      lo=apply(densMat.trap.swkanto,1,quantile,0.025),
                                      hi = apply(densMat.trap.swkanto,1,quantile,0.975),
                                      mean = apply(densMat.trap.swkanto,1,mean))
densSummary.trap.chubukochi = data.frame(CalBP=8000:3000,
                                         lo=apply(densMat.trap.chubukochi,1,quantile,0.025),
                                         hi = apply(densMat.trap.chubukochi,1,quantile,0.975),
                                         mean = apply(densMat.trap.chubukochi,1,mean))
densSummary.trap = data.frame(CalBP=8000:3000,
                              lo=apply(densMat.trap,1,quantile,0.025),
                              hi = apply(densMat.trap,1,quantile,0.975),
                              mean = apply(densMat.trap,1,mean))

densSummary.gauss.swkanto = data.frame(CalBP=8000:3000,
                                       lo=apply(densMat.gauss.swkanto,1,quantile,0.025),
                                       hi = apply(densMat.gauss.swkanto,1,quantile,0.975),
                                       mean = apply(densMat.gauss.swkanto,1,mean))
densSummary.gauss.chubukochi = data.frame(CalBP=8000:3000,
                                          lo=apply(densMat.gauss.chubukochi,1,quantile,0.025),
                                          hi = apply(densMat.gauss.chubukochi,1,quantile,0.975),
                                          mean = apply(densMat.gauss.chubukochi,1,mean))
densSummary.gauss = data.frame(CalBP=8000:3000,
                               lo=apply(densMat.gauss,1,quantile,0.025),
                               hi = apply(densMat.gauss,1,quantile,0.975),
                               mean = apply(densMat.gauss,1,mean))

densSummary.unif.swkanto = data.frame(CalBP=8000:3000,
                                      lo=apply(densMat.unif.swkanto,1,quantile,0.025),
                                      hi = apply(densMat.unif.swkanto,1,quantile,0.975),
                                      mean = apply(densMat.unif.swkanto,1,mean))
densSummary.unif.chubukochi = data.frame(CalBP=8000:3000,
                                         lo=apply(densMat.unif.chubukochi,1,quantile,0.025),
                                         hi = apply(densMat.unif.chubukochi,1,quantile,0.975),
                                         mean = apply(densMat.unif.chubukochi,1,mean))
densSummary.unif = data.frame(CalBP=8000:3000,
                              lo=apply(densMat.unif,1,quantile,0.025),
                              hi = apply(densMat.unif,1,quantile,0.975),
                              mean = apply(densMat.unif,1,mean))

save(tbs,tbs2,tblocks.trap,tblocks.gauss,tblocks.unif,densSummary.trap.swkanto,densSummary.trap.swkanto,densSummary.unif.swkanto,densSummary.trap.chubukochi,densSummary.gauss.swkanto,densSummary.gauss.chubukochi,densSummary.unif.chubukochi,densSummary.gauss,densSummary.unif,densSummary.trap,file="./R_images/simdatesPithouses.RData")






## SPD Analysis ####
tbs = seq(7950,3050,-100)
tbs2 = seq(8000,3000,-100)

#Load C14 data image (generate by running the script in bindC14csv.R)
load("./R_images/spdC14.RData")
spdDataCal = calibrate(spdDataC14$CRA,spdDataC14$Error,normalise=FALSE,calCurves=spdDataC14$ccurve,resOffsets=spdDataC14$dR,resErrors=spdDataC14$dR) #calibrate
spdDataBin = binPrep(spdDataC14$SiteID,spdDataCal,h=200) #bin (200 years)
spdDataSPD = spd(spdDataCal,bins=spdDataBin,timeRange=c(8000,3000)) #generate SPD
spdDataSPD_blocks = spd2rc(spdDataSPD,breaks=seq(8000,3000,-100)) #aggregate by 100yrs blocks
spdDataSPD_sampled = sampleDates(x = spdDataCal,bins = spdDataBin,nsim = nsim) #sample random dates

tblocksCal =matrix(NA,nrow=length(tbs),ncol=nsim)

for (s in 1:nsim)
{
  tblocksCal[,s]=as.numeric(rev(table(cut(t(spdDataSPD_sampled$sdates)[,s],breaks=tbs2))))
}

#save image file
save(spdDataCal,spdDataBin,spdDataSPD,spdDataSPD_blocks,spdDataSPD_sampled,tblocksCal,file="./R_images/spdRes.RData")

## Comparative Analysis ####
load("./R_images/simdatesPithouses.RData")
## compute rolling correlation over 10 blocks (1,000 years)
tblockRoll10.trap = rollCor(tblocks.trap,tblocksCal,rollsize = 10) 
tblockRoll10.gauss = rollCor(tblocks.gauss,tblocksCal,rollsize = 10) 
tblockRoll10.unif = rollCor(tblocks.unif,tblocksCal,rollsize = 10) 
overallCorr.trap = overallCorr.unif = overallCorr.gauss = numeric(length=nsim)

## compute overall correlation
for (s in 1:nsim)
{
  overallCorr.trap[s] = cor(tblocks.trap[,s],tblocksCal[,s])
  overallCorr.gauss[s] = cor(tblocks.gauss[,s],tblocksCal[,s])
  overallCorr.unif[s] = cor(tblocks.unif[,s],tblocksCal[,s])
}

round(quantile(overallCorr.trap,c(0.025,0.5,0.975)),2)


# Comparing to pithouse dynamics
customModel = data.frame(calBP=8000:3000,PrDens=densSummary.trap$mean/sum(densSummary.trap$mean))

res.compare=modelTest(x=spdDataCal,bins = spdDataBin,errors = spdDataC14$Error, timeRange=c(8000,3000), changexpr=expression(100*((t1/t0)^(1/d) - 1)),model='custom',predgrid=customModel,nsim=1000,backsight=100,spdnormalised=TRUE,runm=100)


## Save output in R image file
save(tblockRoll10.trap,tblockRoll10.gauss,tblockRoll10.unif,overallCorr.trap,overallCorr.gauss,overallCorr.unif,res.compare,file="./R_images/comp.RData")





## Compute Annual % Growth Rate ####
tbs = seq(7950,3050,-100)
tbs2 = seq(8000,3000,-100)
# Compute between 5500-5400 and 5100-5000
st = which(tbs==5450)
en = which(tbs==4950)
gr.houses=100*((tblocks.trap[en,]/tblocks.trap[st,])^(1/500)-1)
round(quantile(gr.houses,c(0.025,0.5,0.975)),2)
gr.cal=100*((tblocksCal[en,]/tblocksCal[st,])^(1/500)-1)
round(quantile(gr.cal,c(0.025,0.5,0.975)),2)








