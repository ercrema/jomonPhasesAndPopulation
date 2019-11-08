library(magrittr)
library(dplyr)
library(oxcAAR)
library(rcarbon)

source("./R/utilities.R")
source("./R/oxcalReadjs.R")
source("./R/oxcalScriptCreator.R")
source("./R/outlierAnalysis.R")
source("./R/mcsim.R")


## Read Data ####
c14data = read.csv("./data/c14dates.csv")

## Outlier Analysis ####
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
save.image("./R_images/c14data.RData")

## Generate Summary Table 
phases =c("S0","S1.1","S1.2","S2.1","S2.2",paste0("S",3:8),paste0("Z",1:7),"C1","C234","C56","C78",paste0("C",9:14),paste0("K",1:8),paste0("B",1:6))
kobayashiPhases = c("S0","S1-1","S1-2","S2-1","S2-2","S3-1;S3-2","S4","S5","S6","S7","S8","Z1","Z2","Z3","Z4","Z5","Z6","Z7","C1","C2;C3;C4","C5;C6","C7;C8","C9a;C9bc","C10a;C10b;C10c","C11ab;C11c","C12a;C12b;C12c","C13a;C13b","C14a;C14b","K1-1;K1-2;K1-3","K2","K3","K4","K5","K6","K7","K8-1;K8-2","B1","B2","B3","B4","B5","B6")


table1 = data.frame(phases,kobayashi2017=kobayashiPhases,samples=NA,effsamples=NA,nsites=NA,correspondingphases=NA,correspondingphasesEn=NA)
conversion = read.csv("./data/suzuki/conversion.csv",stringsAsFactors = FALSE)
conversionList = vector("list",length=nrow(conversion))

for (i in 1:length(phases))
{
  tmp = subset(c14data,PhaseAnalysis==phases[i])
  table1$samples[i] = nrow(tmp)
  table1$effsamples[i] = sum(is.na(tmp$CombineGroup)) + length(unique(tmp$CombineGroup[which(!is.na(tmp$CombineGroup))]))
  table1$nsites[i] = length(unique(tmp$SiteID))
  k = which(conversion$Interval_Start==phases[i]|conversion$Interval_End==phases[i])
  table1$correspondingphases[i]=paste(conversion$PotteryPhase_Suzuki1986[k],collapse=";")
  table1$correspondingphasesEn[i]=paste(conversion$PotteryPhase_Suzuki1986_En[k],collapse=";")
}

write.csv(table1,file="./manuscript/tables/table1_base.csv")





## Create Oxcal Scripts ####

oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/gaussian.oxcal",mcname="mcmcGaussian",model="gaussian")
oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/uniform.oxcal",mcname="mcmcUniform",model="uniform")
oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/trapezoid.oxcal",mcname="mcmcTrapezoid",model="trapezoid")


## Retrieve Oxcal Output (from Oxcal Online - notice this requires about 70-90 hours of analysis) ####

df = data.frame(id=as.character(c14data$LabCode),grp=c14data$CombineGroup,stringsAsFactors = FALSE)

# Read js files (this takes some time)
gaussian.agreement = oxcalReadjs(x=df, model='gaussian',path='./oxcal/results/')
uniform.agreement = oxcalReadjs(x=df, model='uniform',path='./oxcal/results/')
trapezoid.agreement = oxcalReadjs(x=df, model='trapezoid',path='./oxcal/results/')

c14data.gaussian.rerun = left_join(c14data,gaussian.agreement$df,by=c("LabCode"="id")) %>%
  subset(agreement>60|combine.agreement>60)
c14data.uniform.rerun = left_join(c14data,uniform.agreement$df,by=c("LabCode"="id")) %>%
  subset(agreement>60|combine.agreement>60)
c14data.trapezoid.rerun = left_join(c14data,trapezoid.agreement$df,by=c("LabCode"="id")) %>%
  subset(agreement>60|combine.agreement>60)

# Resubmission to OxCal 

oxcalScriptGen(id=c14data.gaussian.rerun$LabCode,c14age=c14data.gaussian.rerun$CRA,errors=c14data.gaussian.rerun$Error,group=c14data.gaussian.rerun$CombineGroup,phases=c14data.gaussian.rerun$PhaseAnalysis,fn="./oxcal/oxcalscripts/gaussianR.oxcal",mcname="mcmcGaussianR",model="gaussian")
oxcalScriptGen(id=c14data.uniform.rerun$LabCode,c14age=c14data.uniform.rerun$CRA,errors=c14data.uniform.rerun$Error,group=c14data.uniform.rerun$CombineGroup,phases=c14data.uniform.rerun$PhaseAnalysis,fn="./oxcal/oxcalscripts/uniformR.oxcal",mcname="mcmcUniformR",model="uniform")
oxcalScriptGen(id=c14data.trapezoid.rerun$LabCode,c14age=c14data.trapezoid.rerun$CRA,errors=c14data.trapezoid.rerun$Error,group=c14data.trapezoid.rerun$CombineGroup,phases=c14data.trapezoid.rerun$PhaseAnalysis,fn="./oxcal/oxcalscripts/trapezoidR.oxcal",mcname="mcmcTrapezoidR",model="trapezoid")

# Read re-run results ####

gaussian.agreement = oxcalReadjs(x=df, model='gaussianR',path='../oxcal/results/')
uniform.agreement = oxcalReadjs(x=df, model='uniformR',path='../oxcal/results/')
trapezoid.agreement = oxcalReadjs(x=df, model='trapezoidR',path='../oxcal/results/')

# Read MCMC samples (excluding first (pass number) and last (empty) column)
gaussian.samples = read.csv("../oxcal/oxcalOnlineResults/mcmcGaussianR.csv")[,-c(1,86)] 
uniform.samples = read.csv("../oxcal/oxcalOnlineResults/mcmcUniformR.csv")[,-c(1,86)]
trapezoid.samples = read.csv("../oxcal/oxcalOnlineResults/mcmcTrapezoidR.csv")[,-c(1,170)]

# Convert into an Array
phases =c("S0","S1.1","S1.2","S2.1","S2.2",paste0("S",3:8),paste0("Z",1:7),"C1","C234","C56","C78",paste0("C",9:14),paste0("K",1:8),paste0("B",1:6))
postGaussian=convertToArray(gaussian.samples,type="gaussian",phases)
postUniform=convertToArray(uniform.samples,type="uniform",phases)
postTrapezoid=convertToArray(trapezoid.samples,type="trapezium",phases)


## Prepare Pithouse Data ####
nagano = read.csv("./data/suzuki/nagano.csv",stringsAsFactors = FALSE)
kanagawa = read.csv("./data/suzuki/kanagawa.csv",stringsAsFactors = FALSE)
yamanashi = read.csv("./data/suzuki/yamanashi.csv",stringsAsFactors = FALSE)
tokyo = read.csv("./data/suzuki/tokyo.csv",stringsAsFactors = FALSE)
saitama = read.csv("./data/suzuki/saitama.csv",stringsAsFactors = FALSE)

pthlist=list(Nagano=nagano,Kanagawa=kanagawa,Tokyo=tokyo,Saitama=saitama,Yamanashi=yamanashi)
res=vector("list",length=length(pthlist))

for (i in 1:length(pthlist))
{
  tmp = pthlist[[i]]
  cases = unique(data.frame(st=tmp$Timespan_Start,en=tmp$Timespan_End))
  ncases = nrow(cases)
  cases$counts = NA
  for (j in 1:ncases)
  {
    tmp2=subset(tmp,Timespan_Start==cases$st[j]&Timespan_End==cases$en[j])
    cases$counts[j]=sum(tmp2$Counts)
  }
  cases$pref=names(pthlist)[[i]]
  res[[i]]=cases
  names(res)[[i]]=names(pthlist)[[i]]
}

res=lapply(res,orgTable) #orgTable converts aggregated counts into a data.frame with 1 house per row
pithouseData=rbind.data.frame(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]]) #combine to a single data.frame

# Simulate Pithouse Dates ####
nsim = 1000
simGaussian=mcsim(pithouseData,nsim=1000,posterior=postGaussian,weights="variance")
simUniform=mcsim(pithouseData,nsim=1000,posterior=postUniform,weights="variance")
simTrapezoid=mcsim(pithouseData,nsim=1000,posterior=postTrapezoid,weights="variance")

# count for each 100-year block
posteriors = list(gaussian=simGaussian,uniform=simUniform,trapezoid=simTrapezoid)
tblocks=vector("list",length=3)
tbs=seq(14950,650,-100)
tblocks[[1]] = tblocks[[2]] = tblocks[[3]] = matrix(NA,nrow=144,ncol=nsim)

for (x in 1:3)
{
  for (s in 1:1000)
  {
    tblocks[[x]][,s]=as.numeric(rev(table(cut(posteriors[[x]][,s],breaks=seq(15000,600,-100)))))
  }
}


# Comparison Analysis with Crema 2012
pithouseData_compare = subset(pithouseData,Prefecture!="Nagano" & !StartPhase%in%c("S0","S1.1","S1.2","S2.1","S2.2",paste0("S",3:8),paste0("B",1:6)) & !EndPhase%in%c("S0","S1.1","S1.2","S2.1","S2.2",paste0("S",3:8),paste0("B",1:6)))

# Crema 2012 Re-analysis ####
# There is a small discrepancies in the results between Crema 2012 and this re-analysis.
crema2012_phases = read.csv("./data/suzuki/crema2012_phases.csv",stringsAsFactors = FALSE)
crema2012_counts = read.csv("./data/suzuki/crema2012_counts.csv",stringsAsFactors = FALSE)
crema2012_combined=left_join(crema2012_counts,crema2012_phases,by=c("TimeSpanStart"="Phase")) %>%
  select(TimeSpanStart,TimeSpanEnd,Count,st=Start) %>%
  left_join(crema2012_phases,by=c("TimeSpanEnd"="Phase")) %>%
  select(TimeSpanStart,TimeSpanEnd,Count,st,en=End)
  
nsim = 1000
sim.crema2012=replicate(n = nsim,unlist(apply(crema2012_combined,1,function(x){return(runif(n=x[1],max=x[2],min=x[3]))})))

tbs.crema2012=seq(6950,3250,-100)
tblocks.crema2012=matrix(NA,nrow=length(tbs.crema2012),ncol=nsim)

for (s in 1:1000)
{
  tblocks.crema2012[,s]=as.numeric(rev(table(cut(sim.crema2012[,s],breaks=seq(7000,3200,-100)))))
}


# SPD Analysis
load("./R_images/westKantoC14.RData")
westKantoCal = calibrate(westKantoC14$CRA,westKantoC14$Error,normalise=FALSE)
westKantoBin = binPrep(westKantoC14$SiteID,westKantoC14$CRA,h=200)
westKantoSPD = spd(westKantoCal,timeRange=c(7000,3200))
westKantoSPD2 = spd2rc(westKantoSPD,breaks=seq(7000,3200,-100))
plot(seq(6950,3250,-100),westKantoSPD2$sumblock,type="h",xlim=c(7000,3200))

## Plot Results ####

# Crema 2012 re-analysis
plot(tbs.crema2012,tblocks.crema2012[,1],pch=20,col="lightgrey",xlim=c(7000,3200),ylim=range(tblocks),ylab="Number of Residential Units",xlab="cal BP",type="n")
apply(tblocks.crema2012,2,lines,x=tbs,col="lightgrey")
lines(tbs.crema2012,apply(tblocks.crema2012,1,mean))
points(tbs.crema2012,apply(tblocks.crema2012,1,mean),pch=20)
axis(1,at=seq(16000,600,-100),tck=-0.01,labels=NA)




