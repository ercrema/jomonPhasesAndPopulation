library(magrittr)
library(dplyr)
library(oxcAAR)
library(rcabon)


## Read Data ####
c14data = read.csv("./data/c14dates.csv")

## Outlier Analysis ####
source("./R/outlierAnalysis.R")
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


## Create Oxcal Scripts ####
source("./R/oxcalScriptCreator.R")

oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/gaussian.oxcal",mcname="mcmcGaussian",model="gaussian")
oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/uniform.oxcal",mcname="mcmcUniform",model="uniform")
oxcalScriptGen(id=c14data$LabCode,c14age=c14data$CRA,errors=c14data$Error,group=c14data$CombineGroup,phases=c14data$PhaseAnalysis,fn="./oxcal/oxcalscripts/trapezoid.oxcal",mcname="mcmcTrapezoid",model="trapezoid")


## Retrieve Oxcal Output (from Oxcal Online - notice this requires about 70-90 hours of analysis) ####


## Exclude Low Agreement Index Dates and Create Oxcal Scripts ####


## Plot Posterior Phases ####

## Prepare Pithouse Data ####
source("./R/utilities.R")
phases =c("S0","S1-1","S1-2","S2-1","S2-2",paste0("S",3:8),paste0("Z",1:7),"C1","C234","C56","C78",paste0("C",9:14),paste0("K",1:8),paste0("B",1:6))
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


## Plot Results ####

# Crema 2012 re-analysis
plot(tbs.crema2012,tblocks.crema2012[,1],pch=20,col="lightgrey",xlim=c(7000,3200),ylim=range(tblocks),ylab="Number of Residential Units",xlab="cal BP",type="n")
apply(tblocks.crema2012,2,lines,x=tbs,col="lightgrey")
lines(tbs.crema2012,apply(tblocks.crema2012,1,mean))
points(tbs.crema2012,apply(tblocks.crema2012,1,mean),pch=20)
axis(1,at=seq(16000,600,-100),tck=-0.01,labels=NA)




