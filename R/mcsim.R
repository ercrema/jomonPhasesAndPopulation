#'@title 
#'@description Samples a random date from relative chronologies using posterior samples from Bayesian Chrnological Models.
#'@param x A data.frame with two columns, StartPhase and EndPhase.
#'@param nsim Number of Monte-Carlo simulations
#'@param posterior A posterior class object created using \code{convertToArray()}.
#'@param weights How the probability weight should be assigned to each event in case of multiple phases. One between 'equal', 'variance', or 'custom'. Default is 'variance', and 'custom' requires a list with a vector of weights to each of the multiple phases considered (argument customWeight).
#'@example
#'source("./R/mcsim.R")
#' source("./R/rdirichlet.R")
#' 
#' load("R_images/posteriorSamples.RData")
#' 
#' 
#' nsim=1000
#' # A Sample of three PitHouse with Start/End ceramic phase
#' df=data.frame(StartPhase=c("B1","C1","Z4",),EndPhase=c("B2","C1","Z6"))
#' 
#' # Hypothetical Count Data 
#' B1=90
#' B2=67
#' Z4=120
#' Z5=1
#' Z6=2
#' 
#' df=data.frame(StartPhase=c("B1","C1","Z4","B1"),EndPhase=c("B2","C1","Z6","B2"))
#' weightList=vector('list',length(4))
#' weightList[[1]]=rdirichlet(n=nsim,alpha=c(B1,B2)+1)
#' weightList[[2]]=1
#' weightList[[3]]=rdirichlet(n=nsim,alpha=c(Z4,Z5,Z6)+1)
#' weightList[[4]]=rdirichlet(n=nsim,alpha=c(B1,B2)+1)
#' 
#' eq.prior=simTrapezoid=mcsim(df,nsim=nsim,posterior=postTrapezoid,weights="equal")
#' var.prior=simTrapezoid=mcsim(df,nsim=nsim,posterior=postTrapezoid,weights="variance")
#' empiricalbayes.prior=simTrapezoid=mcsim(df,nsim=nsim,posterior=postTrapezoid,weights="custom",weightList=weightList)
#' 
#' par(mfrow=c(3,1))
#' {
#'   hist(eq.prior[3,],main='Pithouse 3 (Z4-Z6), using equal prior',breaks=seq(5000,7000,50),xlim=c(7000,5000))
#'   hist(var.prior[3,],main='Pithouse 3 Z4-Z6), using variance (duration) prior',breaks=seq(5000,7000,50),xlim=c(7000,5000))
#'   hist(empiricalbayes.prior[3,],main='Pithouse 3 (Z4-Z6), using empricial bayes prior',breaks=seq(5000,7000,50),xlim=c(7000,5000))
#' }

mcsim = function(df,nsim,posterior,weights=c("equal","variance",'custom'),weightList=NULL,verbose=TRUE)
{
  require(trapezoid)
  require(rcarbon)
#Generate output matrix  
  resmat = matrix(nrow=nrow(df),ncol=nsim)
  
#Extract Phase Names
  phases=factor(dimnames(posterior$posterior)[[1]],levels=dimnames(posterior$posterior)[[1]],ordered=TRUE)
  
# Identify unique permutation of Start/Ends and number of cases
  permUnique = unique(df)
  permIndex = sapply(1:nrow(permUnique),function(x,p,complete){return(which(complete$StartPhase==p$StartPhase[x]&complete$EndPhase==p$EndPhase[x]))},p=permUnique,complete=df)

# Start nsim routine
  if (verbose) {pb <- txtProgressBar(min=1, max=nsim, style=3)}
  for (i in 1:nsim)
  {
    if (verbose) {setTxtProgressBar(pb, i)}
  # Select a specific posterior sample
  s = posterior[[1]][,,sample(1:dim(posterior[[1]])[3],size=1)]
  # Create vector of variance for each phase
  if(posterior$type=="trapezium"){v = apply(s,1,function(x){return(varTrap(x[1],x[2],x[3],x[4]))})}
  if(posterior$type=="uniform"){v = apply(s,1,function(x){return(varUni(x[1],x[2]))})}
  if(posterior$type=="gaussian"){v = apply(s,1,function(x){return(x[2]^2)})}
  v = sqrt(v) #convert to SD
  
  
  
  if (weights=="equal") {v=replace(v,1:length(v),1)}
      
  # Create a vector of sampled phases
  sampled.phases = rep(NA,length=nrow(df))
  
  for (p in 1:nrow(permUnique))
  {
    #If single phase just a random draw from the probability distribution
    if(permUnique[p,1]==permUnique[p,2])
    {
      sampled.phases[permIndex[[p]]]=as.character(permUnique[p,1])
    }
    
    #If multiple phase assign randomly with probabilit prob
    if (permUnique[p,1]!=permUnique[p,2])
    {
      
      timespan = which(as.character(phases)==as.character(permUnique[p,1])):which(as.character(phases)==as.character(permUnique[p,2]))
      nsamples = length(permIndex[[p]])
      
      if (weights%in%c('equal','variance'))
      { 
        timespan.weights = v[timespan]
        prob=timespan.weights/(sum(timespan.weights))
      }
      if (weights=='custom')
      {
        tmp.index=permIndex[[p]]
        if (length(tmp.index)>1){tmp.index=sample(tmp.index,size=1)}
        prob = weightList[[tmp.index]][sample(1:nsim,size=1),]
      }
      sampled.phases[permIndex[[p]]] = as.character(phases[sample(timespan,nsamples,prob=prob,replace=TRUE)])
    }
  }
    
  # Identify unique phases
  unique.phases = unique(sampled.phases)
  unique.phases.index = sapply(unique.phases,function(x,y){which(y==x)},y=sampled.phases)
  
  # Sample Dates
  for (k in 1:length(unique.phases))
  {
    n = length(unique.phases.index[[k]])
    if(posterior$type=="trapezium")
    {
      pp = s[unique.phases[k],]
      resmat[unique.phases.index[[k]],i] = (rtrapezoid(n=n,pp[1],pp[2],pp[3],pp[4]))
    }
    if(posterior$type=="uniform")
    {
      pp = s[unique.phases[k],]
      resmat[unique.phases.index[[k]],i] = (runif(n=n,pp[1],pp[2]))
    }
    if(posterior$type=="gaussian")
    {
      pp = s[unique.phases[k],]
      resmat[unique.phases.index[[k]],i] = (rnorm(n=n,pp[1],pp[2]))
    }
    resmat[unique.phases.index[[k]],i] = BCADtoBP(resmat[unique.phases.index[[k]],i])
  }
    
  }
  if (verbose) {close(pb)}
  return(resmat)      
}

  
