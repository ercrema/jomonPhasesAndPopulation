#'@title 
#'@description Samples a random date from relative chronologies using posterior samples from Bayesian Chrnological Models.
#'@param x A data.frame with two columns, StartPhase and EndPhase.
#'@param nsim Number of Monte-Carlo simulations
#'@param posterior A posterior class object created using \code{convertToArray()}.
#'@param weights How the probability weight should be assigned to each event in case of multiple phases. One between 'equal', 'variance', or 'custom'. Default is 'variance', and 'custom' requires a list with a vector of weights to each of the multiple phases considered (argument customWeight).
#'@example
#' load("/Users/enryu/github/ENCOUNTER_github/encounter_jomonphases/results/model0_oxcal.RData")
#' df=data.frame(StartPhase=c("B1","C1","C234","Z4","C234","C234"),EndPhase=c("B2","C1","C56","Z6","C234","C56"),stringsAsFactors = FALSE)
#' posterior.a=convertToArray(model0a,"gaussian")
#' result=sampler(df,nsim=100,posterior=posterior.a,weights="equal")

mcsim = function(df,nsim,posterior,weights=c("equal","variance"))
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
  for (i in 1:nsim)
  {
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
      if (permUnique[p,1]!=permUnique[p,2])
      {
        timespan = which(as.character(phases)==as.character(permUnique[p,1])):which(as.character(phases)==as.character(permUnique[p,2]))
        timespan.weights = v[timespan]
        nsamples = length(permIndex[[p]])
        sampled.phases[permIndex[[p]]] = as.character(phases[sample(timespan,nsamples,prob=timespan.weights/(sum(timespan.weights)),replace=TRUE)])
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
  
  return(resmat)      
}

  
