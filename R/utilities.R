orgTable = function(x)
{
  n = sum(x$counts) #compute total number of houses
  StartPhase = character(0)
  EndPhase = character(0)
  Prefecture = character(0)
  x$st=as.character(x$st)
  x$en=as.character(x$en)
  x$pref=as.character(x$pref)
  
  
  for (i in 1:nrow(x))
  {
    StartPhase=c(StartPhase,rep(x$st[i],x$counts[i]))
    EndPhase=c(EndPhase,rep(x$en[i],x$counts[i]))
  }
  res = data.frame(StartPhase=StartPhase,EndPhase=EndPhase,Prefecture=x$pref[1])
  return(res)
}

#'@title Extracts Variance of a Uniform distribution from model parameters
varUni = function(a,b)
{
  return((1/12)*(b-a)^2)
}


#'@title Extracts Variance of a Trapezoid distribution from model parameters
#'@references 
#'Dorp, J. R. van, & Kotz, S. (2003). Generalized trapezoidal distributions. Metrika, 58(1), 85–97. doi:10.1007/s001840200230



varTrap = function(a,b,c,d)
{
  #Ensure that a<b<c<d
  a = a-0.2
  b = b-0.1
  d = d+0.1
  
  #first moment (mean)
  E=((b-a)*(a+2*b)-3*(b^2-c^2)+(d-c)*(2*c+d))/(3*(d+c-b-a))
  
  #second moment
  E2=((b-a)/(d+c-b-a))*((1/6)*(a+b)^2+(1/3)*b^2) +
    ((1/(d+c-b-a))*((2/3)*(c^3-b^3))) + 
    (((d-c)/(d+c-b-a))*((1/3)*c^2+(1/6)*(c+d)^2))
  
  #Variance
  return((E2-E^2))
}




#'@title Convert Sample Posterior to Array Object with matching phase names and type of posterior

convertToArray= function(x,type,phases)
{
  x = x[which(!apply(x,1,anyNA)),] #subset to posterior with actual samples
  if (type=="uniform")
  {
    res = array(NA,dim=c(length(phases),2,nrow(x)),dimnames = list(phases,c("start","end"),1:nrow(x)))
    for (i in 1:length(phases))
    {
      st=which(colnames(x)==paste0("Start.",phases[i]))
      en=which(colnames(x)==paste0("End.",phases[i]))
      res[i,1,]=x[,st]
      res[i,2,]=x[,en]
    }
  }
  if (type=="gaussian")
  {
    res = array(NA,dim=c(length(phases),2,nrow(x)),dimnames = list(phases,c("mean","sd"),1:nrow(x)))
    for (i in 1:length(phases))
    {
      st=which(colnames(x)==paste0("Start.",phases[i]))
      en=which(colnames(x)==paste0("End.",phases[i]))
      res[i,1,]=x[,st]-c(x[,st]-x[,en])/2
      res[i,2,]=abs(x[,en]-x[,st])/2
    }
    
  }
  if (type=="trapezium")
  {
    res = array(NA,dim=c(length(phases),4,nrow(x)),dimnames = list(phases,c("a","b","c","d"),1:nrow(x)))
    for (i in 1:length(phases))
    {
      a=which(colnames(x)==paste0("Start.of.Start.",phases[i]))
      b=which(colnames(x)==paste0("End.of.Start.",phases[i]))
      c=which(colnames(x)==paste0("Start.of.End.",phases[i]))
      d=which(colnames(x)==paste0("End.of.End.",phases[i]))
      
      #Add/remove 1/3 to avoid synchrony 
      res[i,1,]=x[,a] - 2/3
      res[i,2,]=x[,b] - 1/3
      res[i,3,]=x[,c] + 1/3
      res[i,4,]=x[,d] + 2/3
    }
  }
  return(list(posterior=res,type=type))
}


error.bar <- function(x,upper, lower=upper, length=0.1,...){
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}


rollCor = function(x,y,rollsize=10)
{
 require(TTR)
 rMatrix = matrix(NA,nrow(x),ncol(x))
 for (s in 1:ncol(x))
 {
   rMatrix[,s] = runCor(x[,s],y[,s],n=rollsize)
 }
 return(rMatrix)
}






