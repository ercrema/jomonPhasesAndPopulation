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

