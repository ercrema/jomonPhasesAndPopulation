
oxcalReadjs <- function(x,model,path="")
{
  oxcalinput <- readLines(paste0(path,model,".oxcal"))
  oxcaloutput <- readLines(paste0(path,model,".js"))
  ovAg=as.numeric(strsplit(oxcaloutput[grep("model.overallAgreement",oxcaloutput)],"[=,;]")[[1]][2])
  
  print("Extracting data from oxcal output...")
  
  ii = which(is.na(x$grp)) #dates not part of a R_Combine Group
  
  pb <- txtProgressBar(min=1, max=length(ii), style=3)
  
  x$agreement = NA
  x$combine.agreement = NA
  x$pval = NA
    
  for (i in 1:length(ii))
  {
    setTxtProgressBar(pb, i)
    index=grep(paste0(x$id[ii[i]]," Posterior"),oxcaloutput)
    tmp=oxcaloutput[index]	
    position=strsplit(tmp,"[.]")[[1]][1]
    position=gsub("\\[","\\\\[",position)
    position=gsub("\\]","\\\\]",position)
    index2=grep(paste0(position,".posterior.agreement"),oxcaloutput)
    tmp2=oxcaloutput[index2]
    x$agreement[ii[i]]=as.numeric(strsplit(tmp2,"[=,;]")[[1]][2])	
  }
  close(pb)
  
  combineGroups = unique(x$grp)
  combineGroups = combineGroups[which(!is.na(combineGroups))]
  pb <- txtProgressBar(min=1, max=length(combineGroups), style=3)
  
  
  for (i in 1:length(unique(combineGroups)))
  {
    setTxtProgressBar(pb, i)
    #Extract Agreement
    index=grep(paste0('"Combination',combineGroups[i],' Posterior "'),oxcaloutput)
    tmp=oxcaloutput[index]	
    position=strsplit(tmp,"[.]")[[1]][1]
    position=gsub("\\[","\\\\[",position)
    position=gsub("\\]","\\\\]",position)
    index2=grep(paste0(position,".posterior.agreement"),oxcaloutput)
    tmp2=oxcaloutput[index2]
    x$combine.agreement[which(x$grp==combineGroups[i])] = as.numeric(strsplit(tmp2,"[=,;]")[[1]][2])	
    
    #Extract X-square Test
    index=grep(paste0('"Combination',combineGroups[i],'"'),oxcaloutput)
    tmp=oxcaloutput[index]	
    position=strsplit(tmp,"[.]")[[1]][1]
    position=gsub("\\[","\\\\[",position)
    position=gsub("\\]","\\\\]",position)
    index2=grep(paste0(position,".likelihood.testText"),oxcaloutput)
    tmp2=oxcaloutput[index2]
    df = as.numeric(gsub('^.*df=\\s*|\\s*T.*$', '', tmp2))
    t  = as.numeric(gsub('^.*T=\\s*|\\s*\\(.*$', '', tmp2))
    x$pval[which(x$grp==combineGroups[i])]=1-pchisq(t,df)
  }
  close(pb)
  
  
  print("Done.")
  return(list(ovAg=ovAg,df=x))
}
