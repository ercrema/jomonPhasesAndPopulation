oxcalScriptGen = function(id,c14age,errors,group,phases,fn,interval=100,mcnsim=5000,mcname,model=c("gaussian","uniform","trapezoid"))
{
  export <- file(fn) #create export file
  phaseNames =  unique(phases)
  nphases = length(phaseNames)
  
  
  cat("Plot(){\n",file=fn,append=FALSE) #Start Sequence#
  cat('Outlier_Model("",N(0,2),0,"s");\n',file=fn,append=TRUE)
  cat("Phase(){\n",file=fn,append=TRUE) #Start Sequence#
  
  for (k in 1:nphases)
  {
    i = which(phases == phaseNames[k])
    c14.phase = c14age[i]
    errors.phase = errors[i]
    id.phase = id[i]
    group.phase = group[i]
    
    
    # Phase Boundary Start####
    cat("Sequence(){\n", file = fn, append = TRUE) #Start Sequence#
    if (model=="gaussian")
    {
     cat(paste0('Sigma_Boundary("Start ', phaseNames[k], '");\n'), file = fn, append = TRUE)
    }
    if (model=="uniform")
    {
      cat(paste0('Boundary("Start ', phaseNames[k], '");\n'), file = fn, append = TRUE)
    }
    if (model=="trapezoid")
    {
      cat(paste0('Boundary("Start ', phaseNames[k], '"){\n'), file = fn, append = TRUE)
      cat(paste0('Start("Start of Start ', phaseNames[k], '");\n'), file = fn, append = TRUE)
      cat(paste0('Transition("Period of Start ', phaseNames[k], '");\n'), file = fn, append = TRUE)
      cat(paste0('End("End of Start ', phaseNames[k], '");\n'), file = fn, append = TRUE)
      cat('};\n', file = fn, append = TRUE)
    }
    
    
    # Actual Dates#####
    cat(paste0('Phase("', phaseNames[k], '")\n'), file = fn, append = TRUE)
    cat('{\n', file = fn, append = TRUE)
    
    #start with dates to be combined
    cb = which(!is.na(group.phase))
    
    if (length(cb)>0) 
    {
      combos = unique(group.phase[cb])
      
      for (g in 1:length(combos))
      {
        cat(paste0('R_Combine("Combination',combos[g],'")\n'),file=fn,append=TRUE)
        cat('{\n',file=fn,append=TRUE)
        
        ### Continue from here should loop through all dates within a combine group
        gg = which(group.phase==combos[g])
        
        for (ggg in 1:length(gg))
        {
          cat(paste('R_Date(','\"',id.phase[gg[ggg]],'\",',c14.phase[gg[ggg]],',',errors.phase[gg[ggg]],'){Outlier(0.05);};\n', sep = ""),file = fn, append = TRUE)
        } 
        
        cat('};\n',file=fn,append=TRUE)
      }
      
      id.phase = id.phase[-cb]
      c14.phase = c14.phase[-cb]
      errors.phase = errors.phase[-cb]
    }
    
    #then all other dates
    for (x in 1:length(c14.phase))
    {
      cat(paste('R_Date(','\"',id.phase[x],'\",',c14.phase[x],',',errors.phase[x],');\n', sep = ""),file = fn, append = TRUE)
    }
    cat('};\n', file = fn, append = TRUE)
    
    # Phase Boundary End####
    if (model=="gaussian")
    {
      cat(paste0('Sigma_Boundary("End ', phaseNames[k], '");\n'), file = fn, append = TRUE)
    }
    if (model=="uniform")
    {
      cat(paste0('Boundary("End ', phaseNames[k], '");\n'), file = fn, append = TRUE)
    }
    if (model=="trapezoid")
    {
      cat(paste0('Boundary("End ', phaseNames[k], '"){\n'), file = fn, append = TRUE)
      cat(paste0('Start("Start of End ', phaseNames[k], '");\n'), file = fn, append = TRUE)
      cat(paste0('Transition("Period of End ', phaseNames[k], '");\n'), file = fn, append = TRUE)
      cat(paste0('End("End of End ', phaseNames[k], '");\n'), file = fn, append = TRUE)
      cat('};\n', file = fn, append = TRUE)
    }
    cat('};\n', file = fn, append = TRUE)
  }
  cat('};\n', file = fn, append = TRUE) 
  
  
  # MCMC Samples
  cat(paste("MCMC_Sample(\"", mcname, "\",", interval, ",", mcnsim, "){\n", sep = ""),file = fn,append = TRUE)
  
  if (model=="gaussian"|model=="uniform")
  {
    for (k in 1:nphases)
    {
      cat(paste0('Date("=Start ', phaseNames[k], '");\n'), file = fn, append = TRUE)
      cat(paste0('Date("=End ', phaseNames[k], '");\n'), file = fn, append = TRUE)
    }
  }
  
  if (model=="trapezoid")
  {
    for (k in 1:nphases)
    {
      cat(paste0('Date("=Start of Start ', phaseNames[k], '");\n'),file = fn,append = TRUE)
      cat(paste0('Date("=End of Start ', phaseNames[k], '");\n'),file = fn,append = TRUE)
      cat(paste0('Date("=Start of End ', phaseNames[k], '");\n'),file = fn,append = TRUE)
      cat(paste0('Date("=End of End ', phaseNames[k], '");\n'), file = fn, append = TRUE)
    }
  }
  cat('};\n', file = fn, append = TRUE)
  cat('};\n', file = fn, append = TRUE)
  
  
  close(export)
}




