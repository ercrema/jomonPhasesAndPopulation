outlierProb = function(c14,errors,id=NA,oxpath='~/',fn=NA)
{
  require(oxcAAR)
  require(readr)

  if (anyNA(id))
  {
    id = 1:length(c14)
  }
  
  if(is.na(fn))
  {
    fn=tempfile()
  }	
  export <- file(fn) #create export file
  
  cat("Plot(){\n",file=fn,append=FALSE) #Start Sequence#
  cat('Outlier_Model("",N(0,2),0,"s");\n',file=fn,append=TRUE)
  cat('R_Combine("")\n',file=fn,append=TRUE)
  cat('{\n',file=fn,append=TRUE)
  for (x in 1:length(c14))
  {
    cat(paste('R_Date(','\"',id[x],'\",',c14[x],',',errors[x],'){Outlier(0.05);};\n',sep=""),file=fn,append=TRUE)
  }
  cat('};\n',file=fn,append=TRUE)
  cat('};\n',file=fn,append=TRUE)

  quickSetupOxcal(path=oxpath)
  oxcalscript <- read_file(fn)
  close(export)
  result_file <- executeOxcalScript(oxcalscript)
  oxcaloutput <- readLines(result_file)
  
  index=grep("outlier_post",oxcaloutput)
  tmp=oxcaloutput[index]
  outlier.post=as.numeric(sapply(strsplit(tmp,"[=,;]"),function(x) x[2]))[-1] #extract outlier removing the first value
  
  index2=grep(".likelihood.testText",oxcaloutput)
  tmp2=oxcaloutput[index2]
  df = as.numeric(gsub('^.*df=\\s*|\\s*T.*$', '', tmp2))
  t  = as.numeric(gsub('^.*T=\\s*|\\s*\\(.*$', '', tmp2))
  pval=1-pchisq(t,df)
  
  ovAg=as.numeric(strsplit(oxcaloutput[grep("posterior.overallAgreement=",oxcaloutput)],"[=,;]")[[1]][2])
  
  return(list(agreement=ovAg,pval=pval,df=data.frame(id=id,c14=c14,errors=errors,outlier.prob=outlier.post)))
}



outlierExcluder = function(c14,errors,id=NA,oxpath='~/',fn=NA)
{
	require(dplyr)
	require(magrittr)
	if (anyNA(id)) {id = 1:length(c14)}
	df=data.frame(id=id,c14=c14,errors=errors,exclude=FALSE)
	x = outlierProb(c14=c14,errors=errors,id=id,oxpath=oxpath,fn=fn)
	check = (x$pval<0.05|x$agreement<60)
	tmp=x$df

	if (!check)
	{
		df = data.frame(id=id,c14=c14,errors=errors,outlier.prob=NA,exclude=FALSE)
	}

	if (check&(nrow(tmp)>2))
	{
		while(check)
		{
			i = order(tmp$outlier.prob,decreasing = TRUE)[1]
			df$exclude[which(df$id==tmp$id[i])]=TRUE
			c14=subset(df,!exclude)$c14
			errors=subset(df,!exclude)$errors
			id=subset(df,!exclude)$id

			x = outlierProb(c14=c14,errors=errors,id=id,oxpath=oxpath,fn=fn)
			check = (x$pval<0.05|x$agreement<60)
			tmp=x$df
		}

		df= left_join(df,x$df,by='id') %>% select(id,c14=c14.x,errors=errors.x,outlier.prob,exclude)
	} 
	
	if (check&nrow(tmp)==2)
       	{
		df = data.frame(id=id,c14=c14,errors=errors,outlier.prob=NA,exclude=NA)
	}


	return(list(finaloutput=x,df=df))
}

