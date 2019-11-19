# Read data from each prefecture
nagano = read.csv("./nagano.csv",stringsAsFactors = FALSE)
kanagawa = read.csv("./kanagawa.csv",stringsAsFactors = FALSE)
yamanashi = read.csv("./yamanashi.csv",stringsAsFactors = FALSE)
tokyo = read.csv("./tokyo.csv",stringsAsFactors = FALSE)
saitama = read.csv("./saitama.csv",stringsAsFactors = FALSE)

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
pithouseData$Region = "SWKanto"
pithouseData$Region[which(pithouseData$Prefecture%in%c("Yamanashi","Nagano"))]="ChubuHighlands"
# save data into an R image
save(pithouseData,file="../../R_images/pithouseData.RData")
