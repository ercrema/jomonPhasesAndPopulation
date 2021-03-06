#### This scripts prepare the radiocarbon dates downloaded from the "Database of radiocarbon dates published in Japanese archeological research reports) (URL: https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param). 
#### Queries were carried out on the 5th of November 2019 and 11th of March 2020 by specifying the Prefecture (Kanagawa, Tokyo, Saitama, Nagano, and Yamanashi), and by setting "試料分類" (Material Classification) to "T: 陸産物" (Terrestial) or "M: 海産物" (Marine), and the "時代分類" (Period Classification) to "B:縄文時代" (Jomon). 
#### Downloaded CSV files were originally encoded in Shift-JIS and re-encoded in UTF-8

#### Load Packages ####
library(dplyr)


#### Read CSV Files ####
kanagawaC14 = read.csv("./kanagawa_T_B_5_11_2019.csv",stringsAsFactors = FALSE)
saitamaC14 = read.csv("./saitama_T_B_5_11_2019.csv",stringsAsFactors = FALSE)
tokyoC14 = read.csv("./tokyo_T_B_5_11_2019.csv",stringsAsFactors = FALSE)
naganoC14 = read.csv("./nagano_T_B_5_11_2019.csv",stringsAsFactors = FALSE)
yamanashiC14 = read.csv("./yamanashi_T_B_5_11_2019.csv",stringsAsFactors = FALSE)

#### Marine Dates (no instances for nagano and yamanashi)
kanagawaC14m = read.csv("./kanagawa_M_B_11_3_2020.csv",stringsAsFactors = FALSE)
saitamaC14m = read.csv("./saitama_M_B_11_3_2020.csv",stringsAsFactors = FALSE)
tokyoC14m = read.csv("./tokyo_M_B_11_3_2020.csv",stringsAsFactors = FALSE)

#### Combine into a single data.frame
spdData.terrestrial = rbind.data.frame(kanagawaC14,saitamaC14,tokyoC14,naganoC14,yamanashiC14)
spdData.marine=rbind.data.frame(tokyoC14m,saitamaC14m,kanagawaC14m)
spdData.terrestrial$Marine=FALSE
spdData.marine$Marine=TRUE
spdData=rbind.data.frame(spdData.terrestrial,spdData.marine)

####Define Unique Site Identifier
## No sites unique identifier are provided, but some sites might have different names based on the excavation. 
## Use Address as Unique Identifier
spdData$SiteID = paste0("S",as.numeric(as.factor(spdData$所在地)))

## Eliminate Dates without a Labcode
spdData = subset(spdData,試料番号!="不明") #"不明" is unknown in Japanese

#### Eliminate Dates with the same Lab-Code
spdData$retain=TRUE
any(duplicated(spdData$試料番号)) # Presence of Duplicated Lab-Codes
duplicate.counts = table(spdData$試料番号)
duplicate.counts = duplicate.counts[which(duplicate.counts>1)]
duplicate.names = names(duplicate.counts)

for (i in 1:length(duplicate.names))
{
  tmp.index=which(spdData$試料番号==duplicate.names[i])
  if (length(unique(spdData$C14年代[tmp.index]))>1 |
      length(unique(spdData$C14年代.[tmp.index]))>1 |
      length(unique(spdData$SiteID[tmp.index]))>1)
  {
    spdData$retain[tmp.index]=FALSE
  } else {
    spdData$retain[sample(tmp.index,size=length(tmp.index)-1)]=FALSE
  }
}

spdData = subset(spdData,retain==TRUE)

####Some Dates are not recorded as numeric 
spdData$C14年代=as.numeric(as.character(spdData$C14年代))
spdData$C14年代.=as.numeric(as.character(spdData$C14年代.))
spdData = subset(spdData,!is.na(C14年代)&!is.na(C14年代.))

#### Subset only those within 7500 and 2500 CRA
spdData = subset(spdData,C14年代<7500&C14年代>2500)


####Specify Prefecture
spdData$都道府県=recode(spdData$都道府県,"11：埼玉県"="Saitama","13：東京都"="Tokyo","14：神奈川県"="Kanagawa","19：山梨県"="Yamanashi","20：長野県"="Nagano")

####Remove Dates on bones
length(grep("骨",as.character(spdData$試料の種類))) #only 13 cases
spdData=spdData[-grep("骨",as.character(spdData$試料の種類)),] #Remove Bones (no shells)

####Select only relevant fields
spdDataC14 = select(spdData,都道府県,遺跡名,SiteID,試料番号,C14年代,C14年代.,Marine)
colnames(spdDataC14) = c("Prefecture","SiteName","SiteID","LabCode","CRA","Error","Marine")           
spdDataC14$ccurve='intcal20'  

#### Add Reservoir Effect (following Shishikura et al 2007, DOI:10.1016/j.yqres.2006.09.003)
#### DeltaR calculated for Marine20 with weighted average of all samples
spdDataC14$ccurve[which(spdDataC14$Marine==TRUE)]='marine20'
spdDataC14$dR=0
spdDataC14$dRe=0
spdDataC14$dR[which(spdDataC14$Marine==TRUE)]=-21
spdDataC14$dRe[which(spdDataC14$Marine==TRUE)]=76

save(spdDataC14,file="../../R_images/spdC14.RData")


 


 


