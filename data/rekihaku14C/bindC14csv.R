#### This scripts prepare the radiocarbon dates downloaded from the "Database of radiocarbon dates published in Japanese archeological research reports) (URL: https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param). 
#### Queries were carried out on the 5th of November 2019 by specifying the Prefecture (Kanagawa, Tokyo, Saitama, Nagano, and Yamanashi), and by setting "試料分類" (Material Classification) to "T: 陸産物" (Terrestial), and the "時代分類" (Period Classification) to "B:縄文時代" (Jomon). 

#### Load Packages ####
library(dplyr)


#### Read CSV Files ####
kanagawaC14 = read.csv("./kanagawa_T_B_5_11_2019.csv")
saitamaC14 = read.csv("./saitama_T_B_5_11_2019.csv")
tokyoC14 = read.csv("./tokyo_T_B_5_11_2019.csv")
naganoC14 = read.csv("./nagano_T_B_5_11_2019.csv")
yamanashiC14 = read.csv("./yamanashi_T_B_5_11_2019.csv")


#### Combine into a single data.frame
westKanto = rbind.data.frame(kanagawaC14,saitamaC14,tokyoC14,naganoC14,yamanashiC14)

####Define Unique Site Identifier
## No sites unique identifier are provided, but some sites might have different names based on the excavation. 
## Use Address as Unique Identifier
westKanto$SiteID = paste0("S",as.numeric(westKanto$所在地))

####Remove Dates on bones
westKanto=westKanto[-grep("骨",as.character(westKanto$試料の種類)),] #Remove Bones (no shells)

####Some Dates are not recorded as numeric and subset only those within 17000 and 2500 CRA
westKanto$C14年代=as.numeric(as.character(westKanto$C14年代))
westKanto$C14年代.=as.numeric(as.character(westKanto$C14年代.))
westKanto = subset(westKanto,!is.na(C14年代)&!is.na(C14年代.))
westKanto = subset(westKanto,C14年代<17000&C14年代>2500)


####Specify Prefecture
westKanto$都道府県=recode(westKanto$都道府県,"11：埼玉県"="Saitama","13：東京都"="Tokyo","14：神奈川県"="Kanagawa","19：山梨県"="Yamanashi","20：長野県"="Nagano")


####Select only relevant fields
westKantoC14 = select(westKanto,都道府県,遺跡名,SiteID,試料番号,C14年代,C14年代.)
colnames(westKantoC14) = c("Prefecture","SiteName","SiteID","LabCode","CRA","Error")                   

save(westKantoC14,file="../../R_images/westKantoC14.RData")


 


 


