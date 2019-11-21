# Bayesian phase models, residential frequency data, and multi-proxy inferences of Jomon population dynamics: source code, data, and scripts

This repository contains all data and scripts required to fully reproduce all analyses presented in the paper _Bayesian phase models, residential frequency data, and multi-proxy inferences of Jomon population dynamics_ authored by Crema, E.R., and Kobayashi, K.

The main workflow is recorded in the [log.R](./log.R) file and outputs are stored as R image files located in the [R_images](./R_images) directory.   

## Data Sets and Data Preparation

All raw data used in the paper can be found in the [data](./data) directory. The file [c14dates.csv](./data/c14dates.csv) contains the radiocarbon dates used for the ceramic phase modelling, the folder [rekihaku14C](./data/rekihaku14C) contains CSV files of radiocarbon dates downloaded from the National Museum of Japanese History's [Database of radiocarbon dates published in Japanese archaeological research reports](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param). The R script file [bindC14csv.R](./data/rekihaku14C/bindC14csv.R) in the same directory contains details of the query used to download the data, as well as a script for aggregating the files into an R image file (`westKantoNaganoC14.RData`) contained in the [R_images](./R_images) directory. Te sub-directory [./data/suzuki](./data/suzuki) contains CSV files of Jomon pithouse counts obtained from tables on pages 88 to 93 of 「縄文時代集落の研究」("Research on Jomon Period Settlements"), by 鈴木 保彦  (Suzuki, Yasuhiko), published by 雄山閣 (Yuzankaku), Tokyo, in 2006. The data has been digitised into separate CSV files and aggregated into a long format data.frame containing a pithouse per row. This process is recorded in the R script [pithouseBinder.R](./data/suzuki/pithouseBinder.R) and the resulting data.frame is stored in the R image file [.](./R_images/pithouseData.RData).

## Bayesian ceramic phase modelling

Bayesian modelling have been carried out using [OxCal](https://c14.arch.ox.ac.uk/oxcal/OxCal.html), with the preparation of OxCal scripts done in R. The analyses was conducted in three stages. Firstly, potential outliers were detected within sets of dates associated with the same event (e.g. different organic residues from the same vessel). This was achievied by using the `outlierExcluder()` function (stored [here](./R/outlierAnalysis.R)) which internally calls OxCal from R using the [oxcAAR](https://CRAN.R-project.org/package=oxcAAR) package. The function eliminates potential outliers from the initial set of radiocarbon dates and result of this routine was stored in the R image file [c14data.RData](./R_images/c14data.RData). 

The second step of analysis consisted in creating OxCal scripts for different probability distributions emulating putative _within-phase uncertainty_. This was achieved by using the `oxcalScriptGen()` function (located [in this file](./R/oxcalScriptCreator.R). The output scripts (stored in the directory [./oxcal/oxcalscripts](./oxcal/oxcalscripts)) were then loaded into OxCal for analyses. The results of the Bayesian analyses are stored in the directory [./oxcal/results](./oxcal/results), and include the .csv storing the posterior samples and a JavaScript file (.js) containing key statistics such as individual and overall agreement indices, read in R using the `oxcalReadjs()` function (see source code [here](./R/oxcalReadjs.R)). 

The third and final step consisted of removing samples with low agreement index and generate a new set of OxCal scripts (with suffix "R" to distinguish this second submission to OxCal, e.g. `gaussian.oxcal` and `gaussianR.oxcal`). The result of the OxCal analysis was then processed and the posterior samples organised into a series of objects stored in the R image [posteriorSamples.RData](./R_images/posteriorSamples.RData).

## Monte-Carlo Simulation of Pithouse Dates

Dates of individual pithouse were simulated 5,000 times taking account: 1) the uncertatinty within the phase ( _within phase uncertainty_ ); 2) the uncertainty in defining the membership of the pithouse to a particular phase ( _phase assignement uncertainty_ ); and the uncertainty associated with the parameters of the probability distribution describing the _within phase uncertainty_  (i.e. the _phase boundary uncertainty_). The function `mcsim()` (see source code [here](./R/mcsim.R)) was used for this purpose. The R image fle [simdatesPithouses.RData](./R/simdatesPithouses.RData) contains the 5,000 set of simulated dates, along with counts organised in 100 years bins (between 8,000 and 3,000 cal BP), and the outcome of a composite kernel density estimate analyses (see details in the [log.R](./log.R) file).   

## SPD Analysis

Summed probability distribution (SPD) of radiocarbon dates have been generated using [rcarbon](https://CRAN.R-project.org/package=rcarbon) <!-- specify whether this 1.2 or 1.3 depending on the radiocarbon paper -->. To enable correlation analyses a matrix of 5000 sets of randomly sampled calendar dates was created and aggregated by 100 years intervals between 8,000 and 3,000 cal BP. All outputshave been stored in the R image [./R_images/spdRes.RData](./R_images/spdRes.RData).  


## Comparisons between pithouse data and radiocarbon density

The time-series of residential and radiocarbon density have been compared via correlation analyses and the `modelTest()` function in [rcarbon](https://CRAN.R-project.org/package=rcarbon). The former was carried out by generating 5,000 correlation values by iteratively comparing the the time-series of simulated pithouse dates and randomly sampled calendar dates from the calibrated radiocarbon dates. The latter compared the observed annual growth rate in the SPD against an expectation derived from the average trend obtained from the composite kernel density estimate of pithouse frequencies over time. The results of these analyses are stored in the R image [./R_images/comp.RData](./R_images/comp.RData).

# R Package Used

* trapezoid_2.0-0
* oxcAAR_1.0.0
* rcarbon_1.2.0
* dplyr_0.8.3
* readr_1.3.1
* magrittr_1.5
* TTR_0.23-5
* trapezoid_2.0-0

# Funding
This research was funded by the ERC grant _Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER)_ (Project N. 801953, PI: Enrico Crema). 

# Licence
CC-BY 3.0
