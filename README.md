[![DOI](https://zenodo.org/badge/215105753.svg)](https://zenodo.org/badge/latestdoi/215105753)

# A multi-proxy inference of Jōmon population dynamics using Bayesian phase models, residential data, and summed probability distribution of <sup>14</sup>C dates: source code, data, and scripts

This repository contains an updated version of the data and scripts used in the following paper.

Crema, E.R., Kobayashi, K., 2020. A multi-proxy inference of Jōmon population dynamics using bayesian phase models, residential data, and summed probability distribution of <sup>14</sup>C dates. Journal of Archaeological Science 117, 105136. DOI: https://doi.org/10.1016/j.jas.2020.105136

The original repository contained scripts based on the IntCal13 and Marine13 calibration curves and can be accessed [here](https://github.com/ercrema/jomonPhasesAndPopulation/tree/v.2.0). Analyses and results contained in this repository are based on the IntCal20 and Marine20 curves and contains. This version contains also an updated version of the function `mcsim()` as well as a new rmarkdown file with a [short tutorial](./calibrating_jomon_ceramic_phases.html) on how to use the function with different datasets.

The main workflow is recorded in the [log.R](./log.R) file and outputs are stored as R image files located in the [R_images](./R_images) directory.   

## Data Sets and Data Preparation

All raw data used in the paper can be found in the [data](./data) directory. The file [c14dates.csv](./data/c14dates.csv) contains the radiocarbon dates used for the ceramic phase modelling, the folder [rekihaku14C](./data/rekihaku14C) contains CSV files of radiocarbon dates downloaded from the National Museum of Japanese History's [Database of radiocarbon dates published in Japanese archaeological research reports](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param). The R script file [bindC14csv.R](./data/rekihaku14C/bindC14csv.R) in the same directory contains details of the query used to download the data, as well as a script for aggregating the files into an R image file (`spdC14.RData`) contained in the [R_images](./R_images) directory. The sub-directory [./data/suzuki](./data/suzuki) contains CSV files of Jomon pithouse counts obtained from tables on pages 88 to 93 of 「縄文時代集落の研究」("Research on Jomon Period Settlements"), by 鈴木 保彦  (Suzuki, Yasuhiko), published by 雄山閣 (Yuzankaku), Tokyo, in 2006. The data has been digitised into separate CSV files and aggregated into a long format data.frame containing a pithouse per row. This process is recorded in the R script [pithouseBinder.R](./data/suzuki/pithouseBinder.R) and the resulting data.frame is stored in the R image file (`pithouseData.RData`).

## Bayesian ceramic phase modelling

Bayesian modelling have been carried out using [OxCal v4.4](https://c14.arch.ox.ac.uk/oxcal/OxCal.html), with the preparation of OxCal scripts done in R. The analyses was conducted in three stages. Firstly, potential outliers were detected within sets of dates associated with the same event (e.g. different organic residues from the same vessel). This was achievied by using the `outlierExcluder()` function (stored [here](./R/outlierAnalysis.R)) which internally calls OxCal from R using the [oxcAAR](https://CRAN.R-project.org/package=oxcAAR) package. The function eliminates potential outliers from the initial set of radiocarbon dates and result of this routine was stored in the R image file [c14data.RData](./R_images/c14data.RData). 

The second step of analysis consisted in creating OxCal scripts for different probability distributions emulating putative _within-phase uncertainty_. This was achieved by using the `oxcalScriptGen()` function (located [in this file](./R/oxcalScriptCreator.R). The output scripts (stored in the directory [./oxcal/oxcalscripts](./oxcal/oxcalscripts)) were then loaded into OxCal for analyses. The results of the Bayesian analyses are stored in the directory [./oxcal/results](./oxcal/results), and include the .csv storing the posterior samples and a JavaScript file (.js) containing key statistics such as individual and overall agreement indices, read in R using the `oxcalReadjs()` function (see source code [here](./R/oxcalReadjs.R)). 

The third and final step consisted of removing samples with low agreement index and generate a new set of OxCal scripts (with suffix "R" to distinguish this second submission to OxCal, e.g. `gaussian.oxcal` and `gaussianR.oxcal`). The result of the OxCal analysis was then processed and the posterior samples organised into a series of objects stored in the R image [posteriorSamples.RData](./R_images/posteriorSamples.RData).

## Monte-Carlo Simulation of Pithouse Dates

Dates of individual pithouse were simulated 5,000 times taking account: 1) the uncertatinty within the phase ( _within phase uncertainty_ ); 2) the uncertainty in defining the membership of the pithouse to a particular phase ( _phase assignement uncertainty_ ); and the uncertainty associated with the parameters of the probability distribution describing the _within phase uncertainty_  (i.e. the _phase boundary uncertainty_). The function `mcsim()` (see source code [here](./R/mcsim.R)) was used for this purpose. The R image fle [simdatesPithouses.RData](./R/simdatesPithouses.RData) contains the 5,000 set of simulated dates, along with counts organised in 100 years bins (between 8,000 and 3,000 cal BP), and the outcome of a composite kernel density estimate analyses (see details in the [log.R](./log.R) file).   

## SPD Analysis

Summed probability distribution (SPD) of radiocarbon dates have been generated using [rcarbon](https://CRAN.R-project.org/package=rcarbon) <!-- specify whether this 1.2 or 1.3 depending on the radiocarbon paper -->. To enable correlation analyses with the residential data, a matrix of 5000 sets of randomly sampled calendar dates was created and aggregated by the same 100 years intervals between 8,000 and 3,000 cal BP. All SPD analysis related R objects are stored in the R image [./R_images/spdRes.RData](./R_images/spdRes.RData).  


## Comparisons between pithouse data and radiocarbon density

The time-series of residential and radiocarbon density have been compared via correlation analyses and the `modelTest()` function in [rcarbon](https://CRAN.R-project.org/package=rcarbon). The former was carried out by generating 5,000 correlation values by iteratively comparing the the time-series of simulated pithouse dates and randomly sampled calendar dates from the calibrated radiocarbon dates. The latter compared the observed annual growth rate in the SPD against an expectation derived from the average trend obtained from the composite kernel density estimate of pithouse frequencies over time. The results of these analyses are stored in the R image [./R_images/comp.RData](./R_images/comp.RData).

# File Structure
```
.
├── data
│   ├── c14dates.csv
│   ├── rekihaku14C
│   │   ├── bindC14csv.R
│   │   ├── kanagawa_M_B_11_3_2020.csv
│   │   ├── kanagawa_T_B_5_11_2019.csv
│   │   ├── nagano_T_B_5_11_2019.csv
│   │   ├── saitama_M_B_11_3_2020.csv
│   │   ├── saitama_T_B_5_11_2019.csv
│   │   ├── tokyo_M_B_11_3_2020.csv
│   │   ├── tokyo_M_B_11_3_2020.csv.csv
│   │   ├── tokyo_T_B_5_11_2019.csv
│   │   └── yamanashi_T_B_5_11_2019.csv
│   └── suzuki
│       ├── kanagawa.csv
│       ├── nagano.csv
│       ├── pithouseBinder.R
│       ├── saitama.csv
│       ├── tokyo.csv
│       └── yamanashi.csv
├── esm.pdf
├── esm.Rmd
├── log.R
├── manuscript
│   ├── figures
│   │   ├── figure1.pdf
│   │   ├── figure2.pdf
│   │   ├── figure3.pdf
│   │   ├── figure4.pdf
│   │   └── figurelog.R
│   └── tables
│       └── table1_base.csv
├── oxcal
│   ├── oxcalscripts
│   │   ├── gaussian.oxcal
│   │   ├── gaussianR.oxcal
│   │   ├── trapezoid.oxcal
│   │   ├── trapezoidR.oxcal
│   │   ├── uniform.oxcal
│   │   └── uniformR.oxcal
│   └── results
│       ├── gaussian.js
│       ├── gaussian.log
│       ├── gaussian.oxcal
│       ├── gaussianR.js
│       ├── gaussianR.log
│       ├── gaussianR.oxcal
│       ├── gaussianR.txt
│       ├── gaussian.txt
│       ├── mcmcGaussian.csv
│       ├── mcmcGaussianR.csv
│       ├── mcmcTrapezoid.csv
│       ├── mcmcTrapezoidR.csv
│       ├── mcmcUniform.csv
│       ├── mcmcUniformR.csv
│       ├── trapezoid.js
│       ├── trapezoid.log
│       ├── trapezoid.oxcal
│       ├── trapezoidR.js
│       ├── trapezoidR.log
│       ├── trapezoidR.oxcal
│       ├── trapezoidR.txt
│       ├── trapezoid.txt
│       ├── uniform.js
│       ├── uniform.log
│       ├── uniform.oxcal
│       ├── uniformR.js
│       ├── uniformR.log
│       ├── uniformR.oxcal
│       ├── uniformR.txt
│       └── uniform.txt
├── R
│   ├── mcsim.R
│   ├── outlierAnalysis.R
│   ├── oxcalReadjs.R
│   ├── oxcalScriptCreator.R
│   └── utilities.R
├── README.md
└── R_images
    ├── c14data.RData
    ├── comp.RData
    ├── pithouseData.RData
    ├── posteriorSamples.RData
    ├── simdatesPithouses.RData
    ├── spdC14.RData
    └── spdRes.RData
```

# R Settings
```
attached base packages:
[1] stats     graphics  grDevices utils     methods   base     

other attached packages:
[1] readr_1.3.1     trapezoid_2.0-0 TTR_0.23-5      rcarbon_1.4.1  
[5] oxcAAR_1.0.0    dplyr_0.8.3     magrittr_1.5    nvimcom_0.9-82 

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2            pillar_1.4.2         
 [3] compiler_3.6.1        xts_0.11-2           
 [5] iterators_1.0.12      tools_3.6.1          
 [7] rpart_4.1-15          goftest_1.1-1        
 [9] jsonlite_1.6          tibble_2.1.3         
[11] nlme_3.1-140          lattice_0.20-38      
[13] mgcv_1.8-28           pkgconfig_2.0.3      
[15] rlang_0.4.0           Matrix_1.2-17        
[17] foreach_1.4.7         curl_4.2             
[19] parallel_3.6.1        spatstat.data_1.4-0  
[21] stringr_1.4.0         hms_0.4.2            
[23] spatstat.utils_1.13-0 grid_3.6.1           
[25] tidyselect_0.2.5      glue_1.3.1           
[27] R6_2.4.0              sp_1.3-2             
[29] polyclip_1.10-0       purrr_0.3.2          
[31] deldir_0.1-23         tensor_1.5           
[33] splines_3.6.1         codetools_0.2-16     
[35] assertthat_0.2.1      abind_1.4-5          
[37] spatstat_1.61-0       stringi_1.4.3        
[39] doParallel_1.0.15     crayon_1.3.4         
[41] zoo_1.8-6   
```

# Funding
This research was funded by the ERC grant _Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER)_ (Project N. 801953, PI: Enrico Crema). 

# Licence
CC-BY 3.0
