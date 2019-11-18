# Bayesian phase models, residential frequency data, and multi-proxy inferences of Jomon population dynamics: source code, data, and scripts

This repository contains the C14 Dates, R code and OxCal scripts which allows the full reproduction of all analyses in the paper "Bayesian phase models, residential frequency data, and multi-proxy inferences of Jomon population dynamics" by Crema, E.R., and Kobayashi, K.

The core analyses presented in the paper can be grouped in four stages: 1) data preparation; 2) Bayesian ceramic phase modelling; 3) Monte-Carlo simulation of pithouse frequency time-series; 4) summed probability distribution (SPD) of radiocarbon dates analyses; 5) comparison between pithouse data and SPD; and 6) visualisation. All analytical steps (aside for data preparation) are recorded in the [log.R](./log.R) file.

## Data Sets and Data Preparation

All raw data can be found in the [data](./data) directory. The file [c14dates.csv](./data/c14dates.csv) contains the radiocarbon dates used for the ceramic phase modelling, the folder [rekihaku14C](./data/rekihaku14C) contains CSV files of radiocarbon dates downloaded from the National Museum of Japanese History's [Database of radiocarbon dates published in Japanese archaeological research reports](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param). The file [bindC14csv.R](./data/rekihaku14C/bindC14csv.R) in the same directory contains details of the query used to download the data, as well as an R script for aggregating the files into an R image file (`westKantoNaganoC14.RData`) contained in the [R_images](./R_images) directory. Finally the sub-directory [./data/suzuki](./data/suzuki) contains CSV files of Jomon pithouse counts obtained from pages 88 to 93 of 「縄文時代集落の研究」("Research on Jomon Period Settlements"), by 鈴木 保彦  (Suzuki, Yasuhiko), published by 雄山閣 (Yuzankaku), Tokyo, in 2006.  <!-- Add processing and storage in pithouseData.RData -->

## Bayesian ceramic phase modelling

This stage of the analysis (recorded in the the [log.R](./log.R) file) involves the use of [OxCal](https://c14.arch.ox.ac.uk/oxcal/OxCal.html). The process is carried out in three steps. Firstly, potential outliers are detected within each set of dates associated with the same event (e.g. different organic residues from the same vessel). This is achievied by using the `outlierExcluder()` function (stored [here](./R/outlierAnalysis.R)) which internally calls OxCal from R using the [oxcAAR](https://CRAN.R-project.org/package=oxcAAR) package. The function eliminates potential outliers from our initial sets. The result of this routine is stored in the R image file [c14data.RData](./R_images/c14data.RData). 

The second step of analysis consisted in creating OxCal scripts for different probability distributions emulating putative _within-phase uncertainty_. This is achieved by using the `oxcalScriptGen()` function (located [in this file](./R/oxcalScriptCreator.R), which generates the OxCal scripts from the provided dates. The script (stored in the directory [./oxcal/oxcalscripts](./oxcal/oxcalscripts)) can then be loaded into OxCal for analyses. The results are stored in the directory [./oxcal/results](./oxcal/results), and include the .csv storing the posterior samples and a JavaScript file (.js) containing key statistics such as individual and overall agreement indices that can be read using the `oxcalReadjs()` function (see source code [here](./R/oxcalReadjs.R)). 

The third and final step consist of removing samples with low agreement index for each of the probability distributions and generate a new set of OxCal scripts (with suffix "R" to distinguish this second submission to OxCal, e.g. `gaussian.oxcal` and `gaussianR.oxcal`). The result of the OxCal analysis is then processed and the posterior samples organised into a series of objects stored in R image [posteriorSamples.RData](./R_images/posteriorSamples.RData).




# File Structure

```
./
├── data #Data Folder
│   ├── c14dates.csv #Radiocarbon dates for OxCal Model
│   ├── rekihaku14C #Radiocarbon dates for SPD analysis
│   │   ├── bindC14csv.R #Binder function for downloaded radiocarbon dates
│   │   ├── kanagawa_T_B_5_11_2019.csv #Radiocarbon dates from Kanagawa
│   │   ├── nagano_T_B_5_11_2019.csv #Radiocarbon dates from Nagano
│   │   ├── saitama_T_B_5_11_2019.csv #Radiocarbon dates from Saitama
│   │   ├── tokyo_T_B_5_11_2019.csv #Radiocarbon dates from Tokyo
│   │   └── yamanashi_T_B_5_11_2019.csv #Radiocarbon dates from Yamanashi
│   └── suzuki #Suzuki's Pithouse dataset
│       ├── conversion.csv #Baseline crossreference of ceramic phases 
│       ├── kanagawa.csv #Pithouse counts for Kanagwa
│       ├── nagano.csv #Pithouse counts for Nagano
│       ├── saitama.csv #Pithouse counts for Saitama
│       ├── tokyo.csv #Pithouse counts for Tokyo
│       └── yamanashi.csv #Pithouse counts for Yamanashi
├── log.R #Core Log file
├── manuscript #Manuscript Related Files
│   ├── figures
│   │   ├── figure1.pdf
│   │   ├── figure2.pdf
│   │   └── figurelog.R #R function for generating figures
│   └── tables 
│       └── table1_base.csv
├── oxcal
│   ├── oxcalscripts #OxCal submission scripts. *R.oxcal are resubmissions.
│   │   ├── gaussian.oxcal
│   │   ├── gaussianR.oxcal
│   │   ├── trapezoid.oxcal
│   │   ├── trapezoidR.oxcal
│   │   ├── uniform.oxcal
│   │   └── uniformR.oxcal
│   └── results #OxCal Output Files. mcmc*.csv are posterior samples.
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
│       └── uniform.txt
├── R #R scripts
│   ├── mcsim.R #Monte Carlo simulation routine
│   ├── outlierAnalysis.R #Functions for outlier analysis
│   ├── oxcalReadjs.R #Functions for reading oxcal *.js files
│   ├── oxcalScriptCreator.R #Functions for generating oxcal submission scripts
│   └── utilities.R #Variety of utility functions
├── README.md # This file
├── R_images #R image files
    ├── c14data.RData #C14 data for oxcal analysis
    ├── pithouseData.RData #Re-structured pithouse data
    ├── posteriorSamples.RData #Posterior samples from OxCal model
    ├── simdatesPithouses.RData #Simulated dates for Suzuki's pithouses
    ├── spdRes.RData #SPD analysis
    └── westKantoNaganoC14.RData #C14 data for SPD analysis
```

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

# Licence
