# Bayesian Approaches to Jomon Chronology and Demography: source code, data, and scripts
This repository contains the C14 Dates, R code and OxCal scripts associated with the manuscript "Bayesian Approaches to Jomon Chronology and Demography" by Crema, E.R., and Kobayashi, K.

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

