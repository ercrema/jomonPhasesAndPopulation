# jomonPhasesAndPopulation
This repository contains the C14 Dates, R scripts, and LaTeX source code associated with the manuscript "Combining Bayesian Chronological Models with the Aoristic/Monte-Carlo approach: a case study using a revised Chronology of Jomon Pottery Phases in Japan" by Crema, E.R., and Kobayashi, K.

# File Structure

```
│   README.md ... this document
│   esm.Rmd ... supplementary material with technical details and workflow.
│
|
└───manuscript
│   │   manuscript.tex ... submitted manuscript
│   │   bibfile.bib ... bibTeX file
|   |   model2-names.bst ... LaTeX style file for elsevier
|   |   
│   │
│   └───figures
│       │   limitation_aoristic_mc.pdf
│       │   posterior_jomon_phases.pdf
│       │   pithousecounts_and_roc.pdf
|       |   comparison_crema_imamura.pdf
|       |   figureScript.R ... R script for generating figures
│   
└───data
|   │   c14dates.csv ... C14 dates used for the Bayesian Chronological Models
|   │   kobayashi2008.csv ... Start and End date of Jomon Phases used in Crema 2012
|   |   imamura1996.csv ... Pithouse counts extracted from Imamura 1996 published graph.
|   |   
|   |
|   └───suzuki1986 ... pithouse counts with time-span of existence
|        |   nagano.csv 
|        |   saitama.csv
|        |   tokyo.csv  
|        |   yamanashi.csv
|        |   kanagawa.csv   
|        |   conversion.csv
|
|
└───R
|   |   oxcalScriptGenerator.R ... R script for generating oxcal Scripts
|   |   oxcalOutputReader.R ... Functions for reading oxcal .js files
|   |   doubleMonteCarlo.R ... R script for generating random dates of pithouse constructions
|   |   utilities.R ... several utitlity functions.
|
|
└───oxcal
|   |   oxcalOutput.RData ... Aggregated R image file of oxcal output
|   | 
|   |
|   └───submissionScripts ... oxcal submission scripts
|   |    |   gaussian.oxcal
|   |    |   uniform.oxcal
|   |    |   trapezoid.oxcal
|   | 
|   |
|   └───output ... output data from OxCal
|        |   mcmcGaussian.csv
|        |   mcmcUniform.csv  
|        |   mcmcTrapezoid.csv  
|        |   gaussian.js
|        |   uniform.js  
|        |   trapezoid.js 
|        |   gaussian.txt
|        |   uniform.txt  
|        |   trapezoid.txt  
|        |   gaussian.log
|        |   uniform.log  
|        |   trapezoid.log  
```



# R Package Used

* trapezoid_2.0-0
* oxcAAR_1.0.0
* rcarbon_1.2.0
* dplyr_0.8.3
* readr_1.3.1
* magrittr_1.5


#
