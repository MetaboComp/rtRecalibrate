# rtRecalibrate

**Recalibration of retention time axis in mzML objects (e.g. from LC-HRMS)**  
Associate Professor Carl Brunius  <carl.brunius@chalmers.se>  
Department of Life Sciences
Chalmers University of Technology www.chalmers.se

## General description
The package contains functions to perform alignment of the  
retention time (RT) axis based on metabolite feature landmarks.  
NB: The RT correction requires that metabolite feature landmarks  
have been previously identified (see the Landmark package and tutorial  
at github.com/MetaboComp).

## Installation
- Please ensure that R is installed (https://www.r-project.org/)
- I also recommend installing RStudio for a smoother R experience (https://rstudio.com/) 

mzRecalibrate depends on some packages from BioConductor
```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("mzR", quietly = TRUE)) BiocManager::install("mzR")
if (!require("MSnbase", quietly = TRUE)) BiocManager::install("MSnbase")
if (!require("xcms", quietly = TRUE)) BiocManager::install("xcms")
```

To install the mzRecalibrate package you also need the `remotes` R package
```
if (!require("remotes", quietly = TRUE)) install.packages('remotes')
remotes::install_gitlab('CarlBrunius/rtRecalibrate')
```

## Version history
version | date | comment
:------ | :--- | :------
0.1.03  | 2023-11-30 | Added plotRTAdjust()
0.1.02  | 2023-11-20 | Bug fix addressing lamas in tibble format
0.1.01  | 2023-11-17 | Refactoring by Jo Rainer
0.1.00  | 2023-11-10 | Hello world. Transferred function files from early development.
