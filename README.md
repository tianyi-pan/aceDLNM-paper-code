# Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines

Source code to replicate the results in the paper *Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines*

## Packages
### Install the `R` package `aceDLNM` from [GitHub](https://github.com/tianyi-pan/aceDLNM)
```R
# install.packages('devtools')
devtools::install_github("tianyi-pan/aceDLNM")
```
+ A brief [vignette](https://tianyi-pan.github.io/aceDLNM)

### Other required packages from `CRAN`
```R
install.packages(c(
  "mgcv",
  "parallel",
  "ggplot2",
  "dplyr",
  "Rcpp",
  "RcppEigen",
  "Matrix"
))
```

## Simulation (Section 4)
### Simulation A (Section 4.1)
+ Replicate the simulation: `01-simulation-A.R`
+ Analyse the results: `01-simulation-A-analysis.R`

### Simulation B (Section 4.2)
+ Replicate the simulation: `02-simulation-B.R`
+ Analyse the results: `02-simulation-B-analysis.R`

## Application (Section 5)
+ The datasets are private. 
+ R code for model fitting, plots and tables: `03-application.R`

---- 

## Supplementary Materials
### Web Appendix H.3
+ Replicate the simulation: `04-simulation-appendix.R`
+ Analyse the results: `04-simulation-appendix-analysis.R`
### Web Appendix H.4.1
+ Replicate the simulation: `05-simulation-appendix.R`
+ Analyse the results: `05-simulation-appendix-analysis.R`
### Web Appendix H.4.2
+ Replicate the simulation: `06-simulation-appendix.R`
+ Analyse the results: `06-simulation-appendix-analysis.R`