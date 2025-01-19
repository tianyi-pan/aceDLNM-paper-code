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

## Simulation 
### Simulation A: 
+ Replicate the simulation: `01-simulation-1.R`
+ Analyze the results: `01-simulation-1-analysis.R`

### Simulation B: 
+ Replicate the simulation: `02-simulation-2.R`
+ Analyze the results: `02-simulation-2-analysis.R`

## Application 
+ The datasets are private. 
+ R code for model fitting, plots and tables: `03-application-respiratory.R`,`03-application-circulatory.R`
