###### Supplementary Simulation: Misspecification Scenario B ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Web Appendix H.4.2 Supplementary Simulation: Misspecification Scenario B
## This code is used to generate the tables/figures in the paper. The results are saved in the "results" folder by running the script "06-simulation-appendix.R"
## Tianyi Pan
## 2025
###################################

## load packages
library(dplyr)
# devtools::install_github("tianyi-pan/aceDLNM")
library(aceDLNM)
library(ggplot2)

## ggplot2 settings
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25
GGPLOTSTSIZE <- 24
GGPLOTLEGENDSIZE <- 23
PLOTWIDTH <- 12
PLOTHEIGHT <- 9

## paths
resultpath <- "results"
figurepath <- "figures"
if (!dir.exists(figurepath)) dir.create(figurepath, recursive = TRUE)


## load data
load(file = file.path(resultpath, "06-simulation-appendix-results-1.RData"))
load(file = file.path(resultpath, "06-simulation-appendix-results-2.RData"))
load(file = file.path(resultpath, "06-simulation-appendix-results-3.RData"))
load(file = file.path(resultpath, "06-simulation-appendix-results-4.RData"))
simresults <- c(simresults1, simresults2, simresults3, simresults4)


simresults <- c(simresults1, simresults2, simresults3, simresults4)

### TABLES ##########
datC2_results <- lapply(simresults, function(results){

  out <- list(E = list())

  if(length(results) > 1){
    iter <- as.numeric(results$para[1])
    
    ## distributed lag terms
    result_E <- results$datC2$E
    true.eta_E <- result_E$true.eta
    
    
    out$E$DLNM.E.RMSE <- sqrt(mean((result_E$modDLNM$est - true.eta_E)^2))
    out$E$DLNM.E.Cvg <- mean((result_E$modDLNM$ul >= true.eta_E ) & (result_E$modDLNM$ll <= true.eta_E ))
    
    out$E <- as.data.frame(out$E)
    out$E$iter <- iter
    

    out$E$pDLNM.E.RMSE <- sqrt(mean((result_E$modpDLNM$est - true.eta_E)^2))
    out$E$pDLNM.E.Cvg <- mean((result_E$modpDLNM$ul >= true.eta_E ) & (result_E$modpDLNM$ll <= true.eta_E ))
    
  }
  return(out)
}
)


datC2_E <- lapply(datC2_results, "[[", "E")
datC2_E <- data.table::rbindlist(datC2_E)



## table 
datC2_E.summ <- group_by(datC2_E) %>% summarise(
  ACEDLNM.RMSE.mean = mean(DLNM.E.RMSE),
  ACEDLNM.RMSE.MCSE = sd(DLNM.E.RMSE) / sqrt(n()),
  ACEDLNM.Cvg.mean = mean(DLNM.E.Cvg),
  ACEDLNM.Cvg.MCSE = sd(DLNM.E.Cvg) / sqrt(n()),
  DRFDLNM.RMSE.mean = mean(pDLNM.E.RMSE),
  DRFDLNM.RMSE.MCSE = sd(pDLNM.E.RMSE) / sqrt(n()),
  DRFDLNM.Cvg.mean = mean(pDLNM.E.Cvg),
  DRFDLNM.Cvg.MCSE = sd(pDLNM.E.Cvg) / sqrt(n()),
  n = n()
)

knitr::kable(datC2_E.summ,  digits = 3)
