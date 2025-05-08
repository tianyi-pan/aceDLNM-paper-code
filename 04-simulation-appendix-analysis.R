###### Supplementary Simulation: Linear f ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Web Appendix H.3 Supplementary Simulation: Linear f
## This code is used to generate the tables/figures in the paper. The results are saved in the "results" folder by running the script "04-simulation-appendix.R"
## Tianyi Pan
## 2025
###################################

## load packages
library(dplyr)
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
load(file = file.path(resultpath, "04-simulation-appendix-results-1.RData"))
load(file = file.path(resultpath, "04-simulation-appendix-results-2.RData"))
load(file = file.path(resultpath, "04-simulation-appendix-results-3.RData"))
load(file = file.path(resultpath, "04-simulation-appendix-results-4.RData"))
load(file = file.path(resultpath, "04-simulation-appendix-results-5.RData"))
load(file = file.path(resultpath, "04-simulation-appendix-results-6.RData"))
load(file = file.path(resultpath, "04-simulation-appendix-results-7.RData"))
load(file = file.path(resultpath, "04-simulation-appendix-results-8.RData"))

simresults <- c(simresults1, simresults2, simresults3, simresults4,
                simresults5, simresults6, simresults7, simresults8)


### TABLES ##########
datDLNM_results <- lapply(simresults, function(results){
  
  out <- list(E = list(),
              other = list())
  
  if(length(results) > 1){
    Nt <- as.numeric(results$para[1])
    wltype <- results$para[2]
    fEtype <- results$para[3]
    iter <- as.numeric(results$para[4])
    
    ## distributed lag terms
    result_E <- results$datDLNM$E
    true.eta_E <- result_E$true.eta
    
    
    out$E$DLNM.E.RMSE <- sqrt(mean((result_E$modDLNM$est - true.eta_E)^2))
    out$E$DLNM.E.Cvg <- mean((result_E$modDLNM$ul >= true.eta_E) & (result_E$modDLNM$ll <= true.eta_E))
    out$E$DLNM.E.length <- mean((result_E$modDLNM$ul - result_E$modDLNM$ll))

    out$E$pDLNM.E.RMSE <- sqrt(mean((result_E$modpDLNM$est - true.eta_E)^2))
    out$E$pDLNM.E.Cvg <- mean((result_E$modpDLNM$ul >= true.eta_E) & (result_E$modpDLNM$ll <= true.eta_E))
    out$E$pDLNM.E.length <- mean((result_E$modpDLNM$ul - result_E$modpDLNM$ll))
    
    out$E <- as.data.frame(out$E)
    out$E$Nt <- Nt
    out$E$wltype <- wltype
    out$E$iter <- iter
    
    ## other terms
    result_other <- results$datDLNM$other
    true.eta_other <- result_other$datDLNM$true.eta
    
    out$other$DLNM.other.RMSE <- sqrt(mean((result_other$modDLNM$est - true.eta_other)^2))
    out$other$DLNM.other.Cvg <- mean((result_other$modDLNM$ul >= true.eta_other) & (result_other$modDLNM$ll <= true.eta_other))
    out$other$DLNM.other.length <- mean((result_other$modDLNM$ul - result_other$modDLNM$ll))

    out$other$pDLNM.other.RMSE <- sqrt(mean((result_other$modpDLNM$est - true.eta_other)^2))
    out$other$pDLNM.other.Cvg <- mean((result_other$modpDLNM$ul >= true.eta_other ) & (result_other$modpDLNM$ll <= true.eta_other ))
    out$other$pDLNM.other.length <- mean((result_other$modpDLNM$ul - result_other$modpDLNM$ll))
    
    
    out$other <- as.data.frame(out$other)
    out$other$Nt <- Nt
    out$other$wltype <- wltype
    out$other$fEtype <- fEtype
    out$other$iter <- iter
    
  }
  return(out)
}
)

datDLNM_E <- lapply(datDLNM_results, "[[", "E")
datDLNM_E <- data.table::rbindlist(datDLNM_E)

datDLNM_other <- lapply(datDLNM_results, "[[", "other")
datDLNM_other <- data.table::rbindlist(datDLNM_other)

  


## table 
E.summ <- group_by(datDLNM_E, wltype, Nt) %>% 
  summarise(
    ACEDLNM.RMSE.mean = mean(DLNM.E.RMSE),
    # ACEDLNM.RMSE.MCSE = sd(DLNM.E.RMSE) / sqrt(n()),
    ACEDLNM.Cvg.mean = mean(DLNM.E.Cvg),
    # ACEDLNM.Cvg.MCSE = sd(DLNM.E.Cvg) / sqrt(n()),
    ACEDLNM.length.mean = mean(DLNM.E.length),
    # ACEDLNM.length.MCSE = sd(DLNM.E.length) / sqrt(n()),
    DRFDLNM.RMSE.mean = mean(pDLNM.E.RMSE),
    # DRFDLNM.RMSE.MCSE = sd(pDLNM.E.RMSE) / sqrt(n()),
    DRFDLNM.Cvg.mean = mean(pDLNM.E.Cvg),
    # DRFDLNM.Cvg.MCSE = sd(pDLNM.E.Cvg) / sqrt(n()),
    DRFDLNM.length.mean = mean(pDLNM.E.length),
    # DRFDLNM.length.MCSE = sd(pDLNM.E.length) / sqrt(n()),
    n = n()
  )
knitr::kable(E.summ,  digits = 3)

other.summ <- group_by(datDLNM_other, wltype, Nt) %>% 
  summarise(
    ACEDLNM.RMSE.mean = mean(DLNM.other.RMSE),
    ACEDLNM.RMSE.MCSE = sd(DLNM.other.RMSE) / sqrt(n()),
    ACEDLNM.Cvg.mean = mean(DLNM.other.Cvg),
    ACEDLNM.Cvg.MCSE = sd(DLNM.other.Cvg) / sqrt(n()),
    DRFDLNM.RMSE.mean = mean(pDLNM.other.RMSE),
    DRFDLNM.RMSE.MCSE = sd(pDLNM.other.RMSE) / sqrt(n()),
    DRFDLNM.Cvg.mean = mean(pDLNM.other.Cvg),
    DRFDLNM.Cvg.MCSE = sd(pDLNM.other.Cvg) / sqrt(n()),
    n = n()
  )

E.summ %>% select(wltype, Nt, ACEDLNM.RMSE.mean, ACEDLNM.Cvg.mean, DRFDLNM.RMSE.mean, DRFDLNM.Cvg.mean) %>% 
  knitr::kable(digits = 3)
