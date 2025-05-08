###### Simulation B ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Section 4.2 Simulation B: Relative Performance
## This code is used to generate the tables/figures in the paper. The results are saved in the "results" folder by running the script "02-simulation-B.R"
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
# choose one of them: the data are generated from lag0, lag0-1, lag0-7, or lag0-14
avglag <- 0
avglag <- 1
avglag <- 7
avglag <- 14


load(file = file.path(resultpath, paste0("02-simulation-", 
                                         as.character(avglag), 
                                         "-results-1.RData")))

load(file = file.path(resultpath, paste0("02-simulation-", 
                                         as.character(avglag), 
                                         "-results-2.RData")))

load(file = file.path(resultpath, paste0("02-simulation-", 
                                         as.character(avglag), 
                                         "-results-3.RData")))

load(file = file.path(resultpath, paste0("02-simulation-", 
                                         as.character(avglag), 
                                         "-results-4.RData")))



simresults <- c(simresults1, simresults2, simresults3, simresults4)

### TABLES for scenarios (i)---(iv)##########
datGAM_results <- lapply(simresults, function(results){

  out <- list(E = list(),
              other = list())

  if(length(results) > 1){
    Nt <- as.numeric(results$para[1])
    iter <- as.numeric(results$para[2])
    
    ## distributed lag term f(E)
    result_E <- results$datGAM$E
    true.eta_E <- result_E$true.eta
    
    
    out$E$DLNM.E.RMSE <- sqrt(mean((result_E$modDLNM$est - true.eta_E)^2))
    out$E$DLNM.E.Cvg <- mean((result_E$modDLNM$ul >= true.eta_E ) & (result_E$modDLNM$ll <= true.eta_E ))
    
    out$E$GAM01.E.RMSE <- sqrt(mean((result_E$modGAM01$est - true.eta_E)^2))
    out$E$GAM01.E.Cvg <- mean((result_E$modGAM01$ul >= true.eta_E ) & (result_E$modGAM01$ll <= true.eta_E ))
    
    out$E$GAM0.E.RMSE <- sqrt(mean((result_E$modGAM0$est - true.eta_E)^2))
    out$E$GAM0.E.Cvg <- mean((result_E$modGAM0$ul >= true.eta_E ) & (result_E$modGAM0$ll <= true.eta_E ))
    
    
    out$E$GAMavgmaxL.E.RMSE <- sqrt(mean((result_E$modGAMavgmaxL$est - true.eta_E)^2))
    out$E$GAMavgmaxL.E.Cvg <- mean((result_E$modGAMavgmaxL$ul >= true.eta_E ) & (result_E$modGAMavgmaxL$ll <= true.eta_E ))
    
    out$E$GAMavg7.E.RMSE <- sqrt(mean((result_E$modGAMavg7$est - true.eta_E)^2))
    out$E$GAMavg7.E.Cvg <- mean((result_E$modGAMavg7$ul >= true.eta_E ) & (result_E$modGAMavg7$ll <= true.eta_E ))
    
    out$E$pDLNM.E.RMSE <- sqrt(mean((result_E$modpDLNM$est - true.eta_E)^2))
    out$E$pDLNM.E.Cvg <- mean((result_E$modpDLNM$ul >= true.eta_E ) & (result_E$modpDLNM$ll <= true.eta_E ))
    
    
    out$E <- as.data.frame(out$E)
    out$E$Nt <- Nt
    out$E$iter <- iter
    
    ## other terms
    result_other <- results$datGAM$other
    true.eta_other <- result_other$true.eta
    
    out$other$DLNM.other.RMSE <- sqrt(mean((result_other$modDLNM$est - true.eta_other)^2))
    out$other$DLNM.other.Cvg <- mean((result_other$modDLNM$ul >= true.eta_other) & (result_other$modDLNM$ll <= true.eta_other))
    
    out$other$GAM01.other.RMSE <- sqrt(mean((result_other$modGAM01$est - true.eta_other)^2))
    out$other$GAM01.other.Cvg <- mean((result_other$modGAM01$ul >= true.eta_other) & (result_other$modGAM01$ll <= true.eta_other))
    
    
    out$other$GAM0.other.RMSE <- sqrt(mean((result_other$modGAM0$est - true.eta_other)^2))
    out$other$GAM0.other.Cvg <- mean((result_other$modGAM0$ul >= true.eta_other) & (result_other$modGAM0$ll <= true.eta_other))
    
    
    out$other$GAMavgmaxL.other.RMSE <- sqrt(mean((result_other$modGAMavgmaxL$est - true.eta_other)^2))
    out$other$GAMavgmaxL.other.Cvg <- mean((result_other$modGAMavgmaxL$ul >= true.eta_other) & (result_other$modGAMavgmaxL$ll <= true.eta_other))
    
    out$other$GAMavg7.other.RMSE <- sqrt(mean((result_other$modGAMavg7$est - true.eta_other)^2))
    out$other$GAMavg7.other.Cvg <- mean((result_other$modGAMavg7$ul >= true.eta_other) & (result_other$modGAMavg7$ll <= true.eta_other))
    
    
    out$other$pDLNM.other.RMSE <- sqrt(mean((result_other$modpDLNM$est - true.eta_other)^2))
    out$other$pDLNM.other.Cvg <- mean((result_other$modpDLNM$ul >= true.eta_other ) & (result_other$modpDLNM$ll <= true.eta_other ))
    
    
    out$other <- as.data.frame(out$other)
    out$other$Nt <- Nt
    out$other$iter <- iter
    
  }
  return(out)
}
)



datGAM_E <- lapply(datGAM_results, "[[", "E")
datGAM_E <- data.table::rbindlist(datGAM_E)

datGAM_other <- lapply(datGAM_results, "[[", "other")
datGAM_other <- data.table::rbindlist(datGAM_other)



### TABLES for scenario (v)##########
if(avglag == 0) {
  datDLNM_results <- lapply(simresults, function(results){
    
    out <- list(E = list(),
                other = list())
    
    if(length(results) > 1){
      Nt <- as.numeric(results$para[1])
      iter <- as.numeric(results$para[2])
      
      ## distributed lag term f(E)
      result_E <- results$datDLNM$E
      true.eta_E <- result_E$true.eta
      
      
      out$E$DLNM.E.RMSE <- sqrt(mean((result_E$modDLNM$est - true.eta_E)^2))
      out$E$DLNM.E.Cvg <- mean((result_E$modDLNM$ul >= true.eta_E) & (result_E$modDLNM$ll <= true.eta_E))
      
      out$E$GAM01.E.RMSE <- sqrt(mean((result_E$modGAM01$est - true.eta_E)^2))
      out$E$GAM01.E.Cvg <- mean((result_E$modGAM01$ul >= true.eta_E) & (result_E$modGAM01$ll <= true.eta_E))
      
      out$E$GAM0.E.RMSE <- sqrt(mean((result_E$modGAM0$est - true.eta_E)^2))
      out$E$GAM0.E.Cvg <- mean((result_E$modGAM0$ul >= true.eta_E) & (result_E$modGAM0$ll <= true.eta_E))
      
      
      out$E$GAMavgmaxL.E.RMSE <- sqrt(mean((result_E$modGAMavgmaxL$est - true.eta_E)^2))
      out$E$GAMavgmaxL.E.Cvg <- mean((result_E$modGAMavgmaxL$ul >= true.eta_E) & (result_E$modGAMavgmaxL$ll <= true.eta_E))
      
      out$E$GAMavg7.E.RMSE <- sqrt(mean((result_E$modGAMavg7$est - true.eta_E)^2))
      out$E$GAMavg7.E.Cvg <- mean((result_E$modGAMavg7$ul >= true.eta_E) & (result_E$modGAMavg7$ll <= true.eta_E))
      
      out$E$pDLNM.E.RMSE <- sqrt(mean((result_E$modpDLNM$est - true.eta_E)^2))
      out$E$pDLNM.E.Cvg <- mean((result_E$modpDLNM$ul >= true.eta_E) & (result_E$modpDLNM$ll <= true.eta_E))
      
      
      out$E <- as.data.frame(out$E)
      out$E$Nt <- Nt
      out$E$iter <- iter
      
      ## other
      result_other <- results$datDLNM$other
      true.eta_other <- result_other$true.eta
      
      out$other$DLNM.other.RMSE <- sqrt(mean((result_other$modDLNM$est - true.eta_other)^2))
      out$other$DLNM.other.Cvg <- mean((result_other$modDLNM$ul >= true.eta_other) & (result_other$modDLNM$ll <= true.eta_other))
      
  
      out$other$GAM01.other.RMSE <- sqrt(mean((result_other$modGAM01$est - true.eta_other)^2))
      out$other$GAM01.other.Cvg <- mean((result_other$modGAM01$ul >= true.eta_other) & (result_other$modGAM01$ll <= true.eta_other))
      
      out$other$GAM0.other.RMSE <- sqrt(mean((result_other$modGAM0$est - true.eta_other)^2))
      out$other$GAM0.other.Cvg <- mean((result_other$modGAM0$ul >= true.eta_other) & (result_other$modGAM0$ll <= true.eta_other))
      
      
      out$other$GAMavgmaxL.other.RMSE <- sqrt(mean((result_other$modGAMavgmaxL$est - true.eta_other)^2))
      out$other$GAMavgmaxL.other.Cvg <- mean((result_other$modGAMavgmaxL$ul >= true.eta_other) & (result_other$modGAMavgmaxL$ll <= true.eta_other))
      
      out$other$GAMavg7.other.RMSE <- sqrt(mean((result_other$modGAMavg7$est - true.eta_other)^2))
      out$other$GAMavg7.other.Cvg <- mean((result_other$modGAMavg7$ul >= true.eta_other) & (result_other$modGAMavg7$ll <= true.eta_other))
      
      out$other$pDLNM.other.RMSE <- sqrt(mean((result_other$modpDLNM$est - true.eta_other)^2))
      out$other$pDLNM.other.Cvg <- mean((result_other$modpDLNM$ul >= true.eta_other ) & (result_other$modpDLNM$ll <= true.eta_other ))
      
      
      
      out$other <- as.data.frame(out$other)
      out$other$Nt <- Nt
      out$other$iter <- iter
      
    }
    return(out)
  }
  )
  
  datDLNM_E <- lapply(datDLNM_results, "[[", "E")
  datDLNM_E <- data.table::rbindlist(datDLNM_E)
  
  datDLNM_other <- lapply(datDLNM_results, "[[", "other")
  datDLNM_other <- data.table::rbindlist(datDLNM_other)
  
  
}

### report tables ##########
n <- nrow(datGAM_E)
MCSE <- function(xx) sd(xx) / sqrt(n)

cat("distributed lag terms: \n")
knitr::kable(apply(filter(datGAM_E)%>%select(-iter), 2, mean), digits = 3)
knitr::kable(apply(filter(datGAM_E, Nt == 2000)%>%select(-iter), 2, MCSE), digits = 3)

cat("time trend: \n")
knitr::kable(apply(filter(datGAM_other, Nt == 2000)%>%select(-iter), 2, mean), digits = 3)
knitr::kable(apply(filter(datGAM_other, Nt == 2000)%>%select(-iter), 2, MCSE), digits = 3)



if(avglag == 0) {
  ## scenario (v)
  knitr::kable(apply(filter(datDLNM_E, Nt == 2000)%>%select(-iter), 2, mean), digits = 3)
  knitr::kable(apply(filter(datDLNM_E, Nt == 2000)%>%select(-iter), 2, MCSE), digits = 3)

  knitr::kable(apply(filter(datDLNM_other, Nt == 2000)%>%select(-iter), 2, mean), digits = 3)
  knitr::kable(apply(filter(datDLNM_other, Nt == 2000)%>%select(-iter), 2, MCSE), digits = 3)
}
