#!/usr/bin/env Rscript

###### Supplementary Simulation: Misspecification Scenario B ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Web Appendix H.4.2 Supplementary Simulation: Misspecification Scenario B
## This code is used to run the simulation study. To obtain the tables/figures in the paper, please see: 06-simulation-appendix-analysis.R
## Tianyi Pan
## 2025
###################################

## load packages
# devtools::install_github("tianyi-pan/aceDLNM")
library(aceDLNM)
library(parallel)
library(mgcv)

library(dlnm); library(splines) ; library(tsModel) # packages for dlnm (Gasparrini et.al. 2017)

options(mc.cores = parallel::detectCores()-1)

## paths
resultpath <- "results"
if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)


R <- 10000  # number of iterations; you could adjust this numbder to control the running time of this simulation
para <- expand.grid(iter = 1:R)
para_run_list <- as.list(data.frame(t(para))) # convert matrix to list by row



##### generate data from DRF-DLNM ########
### Code from https://onlinelibrary.wiley.com/doi/10.1111/biom.12645
x <- chicagoNMMAPS$temp
x <- (x-min(x))/diff(range(x))*10
Q <- Lag(x,0:40)
# BASIC FUNCTIONS TO SIMULATE UNIDIMENSIONAL SHAPES
flin <- function(x) 0.1*x
fflex <- function(x) {
  coef <- c(0.2118881,0.1406585,-0.0982663,0.0153671,-0.0006265)
  as.numeric(outer(x,0:4,'^')%*%coef)
}
fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))
wconst <- function(lag) lag-lag+0.20
wdecay <- function(lag) exp(-lag/2)
wpeak1 <- function(lag) 12*dnorm(lag,8,5)
wpeak2 <- function(lag) 15*dnorm(lag,8,10)
wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))
# FUNCTIONS TO SIMULATE THE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE
fplane <- function(x,lag) 0.1 * (flin(x)-flin(2)) * wconst(lag)
ftemp <- function(x,lag) 0.1 * (fflex(x)-fflex(5)) * 
  ifelse(is.na(x),NA,ifelse(x>=5,wdecay(lag),wpeak1(lag)))
fcomplex <- function(x,lag) 0.1 * (fdnorm(x)-fdnorm(5)) * 
  ifelse(is.na(x),NA,ifelse(x>=5,wdnorm(lag),wpeak2(lag)))
# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim <- c("fplane","ftemp","fcomplex")
names(combsim) <- c("Plane","Temperature","Complex")
cen <- c(2,5,5)
# LIST WITH TRUE EFFECT SURFACES OF FOR EACH COMBINATIONS
trueeff <- lapply(combsim, function(fun) {
  temp <- outer(seq(0,10,0.25),0:40,fun)
  dimnames(temp) <- list(seq(0,10,0.25),paste("lag",0:40,sep=""))
  return(temp)
})
names(trueeff) <- names(combsim)
# FUNCTION TO COMPUTE THE CUMULATIVE EFFECT GIVEN AN EXPOSURE HISTORY
fcumeff <- function(hist,lag,fun) sum(do.call(fun,list(hist,lag)))
# BASELINE
base <- c(15,150,15)
# NOMINAL VALUE
qn <- qnorm(0.975)
j <- 2 ## for ftemp
# combsim[2]: ftemp
# SIMULATE THE DATA
cumeff <- apply(Q,1,fcumeff,0:40,combsim[j])


######### run simulation #########
run_lst <- function(lst){
  tryCatch(
    expr = {
      
      i <- as.numeric(lst[1])
      
      ## random seed from: https://onlinelibrary.wiley.com/doi/10.1111/biom.12645
      seed <- 13041975 + i 
      set.seed(seed)

      ## generate data
      cat("start: ", lst, "\n")
      
      eta.C2 <- log(base[j])+cumeff
      suppressWarnings(y <- rpois(length(x),exp(eta.C2)))
      
      eta_E.C2 <- eta.C2[-(1:40)]
      ### prepare data
      dat.C2 <- data.frame(
        t = 1:length(x),
        x = x,
        y = y
      )
      
      ##### MODEL FITTING #######
      ##### 1. C2 ##########
      # true: eta.C2
      ## 1.1 ACE-DLNM
    
      datC2.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                   dat = dat.C2,
                                   pc = NULL,
                                   kw = 20, 
                                   kE = 20,
                                   maxL = 40,
                                   eta = TRUE,
                                   verbose = FALSE)
    
      
      datC2.modDLNM.summ <- summary(datC2.modDLNM, plot = FALSE)
      
      datC2.modDLNM.eta_E.est <- datC2.modDLNM$eta_E$est
      datC2.modDLNM.eta_E.ul <- datC2.modDLNM$eta_E$ul
      datC2.modDLNM.eta_E.ll <- datC2.modDLNM$eta_E$ll
      
      datC2.modDLNM.eta_other.est <- datC2.modDLNM$eta_other$est
      datC2.modDLNM.eta_other.ul <- datC2.modDLNM$eta_other$ul
      datC2.modDLNM.eta_other.ll <- datC2.modDLNM$eta_other$ll
      
      ## 1.2 DRF-DLNM
      cb <- crossbasis(Q[,1],lag=c(0,40),argvar=list(fun='ps',df=9),
                       arglag=list(fun='ps'))
      
      cbgamPen <- cbPen(cb)
      datC2.modepDLNM <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
      
      
      pred.lpmatrix <- predict(datC2.modepDLNM, type = "lpmatrix")
      
      datC2.modepDLNM.eta_E.est <- as.vector(pred.lpmatrix %*% datC2.modepDLNM$coefficients)
      
      datC2.modepDLNM.eta_E.sd <- sqrt(diag(pred.lpmatrix %*% vcov(datC2.modepDLNM) %*% t(pred.lpmatrix)))
      
      
      datC2.modepDLNM.eta_E.ul <- as.vector(datC2.modepDLNM.eta_E.est + 1.96*datC2.modepDLNM.eta_E.sd)
      datC2.modepDLNM.eta_E.ll <- as.vector(datC2.modepDLNM.eta_E.est - 1.96*datC2.modepDLNM.eta_E.sd)
      

      cat("finish: ", lst, "\n")
      #### RETURN ###########
      out = list(para = c(lst),
                  datC2 = list(E = list(true.eta = eta_E.C2,
                                        modDLNM = list(funs = datC2.modDLNM.summ$est,
                                                       est = datC2.modDLNM.eta_E.est,
                                                       ll = datC2.modDLNM.eta_E.ll,
                                                       ul = datC2.modDLNM.eta_E.ul),
                                        modpDLNM = list(est = datC2.modepDLNM.eta_E.est,
                                                        ll = datC2.modepDLNM.eta_E.ll,
                                                        ul = datC2.modepDLNM.eta_E.ul))
                               ))
      return(out)
    },
    error = function(e){
      cat("error in ", lst, ".\n")
      warning(e)
    }
  )
}

  
  
## run simulation
len <- length(para_run_list)
para_run_list1 <- para_run_list[0+(1:(len/4))]
para_run_list2 <- para_run_list[(1*len/4)+(1:(len/4))]
para_run_list3 <- para_run_list[(2*len/4)+(1:(len/4))]
para_run_list4 <- para_run_list[(3*len/4)+(1:(len/4))]



simresults1 <- mclapply(para_run_list1, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "06-simulation-appendix-results-1.RData"))
rm(simresults1)

simresults2 <- mclapply(para_run_list2, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "06-simulation-appendix-results-2.RData"))
rm(simresults2)

simresults3 <- mclapply(para_run_list3, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "06-simulation-appendix-results-3.RData"))
rm(simresults3)

simresults4 <- mclapply(para_run_list4, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "06-simulation-appendix-results-4.RData"))
rm(simresults4)


cat("DONE!")
