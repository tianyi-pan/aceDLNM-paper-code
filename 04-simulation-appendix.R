#!/usr/bin/env Rscript

###### Supplementary Simulation: Linear f ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Web Appendix H.3 Supplementary Simulation: Linear f
## This code is used to run the simulation study. To obtain the tables/figures in the paper, please see: 04-simulation-appendix-analysis.R
## Tianyi Pan
## 2025
###################################

## load packages
library(aceDLNM)
library(parallel)
library(mgcv)
library(dlnm); library(splines) ; library(tsModel) # packages for dlm (Gasparrini et.al. 2017)

options(mc.cores = parallel::detectCores()-1)

## paths
resultpath <- "results"
if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)


Nt.list <- c(1000, 2000)
wltype.list <- c("type1", "type2", "type3")
fEtype.list <- c("linear")
R <- 10000 # number of iterations; you could adjust this numbder to control the running time of this simulation

para <- expand.grid(Nt = Nt.list,
                    wltype = wltype.list,
                    fEtype = fEtype.list,
                    iter = 1:R)

para_run_list <- as.list(data.frame(t(para))) # convert matrix to list by row


## function to run the simulation
run_lst <- function(lst){
  tryCatch(
    expr = {
      set.seed(as.numeric(lst[4]))
      
      ## generate data
      cat("start: ", lst, "\n")
      
      gt <- function(x) 0.5*sin(x/150)
      
      maxL <- 14
      theta <- 8
      
      Nt <- as.numeric(lst[1])
      wltype <- lst[2]
      fEtype <- lst[3]
      
      ttmp <- 1:(maxL+Nt)
      ttmp <- ttmp[-(1:maxL)]
      meangt <- mean(gt(ttmp))
      otherterm <- data.frame(trend = gt(ttmp) - meangt,
                              intercept = 0)
      
      dat.gen <- GenerateData(fEtype = fEtype, wltype = wltype,
                              Nt = Nt,
                              theta = theta,
                              maxL = maxL,
                              other = otherterm,
                              interpolate = TRUE)
      
      eta_other.DLNM <- dat.gen$eta_other.sim
      eta_E.DLNM <- dat.gen$eta.sim - dat.gen$eta_other.sim
      
      dat <- data.frame(x = dat.gen$x,
                        t = dat.gen$t,
                        y = dat.gen$y)
      
      ### 1. ACE-DLNM
      if(Nt == 1000) {
        datDLNM.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                   smooth = ~s(t, bs = "bs", k = 10),
                                   dat = dat,
                                   kw = 20,
                                   kE = 20,
                                   maxL = maxL,
                                   eta = TRUE,
                                   verbose = FALSE)
      } else {
        datDLNM.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                   smooth = ~s(t, bs = "bs", k = 20),
                                   dat = dat,
                                   kw = 20,
                                   kE = 20,
                                   maxL = maxL,
                                   eta = TRUE,
                                   verbose = FALSE)
      }
      
      
      datDLNM.modDLNM.eta_E.est <- datDLNM.modDLNM$eta_E$est
      datDLNM.modDLNM.eta_E.ul <- datDLNM.modDLNM$eta_E$ul
      datDLNM.modDLNM.eta_E.ll <- datDLNM.modDLNM$eta_E$ll
      
      datDLNM.modDLNM.eta_other.est <- datDLNM.modDLNM$eta_other$est
      datDLNM.modDLNM.eta_other.ul <- datDLNM.modDLNM$eta_other$ul
      datDLNM.modDLNM.eta_other.ll <- datDLNM.modDLNM$eta_other$ll
      

      
      ### 2. DLM
      dat.DLM <- dat[-c(1:maxL),]
      Xt.deBoor <- function(tnew) as.numeric(aceDLNM:::deBoor(tnew, dat.gen$knots_x, dat.gen$alpha_x, 4))
      
      Q <- matrix(NA, nrow = Nt, ncol = 1000)
      for (i in 1:Nt) {
        Q[i,] <- sapply(seq(i+maxL+0.5, i-0.5, length.out = 1000), Xt.deBoor)
      }
      cb <- crossbasis(Q, lag=c(0,999), argvar=list("lin"),
                       arglag=list(fun='ps', df=19))
      cbgamPen <- cbPen(cb)
      
      
      if(Nt == 1000) {
        datDLNM.modepDLNM <- gam(dat.DLM$y~cb+s(dat.DLM$t, bs = "bs", k = 10),family=nb(),paraPen=list(cb=cbgamPen),method='REML')  
      } else {
        datDLNM.modepDLNM <- gam(dat.DLM$y~cb+s(dat.DLM$t, bs = "bs", k = 20),family=nb(),paraPen=list(cb=cbgamPen),method='REML')
      }
      
      
      
      pred.lpmatrix <- predict(datDLNM.modepDLNM, type = "lpmatrix")
      datDLNM.modepDLNM.eta_E.est <- as.vector(pred.lpmatrix[, 1:20] %*% datDLNM.modepDLNM$coefficients[1:20])
      datDLNM.modepDLNM.eta_other.est <- as.vector(pred.lpmatrix[, -(1:20)] %*% datDLNM.modepDLNM$coefficients[-(1:20)])
      
      datDLNM.modepDLNM.eta_E.sd <- sqrt(diag(pred.lpmatrix[, 1:20] %*% vcov(datDLNM.modepDLNM)[1:20, 1:20] %*% t(pred.lpmatrix[, 1:20])))
      datDLNM.modepDLNM.eta_other.sd <- sqrt(diag(pred.lpmatrix[, -(1:20)] %*% vcov(datDLNM.modepDLNM)[-(1:20), -(1:20)] %*% t(pred.lpmatrix[, -(1:20)])))
      
      datDLNM.modepDLNM.eta_E.ul <- as.vector(datDLNM.modepDLNM.eta_E.est + 1.96*datDLNM.modepDLNM.eta_E.sd)
      datDLNM.modepDLNM.eta_E.ll <- as.vector(datDLNM.modepDLNM.eta_E.est - 1.96*datDLNM.modepDLNM.eta_E.sd)
      
      datDLNM.modepDLNM.eta_other.ul <- as.vector(datDLNM.modepDLNM.eta_other.est + 1.96*datDLNM.modepDLNM.eta_other.sd)
      datDLNM.modepDLNM.eta_other.ll <- as.vector(datDLNM.modepDLNM.eta_other.est - 1.96*datDLNM.modepDLNM.eta_other.sd)
      
      cat("finish: ", lst, "\n")
      
      #### RETURN ###########
      out = list(para = c(lst), 
                 datDLNM = list(E = list(true.eta = eta_E.DLNM,
                                  modDLNM = list(est = datDLNM.modDLNM.eta_E.est,
                                                 ll = datDLNM.modDLNM.eta_E.ll,
                                                 ul = datDLNM.modDLNM.eta_E.ul),
                                  modpDLNM = list(est = datDLNM.modepDLNM.eta_E.est,
                                                  ll = datDLNM.modepDLNM.eta_E.ll,
                                                  ul = datDLNM.modepDLNM.eta_E.ul)),
                         other = list(true.eta = eta_other.DLNM,
                                      modDLNM = list(est = datDLNM.modDLNM.eta_other.est,
                                                     ll = datDLNM.modDLNM.eta_other.ll,
                                                     ul = datDLNM.modDLNM.eta_other.ul),
                                      modpDLNM = list(est = datDLNM.modepDLNM.eta_other.est,
                                                      ll = datDLNM.modepDLNM.eta_other.ll,
                                                      ul = datDLNM.modepDLNM.eta_other.ul))
                         )
                 )
      return(out)
    },
    error = function(e){
      cat("error in ", lst, ".\n")
      warning(e)
    })
}



## run the simulation
len <- length(para_run_list)
para_run_list1 <- para_run_list[0+(1:(len/8))]
para_run_list2 <- para_run_list[(1*len/8)+(1:(len/8))]
para_run_list3 <- para_run_list[(2*len/8)+(1:(len/8))]
para_run_list4 <- para_run_list[(3*len/8)+(1:(len/8))]
para_run_list5 <- para_run_list[(4*len/8)+(1:(len/8))]
para_run_list6 <- para_run_list[(5*len/8)+(1:(len/8))]
para_run_list7 <- para_run_list[(6*len/8)+(1:(len/8))]
para_run_list8 <- para_run_list[(7*len/8)+(1:(len/8))]



simresults1 <- mclapply(para_run_list1, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-1.RData"))
rm(simresults1)

simresults2 <- mclapply(para_run_list2, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-2.RData"))
rm(simresults2)

simresults3 <- mclapply(para_run_list3, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-3.RData"))
rm(simresults3)

simresults4 <- mclapply(para_run_list4, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-4.RData"))
rm(simresults4)


simresults5 <- mclapply(para_run_list5, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-5.RData"))
rm(simresults5)


simresults6 <- mclapply(para_run_list6, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-6.RData"))
rm(simresults6)


simresults7 <- mclapply(para_run_list7, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-7.RData"))
rm(simresults7)


simresults8 <- mclapply(para_run_list8, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "04-simulation-appendix-results-8.RData"))
rm(simresults8)

cat("DONE!")