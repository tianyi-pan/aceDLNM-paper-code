#!/usr/bin/env Rscript

###### Supplementary Simulation: Misspecification Scenario A ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Web Appendix H.4.1 Supplementary Simulation: Misspecification Scenario A
## This code is used to run the simulation study. To obtain the tables/figures in the paper, please see: 05-simulation-appendix-analysis.R
## Tianyi Pan
## 2025
###################################

## load packages
library(aceDLNM)
library(parallel)
library(mgcv)

library(dlnm); library(splines) ; library(tsModel) # packages for DRF-DLNM (Gasparrini et.al. 2017)


options(mc.cores = parallel::detectCores()-1)

## paths
resultpath <- "results"
if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)

R <- 10000 # number of iterations; you could adjust this numbder to control the running time of this simulation


a.list <- sqrt(seq(0,1, by = 0.1)) # relative roughness: control the shape of the association function. If the association if linear, a = 0. 
para <- expand.grid(Nt = c(2000), 
                    wltype = c("type1"), 
                    fxtype = a.list, 
                    iter = 1:R)
para_run_list <- as.list(data.frame(t(para))) # convert matrix to list by row


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
      fxtype <- lst[3]
      a <- as.numeric(fxtype)
      
      ttmp <- 1:(maxL+Nt)
      ttmp <- ttmp[-(1:maxL)]
      meangt <- mean(gt(ttmp))
      otherterm <- data.frame(trend = gt(ttmp) - meangt,
                              intercept = 0)
      
      gt <- function(x) 0.5*sin(x/150)
      theta <- 8
      maxL <- 14
      
      
      ## obtain x
      t <- 1:(Nt + maxL)
      x <- PM25.waterloo$PM25[t]
      
      x.min <- min(x); x.max <- max(x)

        
      
      switch (wltype,
              type1 = {
                wlc <- function(l) dnorm(l, mean = 3, sd = 3.5)^2
                wl_de <- sqrt(integrate(wlc, lower = 0, upper = 15)$value)
                wl <- function(l) dnorm(l, mean = 3, sd = 3.5)/wl_de
              },
              type2 = {
                wlc <- function(l) (1/(1+exp(l-8)))^2
                wl_de <- sqrt(integrate(wlc, lower = 0, upper = 15)$value)
                wl <- function(l) 1/(1+exp(l-8))/wl_de
              },
              type3 = {
                wlc <- function(l) (1/(1+exp(l-1)))^2
                wl_de <- sqrt(integrate(wlc, lower = 0, upper = 15)$value)
                wl <- function(l) 1/(1+exp(l-1))/wl_de
              }
      )
      
      wl.discrete <- sapply(0:maxL, wl)

      
      
      ## association function. parameter "a" controls the roughness of the function
      x.med <- (x.min + x.max)/2
      fx <- function(x) ((x-(x.max+x.min)/2) / ((x.max - x.min)/3.5) + 3)/2 - a*sqrt(0.5)*(x-x.med)^2/((x.max-x.med)^2) + a*sqrt(0.5)
      

      eta_E.C1 <- sapply( (1:Nt)+14, function(s) 
          sum(fx(x[s-c(0:14)]) * wl.discrete)
        )
      eta_other.C1 <- gt(ttmp) - meangt

      eta.C1 <- eta_E.C1 + gt(ttmp) - meangt
      y.C1 <- sapply(eta.C1, function(eta.) rnbinom(n = 1, size = theta, mu = exp(eta.)))
      
      dat.C1 <- data.frame(x = x, 
                           t = t,
                           y = c(rep(0,maxL), y.C1))
      
      
      ##### MODEL FITTING #######
      ##### 1. C1 ##########
      # true: eta.C1
      ## 1.1 ACE-DLNM #########
      if(Nt == 1000) {
        ## 10 knots for s(t)
        datC1.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                 smooth = ~s(t, bs = "bs", k = 10),
                                 dat = dat.C1,
                                 pc = NULL,
                                 kw = 20, 
                                 kE = 20,
                                 maxL = maxL,
                                 eta = TRUE,
                                 verbose = FALSE)
      } else if (Nt == 2000) {
        ## 20 knots for s(t)
        datC1.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                 smooth = ~s(t, bs = "bs", k = 20),
                                 dat = dat.C1,
                                 kw = 20, 
                                 kE = 20,
                                 maxL = maxL,
                                 eta = TRUE,
                                 verbose = FALSE)
      }
      
      
      datC1.modDLNM.summ <- summary(datC1.modDLNM, plot = FALSE)
      
      datC1.modDLNM.eta_E.est <- datC1.modDLNM$eta_E$est
      datC1.modDLNM.eta_E.ul <- datC1.modDLNM$eta_E$ul
      datC1.modDLNM.eta_E.ll <- datC1.modDLNM$eta_E$ll
      
      datC1.modDLNM.eta_other.est <- datC1.modDLNM$eta_other$est
      datC1.modDLNM.eta_other.ul <- datC1.modDLNM$eta_other$ul
      datC1.modDLNM.eta_other.ll <- datC1.modDLNM$eta_other$ll
      
      ## 1.2 DRF-DLNM ###############
      # MATRIX Q OF EXPOSURE HISTORIES
      Q <- Lag(x[1:(Nt+maxL)],0:maxL)
      cb <- crossbasis(Q[,1],lag=c(0,maxL),argvar=list(fun='ps',df=19),
                       arglag=list(fun='ps', df=19))
      cbgamPen <- cbPen(cb)
      
      if (Nt == 1000) {
        datC1.modepDLNM <- gam(dat.C1$y~cb+s(dat.C1$t, bs = "bs", k = 10),family=nb(),paraPen=list(cb=cbgamPen),method='REML')  
      } else if (Nt == 2000) {
        datC1.modepDLNM <- gam(dat.C1$y~cb+s(dat.C1$t, bs = "bs", k = 20),family=nb(),paraPen=list(cb=cbgamPen),method='REML')
      }
      
      
      pred.lpmatrix <- predict(datC1.modepDLNM, type = "lpmatrix")
      
      datC1.modepDLNM.eta_E.est <- as.vector(pred.lpmatrix[, 1:362] %*% datC1.modepDLNM$coefficients[1:362])
      datC1.modepDLNM.eta_other.est <- as.vector(pred.lpmatrix[, -(1:362)] %*% datC1.modepDLNM$coefficients[-(1:362)])
      
      datC1.modepDLNM.eta_E.sd <- sqrt(diag(pred.lpmatrix[, 1:362] %*% vcov(datC1.modepDLNM)[1:362, 1:362] %*% t(pred.lpmatrix[, 1:362])))
      datC1.modepDLNM.eta_other.sd <- sqrt(diag(pred.lpmatrix[, -(1:362)] %*% vcov(datC1.modepDLNM)[-(1:362), -(1:362)] %*% t(pred.lpmatrix[, -(1:362)])))
      
      datC1.modepDLNM.eta_E.ul <- as.vector(datC1.modepDLNM.eta_E.est + 1.96*datC1.modepDLNM.eta_E.sd)
      datC1.modepDLNM.eta_E.ll <- as.vector(datC1.modepDLNM.eta_E.est - 1.96*datC1.modepDLNM.eta_E.sd)
      
      datC1.modepDLNM.eta_other.ul <- as.vector(datC1.modepDLNM.eta_other.est + 1.96*datC1.modepDLNM.eta_other.sd)
      datC1.modepDLNM.eta_other.ll <- as.vector(datC1.modepDLNM.eta_other.est - 1.96*datC1.modepDLNM.eta_other.sd)
      
      
      
      cat("finish: ", lst, "\n")
      #### RETURN ###########
      out = list(para = c(lst),
                 datC1 = list(E = list(true.eta = eta_E.C1,
                                       modDLNM = list(funs = datC1.modDLNM.summ$est,
                                                      est = datC1.modDLNM.eta_E.est,
                                                      ll = datC1.modDLNM.eta_E.ll,
                                                      ul = datC1.modDLNM.eta_E.ul),
                                       modpDLNM = list(est = datC1.modepDLNM.eta_E.est,
                                                       ll = datC1.modepDLNM.eta_E.ll,
                                                       ul = datC1.modepDLNM.eta_E.ul)),
                              other = list(true.eta = eta_other.C1,
                                           modDLNM = list(est = datC1.modDLNM.eta_other.est,
                                                          ll = datC1.modDLNM.eta_other.ll,
                                                          ul = datC1.modDLNM.eta_other.ul),
                                           modpDLNM = list(est = datC1.modepDLNM.eta_E.est,
                                                           ll = datC1.modepDLNM.eta_E.ll,
                                                           ul = datC1.modepDLNM.eta_E.ul))))
      return(out)
    },
    error = function(e){
      cat("error in ", lst, ".\n")
      warning(e)
    }
  )
}





## run the simulation
len <- length(para_run_list)
para_run_list1 <- para_run_list[0+(1:(len/4))]
para_run_list2 <- para_run_list[(1*len/4)+(1:(len/4))]
para_run_list3 <- para_run_list[(2*len/4)+(1:(len/4))]
para_run_list4 <- para_run_list[(3*len/4)+(1:(len/4))]



simresults1 <- mclapply(para_run_list1, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "05-simulation-appendix-results-1.RData"))
rm(simresults1)

simresults2 <- mclapply(para_run_list2, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "05-simulation-appendix-results-2.RData"))
rm(simresults2)

simresults3 <- mclapply(para_run_list3, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "05-simulation-appendix-results-3.RData"))
rm(simresults3)

simresults4 <- mclapply(para_run_list4, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "05-simulation-appendix-results-4.RData"))
rm(simresults4)



cat("DONE!")
