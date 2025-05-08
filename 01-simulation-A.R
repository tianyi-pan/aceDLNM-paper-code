#!/usr/bin/env Rscript

###### Simulation A ###############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Section 4.1 Simulation A: Absolute Performance
## This code is used to run the simulation study. To obtain the tables/figures in the paper, please see: 01-simulation-A-analysis.R
## Tianyi Pan
## 2025
###################################

## load packages
library(aceDLNM)
library(parallel)
library(mgcv)


Nt.list <- c(1000, 2000) 
wltype.list <- c("type1", "type2", "type3")
fEtype.list <- c("cubic", "linear", "quadratic")
R <- 10000 # number of iterations; you could adjust this numbder to control the running time of this simulation

para <- expand.grid(Nt = Nt.list,
                    wltype = wltype.list,
                    fEtype = fEtype.list,
                    iter = (1:R)+5000)
                    # iter = 1:R)

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


      dat <- data.frame(x = dat.gen$x,
                        t = dat.gen$t,
                        y = dat.gen$y)
      
      ## fit ACE-DLNM
      ## sampling method (default)
      if(Nt == 1000) {
        results <- aceDLNM(formula = y~sX(t, x),
                              smooth = ~s(t, bs = "bs", k = 10),
                              dat = dat,
                              pc = NULL,
                              kw = 20,
                              kE = 20,
                              maxL = maxL,
                              verbose = FALSE)
      } else {
        results <- aceDLNM(formula = y~sX(t, x),
                              smooth = ~s(t, bs = "bs", k = 20),
                              dat = dat,
                              pc = NULL,
                              kw = 20,
                              kE = 20,
                              maxL = maxL,
                              verbose = FALSE)
      }
      ## delta method (for Web Appendix G)
      if(Nt == 1000) {
        results.delta <- aceDLNM(formula = y~sX(t, x),
                              smooth = ~s(t, bs = "bs", k = 10),
                              dat = dat,
                              pc = NULL,
                              kw = 20,
                              kE = 20,
                              maxL = maxL,
                              delta.method = TRUE,
                              verbose = FALSE)
      } else {
        results.delta <- aceDLNM(formula = y~sX(t, x),
                              smooth = ~s(t, bs = "bs", k = 20),
                              dat = dat,
                              pc = NULL,
                              kw = 20,
                              kE = 20,
                              maxL = maxL,
                              delta.method = TRUE,
                              verbose = FALSE)
      }

      true.function <- dat.gen$true.f
      true.function$smooth <- function(x,var){
        if(var == "t") gt(x) - meangt
      }


      results.summary <- summary(results,
                                 E.eval = seq(dat.gen$true.f$Emin, dat.gen$true.f$Emax, length.out = 500),
                                 true.function = true.function, plot = FALSE)
      results.delta.summary <- summary(results.delta,
                                 E.eval = seq(dat.gen$true.f$Emin, dat.gen$true.f$Emax, length.out = 500),
                                 true.function = true.function, plot = FALSE)


      #### RETURN ###########
      cat("finish: ", lst, "\n")
      return(list(para = c(lst),
                  model = list(wl = results.summary$est$wl,
                               fE = results.summary$est$fE,
                               gt = results.summary$est$smooth,
                               opt = results$opt,
                               convergence = results$opt$convergence,
                               opttime = results$opttime),
                  model.delta = list(wl = results.delta.summary$est$wl,
                               fE = results.delta.summary$est$fE,
                               gt = results.delta.summary$est$smooth,
                               opt = results.delta$opt,
                               convergence = results.delta$opt$convergence,
                               opttime = results.delta$opttime))
      )
    },
    error = function(e){
      cat("error in ", lst, ".\n")
      warning(e)
    })
}


## set paths
resultpath <- "results"
if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)


options(mc.cores = parallel::detectCores()-1)

## run the simulation
## split the para_run_list into 16 parts, to avoid memory issues in parallel computing
len <- length(para_run_list)
para_run_list1 <- para_run_list[0+(1:(len/16))]
para_run_list2 <- para_run_list[(1*len/16)+(1:(len/16))]
para_run_list3 <- para_run_list[(2*len/16)+(1:(len/16))]
para_run_list4 <- para_run_list[(3*len/16)+(1:(len/16))]
para_run_list5 <- para_run_list[(4*len/16)+(1:(len/16))]
para_run_list6 <- para_run_list[(5*len/16)+(1:(len/16))]
para_run_list7 <- para_run_list[(6*len/16)+(1:(len/16))]
para_run_list8 <- para_run_list[(7*len/16)+(1:(len/16))]
para_run_list9 <- para_run_list[(8*len/16)+(1:(len/16))]
para_run_list10 <- para_run_list[(9*len/16)+(1:(len/16))]
para_run_list11 <- para_run_list[(10*len/16)+(1:(len/16))]
para_run_list12 <- para_run_list[(11*len/16)+(1:(len/16))]
para_run_list13 <- para_run_list[(12*len/16)+(1:(len/16))]
para_run_list14 <- para_run_list[(13*len/16)+(1:(len/16))]
para_run_list15 <- para_run_list[(14*len/16)+(1:(len/16))]
para_run_list16 <- para_run_list[(15*len/16)+(1:(len/16))]



simresults1 <- mclapply(para_run_list1, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-1.RData"))
rm(simresults1)

simresults2 <- mclapply(para_run_list2, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-2.RData"))
rm(simresults2)

simresults3 <- mclapply(para_run_list3, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-3.RData"))
rm(simresults3)

simresults4 <- mclapply(para_run_list4, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-4.RData"))
rm(simresults4)

simresults5 <- mclapply(para_run_list5, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-5.RData"))
rm(simresults5)

simresults6 <- mclapply(para_run_list6, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-6.RData"))
rm(simresults6)

simresults7 <- mclapply(para_run_list7, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-7.RData"))
rm(simresults7)

simresults8 <- mclapply(para_run_list8, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-8.RData"))
rm(simresults8)

simresults9 <- mclapply(para_run_list9, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-9.RData"))
rm(simresults9)

simresults10 <- mclapply(para_run_list10, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-10.RData"))
rm(simresults10)

simresults11 <- mclapply(para_run_list11, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-11.RData"))
rm(simresults11)

simresults12 <- mclapply(para_run_list12, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-12.RData"))
rm(simresults12)

simresults13 <- mclapply(para_run_list13, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-13.RData"))
rm(simresults13)

simresults14 <- mclapply(para_run_list14, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-14.RData"))
rm(simresults14)

simresults15 <- mclapply(para_run_list15, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-15.RData"))
rm(simresults15)

simresults16 <- mclapply(para_run_list16, run_lst, mc.preschedule = FALSE)
save.image(file = file.path(resultpath, "01-simulation-results-16.RData"))
rm(simresults16)