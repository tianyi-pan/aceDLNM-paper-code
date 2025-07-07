#!/usr/bin/env Rscript

###### Simulation A ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Section 4.2 Simulation B: Relative Performance
## This code is used to run the simulation study. To obtain the tables/figures in the paper, please see: 02-simulation-B-analysis.R
## Tianyi Pan
## 2025
###################################


## load packages
# devtools::install_github("tianyi-pan/aceDLNM")
library(aceDLNM)
library(parallel)
library(mgcv)
library(dlnm); library(splines) ; library(tsModel) # packages for RDF-DLNM (Gasparrini et.al. 2017)

## paths
options(mc.cores = parallel::detectCores()-1)
resultpath <- "results"
if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)



R <- 10000 # number of iterations; you could adjust this numbder to control the running time of this simulation

Nt.list <- c(2000)
para <- expand.grid(Nt = Nt.list,
                    iter = 1:R)
para_run_list <- as.list(data.frame(t(para))) # convert matrix to list by row

avglag.list <- c(0, 1, 7, 14) # generata data with lag0, lag0-1, lag0-7, or lag0-14


for(avglag in avglag.list){
  cat("Start avglag: ", avglag, ". \n")

  ## function to run the simulation
  run_lst <- function(lst){
    tryCatch(
      expr = {
        set.seed(as.numeric(lst[2]))

        ## generate data
        cat("start: ", lst, "\n")

        gt <- function(x) 0.5*sin(x/150)
        maxL <- 14
        theta <- 8

        Nt <- as.numeric(lst[1])

        ttmp <- 1:(maxL+Nt)
        ttmp <- ttmp[-(1:maxL)]
        meangt <- mean(gt(ttmp))
        otherterm <- data.frame(trend = gt(ttmp) - meangt,
                                intercept = 0)

        ## extract data from the pacakge
        t <- 1:(Nt + maxL)
        x <- PM25.waterloo$PM25[t]


        ### GENERATE DATA ########
        ## 1. generate data from GAM with unequal weights
        ## type (i) w
        wl <- function(l) dnorm(l, mean = 3, sd = 3.5)

        wl.discrete <- sapply(0:maxL, wl)
        wl.discrete <- wl.discrete/sum(wl.discrete)

        x.DLNM.true <- sapply(15:2014, function(i) sum(x[i - c(0:14)] * wl.discrete))

        x.DLNM.true.max <- ceiling(max(x.DLNM.true, na.rm = TRUE))
        x.DLNM.true.min <- floor(min(x.DLNM.true, na.rm = TRUE))

        ## type (i) f cubic 
        fEtmp.DLNM <- function(x) (x-(x.DLNM.true.min + (x.DLNM.true.max-x.DLNM.true.min)*0.3))*(x-(x.DLNM.true.min + (x.DLNM.true.max-x.DLNM.true.min)*0.2))*(x-(x.DLNM.true.min + (x.DLNM.true.max-x.DLNM.true.min)*0.9))
        fE.DLNM <- function(x) fEtmp.DLNM(x) / (fEtmp.DLNM(x.DLNM.true.max)/2.5) + 2.5

        eta.DLNM <- fE.DLNM(x.DLNM.true) + gt(ttmp) - meangt
        y.DLNM <- sapply(eta.DLNM, function(eta.) rnbinom(n = 1, size = theta, mu = exp(eta.)))

        eta_E.DLNM <- fE.DLNM(x.DLNM.true)
        eta_other.DLNM <- gt(ttmp) - meangt


        dat.DLNM <- data.frame(x = x,
                               t = t,
                               y = c(rep(0,maxL), y.DLNM))

        ## 2. generate data from GAM with chosen avglag
        # avglag: 0: lag0; 1: lag0-1; 7: lago-7; 14: lag0-14

        x.gam.true <- sapply(15:2014, function(i) mean(x[i - c(0:avglag)]))

        x.gam.true.max <- ceiling(max(x.gam.true, na.rm = TRUE))
        x.gam.true.min <- floor(min(x.gam.true, na.rm = TRUE))

        ## type (i) f cubic
        fEtmp.gam <- function(x) (x-(x.gam.true.min + (x.gam.true.max-x.gam.true.min)*0.3))*(x-(x.gam.true.min + (x.gam.true.max-x.gam.true.min)*0.2))*(x-(x.gam.true.min + (x.gam.true.max-x.gam.true.min)*0.9))
        fE.gam <- function(x) fEtmp.gam(x) / (fEtmp.gam(x.gam.true.max)/2.5) + 2.5


        eta.gam <- fE.gam(x.gam.true) + gt(ttmp) - meangt

        eta_E.gam <- fE.gam(x.gam.true)
        eta_other.gam <- gt(ttmp) - meangt

        y.gam <- sapply(eta.gam, function(eta.) rnbinom(n = 1, size = theta, mu = exp(eta.)))

        dat.GAM <- data.frame(x = x,
                              t = t,
                              y = c(rep(0,maxL), y.gam))



        ##### MODEL FITTING #######
        ##### 1. data from GAM ##########
        # true: eta.gam
        ## 1.1 DLNM with delta = FALSE

        ## 20 knots for s(t)
        # kw = 20
        datGAM.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                  smooth = ~s(t, bs = "bs", k = 20),
                                  dat = dat.GAM,
                                  kw = 20,
                                  kE = 20,
                                  maxL = maxL,
                                  eta = TRUE,
                                  verbose = FALSE)



        datGAM.modDLNM.summ <- summary(datGAM.modDLNM, plot = FALSE)

        datGAM.modDLNM.eta_E.est <- datGAM.modDLNM$eta_E$est
        datGAM.modDLNM.eta_E.ul <- datGAM.modDLNM$eta_E$ul
        datGAM.modDLNM.eta_E.ll <- datGAM.modDLNM$eta_E$ll

        datGAM.modDLNM.eta_other.est <- datGAM.modDLNM$eta_other$est
        datGAM.modDLNM.eta_other.ul <- datGAM.modDLNM$eta_other$ul
        datGAM.modDLNM.eta_other.ll <- datGAM.modDLNM$eta_other$ll

        rm("datGAM.modDLNM")


        ## 1.2 GAM with mean of lag-0 and lag-1
        dat1 <- dat.GAM
        x.lagged <- (x + dplyr::lag(x))/2
        dat1$x.mean01 <- x.lagged[1:(Nt+maxL)]
        dat1 <- dat1[-(1:maxL),]


        datGAM.modeGAM01 <- mgcv::gam(y ~ s(x.mean01, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                      data = dat1,
                                      family = nb())

        datGAM.modeGAM01.lpmatrix <- predict(datGAM.modeGAM01, type = "lpmatrix")
        datGAM.modeGAM01.beta <- datGAM.modeGAM01$coefficients


        datGAM.modeGAM01.eta_E.est <- datGAM.modeGAM01.lpmatrix[,1:20] %*% datGAM.modeGAM01.beta[1:20]
        datGAM.modeGAM01.eta_other.est <- datGAM.modeGAM01.lpmatrix[,21:39] %*% datGAM.modeGAM01.beta[21:39]


        datGAM.modeGAM01.eta_E.sd <- sqrt(diag(datGAM.modeGAM01.lpmatrix[,1:20] %*% vcov(datGAM.modeGAM01)[1:20,1:20] %*% t(datGAM.modeGAM01.lpmatrix[,1:20])))
        datGAM.modeGAM01.eta_other.sd <- sqrt(diag(datGAM.modeGAM01.lpmatrix[,21:39] %*% vcov(datGAM.modeGAM01)[21:39, 21:39] %*% t(datGAM.modeGAM01.lpmatrix[,21:39])))

        datGAM.modeGAM01.eta_E.ul <- datGAM.modeGAM01.eta_E.est + 1.96*datGAM.modeGAM01.eta_E.sd
        datGAM.modeGAM01.eta_E.ll <- datGAM.modeGAM01.eta_E.est - 1.96*datGAM.modeGAM01.eta_E.sd

        datGAM.modeGAM01.eta_other.ul <- datGAM.modeGAM01.eta_other.est + 1.96*datGAM.modeGAM01.eta_other.sd
        datGAM.modeGAM01.eta_other.ll <- datGAM.modeGAM01.eta_other.est - 1.96*datGAM.modeGAM01.eta_other.sd


        ## 1.3 GAM with mean over 0:14
        x.avgmaxL <- sapply((1:Nt)+maxL, function(i) mean(x[i - c(0:maxL)]) )
        dat2 <- dat.GAM
        dat2 <- dat2[-(1:maxL),]
        dat2$x.avgmaxL <- x.avgmaxL[1:Nt]



        datGAM.modeGAMavgmaxL <- mgcv::gam(y ~ s(x.avgmaxL, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                           data = dat2,
                                           family = nb())

        datGAM.modeGAMavgmaxL.lpmatrix <- predict(datGAM.modeGAMavgmaxL, type = "lpmatrix")
        datGAM.modeGAMavgmaxL.beta <- datGAM.modeGAMavgmaxL$coefficients


        datGAM.modeGAMavgmaxL.eta_E.est <- datGAM.modeGAMavgmaxL.lpmatrix[,1:20] %*% datGAM.modeGAMavgmaxL.beta[1:20]
        datGAM.modeGAMavgmaxL.eta_other.est <- datGAM.modeGAMavgmaxL.lpmatrix[,21:39] %*% datGAM.modeGAMavgmaxL.beta[21:39]


        datGAM.modeGAMavgmaxL.eta_E.sd <- sqrt(diag(datGAM.modeGAMavgmaxL.lpmatrix[,1:20] %*% vcov(datGAM.modeGAMavgmaxL)[1:20,1:20] %*% t(datGAM.modeGAMavgmaxL.lpmatrix[,1:20])))
        datGAM.modeGAMavgmaxL.eta_other.sd <- sqrt(diag(datGAM.modeGAMavgmaxL.lpmatrix[,21:39] %*% vcov(datGAM.modeGAMavgmaxL)[21:39, 21:39] %*% t(datGAM.modeGAMavgmaxL.lpmatrix[,21:39])))

        datGAM.modeGAMavgmaxL.eta_E.ul <- datGAM.modeGAMavgmaxL.eta_E.est + 1.96*datGAM.modeGAMavgmaxL.eta_E.sd
        datGAM.modeGAMavgmaxL.eta_E.ll <- datGAM.modeGAMavgmaxL.eta_E.est - 1.96*datGAM.modeGAMavgmaxL.eta_E.sd

        datGAM.modeGAMavgmaxL.eta_other.ul <- datGAM.modeGAMavgmaxL.eta_other.est + 1.96*datGAM.modeGAMavgmaxL.eta_other.sd
        datGAM.modeGAMavgmaxL.eta_other.ll <- datGAM.modeGAMavgmaxL.eta_other.est - 1.96*datGAM.modeGAMavgmaxL.eta_other.sd


        ## 1.4 GAM with mean over 0:7
        x.avg7 <- sapply((1:Nt)+maxL, function(i) mean(x[i - c(0:7)]) )
        dat3 <- dat.GAM
        dat3 <- dat3[-(1:maxL),]
        dat3$x.avg7 <- x.avg7[1:Nt]



        datGAM.modeGAMavg7 <- mgcv::gam(y ~ s(x.avg7, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                        data = dat3,
                                        family = nb())



        datGAM.modeGAMavg7.lpmatrix <- predict(datGAM.modeGAMavg7, type = "lpmatrix")
        datGAM.modeGAMavg7.beta <- datGAM.modeGAMavg7$coefficients


        datGAM.modeGAMavg7.eta_E.est <- datGAM.modeGAMavg7.lpmatrix[,1:20] %*% datGAM.modeGAMavg7.beta[1:20]
        datGAM.modeGAMavg7.eta_other.est <- datGAM.modeGAMavg7.lpmatrix[,21:39] %*% datGAM.modeGAMavg7.beta[21:39]


        datGAM.modeGAMavg7.eta_E.sd <- sqrt(diag(datGAM.modeGAMavg7.lpmatrix[,1:20] %*% vcov(datGAM.modeGAMavg7)[1:20,1:20] %*% t(datGAM.modeGAMavg7.lpmatrix[,1:20])))
        datGAM.modeGAMavg7.eta_other.sd <- sqrt(diag(datGAM.modeGAMavg7.lpmatrix[,21:39] %*% vcov(datGAM.modeGAMavg7)[21:39, 21:39] %*% t(datGAM.modeGAMavg7.lpmatrix[,21:39])))

        datGAM.modeGAMavg7.eta_E.ul <- datGAM.modeGAMavg7.eta_E.est + 1.96*datGAM.modeGAMavg7.eta_E.sd
        datGAM.modeGAMavg7.eta_E.ll <- datGAM.modeGAMavg7.eta_E.est - 1.96*datGAM.modeGAMavg7.eta_E.sd

        datGAM.modeGAMavg7.eta_other.ul <- datGAM.modeGAMavg7.eta_other.est + 1.96*datGAM.modeGAMavg7.eta_other.sd
        datGAM.modeGAMavg7.eta_other.ll <- datGAM.modeGAMavg7.eta_other.est - 1.96*datGAM.modeGAMavg7.eta_other.sd

        # mean((datGAM.modeGAMavg7.eta_E.est - eta_E.gam)^2)
        # mean((datGAM.modeGAMavg7.eta_E.ul >= eta_E.gam) & (datGAM.modeGAMavg7.eta_E.ll <= eta_E.gam))


        ## 1.5 GAM without lag
        dat4 <- dat.GAM
        dat4 <- dat4[-(1:maxL),]

        datGAM.modeGAM0 <- mgcv::gam(y ~ s(x, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                     data = dat4,
                                     family = nb())


        datGAM.modeGAM0.lpmatrix <- predict(datGAM.modeGAM0, type = "lpmatrix")
        datGAM.modeGAM0.beta <- datGAM.modeGAM0$coefficients


        datGAM.modeGAM0.eta_E.est <- datGAM.modeGAM0.lpmatrix[,1:20] %*% datGAM.modeGAM0.beta[1:20]
        datGAM.modeGAM0.eta_other.est <- datGAM.modeGAM0.lpmatrix[,21:39] %*% datGAM.modeGAM0.beta[21:39]


        datGAM.modeGAM0.eta_E.sd <- sqrt(diag(datGAM.modeGAM0.lpmatrix[,1:20] %*% vcov(datGAM.modeGAM0)[1:20,1:20] %*% t(datGAM.modeGAM0.lpmatrix[,1:20])))
        datGAM.modeGAM0.eta_other.sd <- sqrt(diag(datGAM.modeGAM0.lpmatrix[,21:39] %*% vcov(datGAM.modeGAM0)[21:39, 21:39] %*% t(datGAM.modeGAM0.lpmatrix[,21:39])))

        datGAM.modeGAM0.eta_E.ul <- datGAM.modeGAM0.eta_E.est + 1.96*datGAM.modeGAM0.eta_E.sd
        datGAM.modeGAM0.eta_E.ll <- datGAM.modeGAM0.eta_E.est - 1.96*datGAM.modeGAM0.eta_E.sd

        datGAM.modeGAM0.eta_other.ul <- datGAM.modeGAM0.eta_other.est + 1.96*datGAM.modeGAM0.eta_other.sd
        datGAM.modeGAM0.eta_other.ll <- datGAM.modeGAM0.eta_other.est - 1.96*datGAM.modeGAM0.eta_other.sd


        ## 1.6 DRF-DLNM (Gasparrini et.al. 2017)

        # MATRIX Q OF EXPOSURE HISTORIES
        Q <- Lag(x[1:(Nt+maxL)],0:maxL)
        cb <- crossbasis(Q[,1],lag=c(0,maxL),argvar=list(fun='ps',df=19),
                         arglag=list(fun='ps', df=19))
        cbgamPen <- cbPen(cb)


        datGAM.modepDLNM <- gam(dat.GAM$y~cb+s(dat.GAM$t, bs = "bs", k = 20),family=nb(),paraPen=list(cb=cbgamPen),method='REML')



        pred.lpmatrix <- predict(datGAM.modepDLNM, type = "lpmatrix")

        datGAM.modepDLNM.eta_E.est <- as.vector(pred.lpmatrix[, 1:362] %*% datGAM.modepDLNM$coefficients[1:362])
        datGAM.modepDLNM.eta_other.est <- as.vector(pred.lpmatrix[, -(1:362)] %*% datGAM.modepDLNM$coefficients[-(1:362)])

        datGAM.modepDLNM.eta_E.sd <- sqrt(diag(pred.lpmatrix[, 1:362] %*% vcov(datGAM.modepDLNM)[1:362, 1:362] %*% t(pred.lpmatrix[, 1:362])))
        datGAM.modepDLNM.eta_other.sd <- sqrt(diag(pred.lpmatrix[, -(1:362)] %*% vcov(datGAM.modepDLNM)[-(1:362), -(1:362)] %*% t(pred.lpmatrix[, -(1:362)])))

        datGAM.modepDLNM.eta_E.ul <- as.vector(datGAM.modepDLNM.eta_E.est + 1.96*datGAM.modepDLNM.eta_E.sd)
        datGAM.modepDLNM.eta_E.ll <- as.vector(datGAM.modepDLNM.eta_E.est - 1.96*datGAM.modepDLNM.eta_E.sd)

        datGAM.modepDLNM.eta_other.ul <- as.vector(datGAM.modepDLNM.eta_other.est + 1.96*datGAM.modepDLNM.eta_other.sd)
        datGAM.modepDLNM.eta_other.ll <- as.vector(datGAM.modepDLNM.eta_other.est - 1.96*datGAM.modepDLNM.eta_other.sd)

        if(avglag == 0) {
          ##### 2. data from GAM with unequal weights ##########

          ## 2.1 DLNM with delta = FALSE

          ## 20 knots for s(t)
          datDLNM.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                     smooth = ~s(t, bs = "bs", k = 20),
                                     dat = dat.DLNM,
                                     pc = NULL,
                                     kw = 20,
                                     kE = 20,
                                     maxL = maxL,
                                     eta = TRUE,
                                     verbose = FALSE)


          datDLNM.modDLNM.summ <- summary(datDLNM.modDLNM, plot = FALSE)

          datDLNM.modDLNM.eta_E.est <- datDLNM.modDLNM$eta_E$est
          datDLNM.modDLNM.eta_E.ul <- datDLNM.modDLNM$eta_E$ul
          datDLNM.modDLNM.eta_E.ll <- datDLNM.modDLNM$eta_E$ll

          datDLNM.modDLNM.eta_other.est <- datDLNM.modDLNM$eta_other$est
          datDLNM.modDLNM.eta_other.ul <- datDLNM.modDLNM$eta_other$ul
          datDLNM.modDLNM.eta_other.ll <- datDLNM.modDLNM$eta_other$ll

          rm("datDLNM.modDLNM")

          ## 2.2 GAM with mean of lag-0 and lag-1
          dat1 <- dat.DLNM
          x.lagged <- (x + dplyr::lag(x))/2
          dat1$x.mean01 <- x.lagged[1:(Nt+maxL)]
          dat1 <- dat1[-(1:maxL),]



          datDLNM.modeGAM01 <- mgcv::gam(y ~ s(x.mean01, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                         data = dat1,
                                         family = nb())


          datDLNM.modeGAM01.lpmatrix <- predict(datDLNM.modeGAM01, type = "lpmatrix")
          datDLNM.modeGAM01.beta <- datDLNM.modeGAM01$coefficients


          datDLNM.modeGAM01.eta_E.est <- datDLNM.modeGAM01.lpmatrix[,1:20] %*% datDLNM.modeGAM01.beta[1:20]
          datDLNM.modeGAM01.eta_other.est <- datDLNM.modeGAM01.lpmatrix[,21:39] %*% datDLNM.modeGAM01.beta[21:39]


          datDLNM.modeGAM01.eta_E.sd <- sqrt(diag(datDLNM.modeGAM01.lpmatrix[,1:20] %*% vcov(datDLNM.modeGAM01)[1:20,1:20] %*% t(datDLNM.modeGAM01.lpmatrix[,1:20])))
          datDLNM.modeGAM01.eta_other.sd <- sqrt(diag(datDLNM.modeGAM01.lpmatrix[,21:39] %*% vcov(datDLNM.modeGAM01)[21:39, 21:39] %*% t(datDLNM.modeGAM01.lpmatrix[,21:39])))

          datDLNM.modeGAM01.eta_E.ul <- datDLNM.modeGAM01.eta_E.est + 1.96*datDLNM.modeGAM01.eta_E.sd
          datDLNM.modeGAM01.eta_E.ll <- datDLNM.modeGAM01.eta_E.est - 1.96*datDLNM.modeGAM01.eta_E.sd

          datDLNM.modeGAM01.eta_other.ul <- datDLNM.modeGAM01.eta_other.est + 1.96*datDLNM.modeGAM01.eta_other.sd
          datDLNM.modeGAM01.eta_other.ll <- datDLNM.modeGAM01.eta_other.est - 1.96*datDLNM.modeGAM01.eta_other.sd



          ## 2.3 GAM with mean over 0:14
          x.avgmaxL <- sapply((1:Nt)+maxL, function(i) mean(x[i - c(0:maxL)]) )
          dat2 <- dat.DLNM
          dat2 <- dat2[-(1:maxL),]
          dat2$x.avgmaxL <- x.avgmaxL[1:Nt]



          datDLNM.modeGAMavgmaxL <- mgcv::gam(y ~ s(x.avgmaxL, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                              data = dat2,
                                              family = nb())


          datDLNM.modeGAMavgmaxL.lpmatrix <- predict(datDLNM.modeGAMavgmaxL, type = "lpmatrix")
          datDLNM.modeGAMavgmaxL.beta <- datDLNM.modeGAMavgmaxL$coefficients


          datDLNM.modeGAMavgmaxL.eta_E.est <- datDLNM.modeGAMavgmaxL.lpmatrix[,1:20] %*% datDLNM.modeGAMavgmaxL.beta[1:20]
          datDLNM.modeGAMavgmaxL.eta_other.est <- datDLNM.modeGAMavgmaxL.lpmatrix[,21:39] %*% datDLNM.modeGAMavgmaxL.beta[21:39]


          datDLNM.modeGAMavgmaxL.eta_E.sd <- sqrt(diag(datDLNM.modeGAMavgmaxL.lpmatrix[,1:20] %*% vcov(datDLNM.modeGAMavgmaxL)[1:20,1:20] %*% t(datDLNM.modeGAMavgmaxL.lpmatrix[,1:20])))
          datDLNM.modeGAMavgmaxL.eta_other.sd <- sqrt(diag(datDLNM.modeGAMavgmaxL.lpmatrix[,21:39] %*% vcov(datDLNM.modeGAMavgmaxL)[21:39, 21:39] %*% t(datDLNM.modeGAMavgmaxL.lpmatrix[,21:39])))

          datDLNM.modeGAMavgmaxL.eta_E.ul <- datDLNM.modeGAMavgmaxL.eta_E.est + 1.96*datDLNM.modeGAMavgmaxL.eta_E.sd
          datDLNM.modeGAMavgmaxL.eta_E.ll <- datDLNM.modeGAMavgmaxL.eta_E.est - 1.96*datDLNM.modeGAMavgmaxL.eta_E.sd

          datDLNM.modeGAMavgmaxL.eta_other.ul <- datDLNM.modeGAMavgmaxL.eta_other.est + 1.96*datDLNM.modeGAMavgmaxL.eta_other.sd
          datDLNM.modeGAMavgmaxL.eta_other.ll <- datDLNM.modeGAMavgmaxL.eta_other.est - 1.96*datDLNM.modeGAMavgmaxL.eta_other.sd



          ## 2.4 GAM with mean over 0:7
          x.avg7 <- sapply((1:Nt)+maxL, function(i) mean(x[i - c(0:7)]) )
          dat3 <- dat.DLNM
          dat3 <- dat3[-(1:maxL),]
          dat3$x.avg7 <- x.avg7[1:Nt]


          datDLNM.modeGAMavg7 <- mgcv::gam(y ~ s(x.avg7, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                           data = dat3,
                                           family = nb())


          datDLNM.modeGAMavg7.lpmatrix <- predict(datDLNM.modeGAMavg7, type = "lpmatrix")
          datDLNM.modeGAMavg7.beta <- datDLNM.modeGAMavg7$coefficients


          datDLNM.modeGAMavg7.eta_E.est <- datDLNM.modeGAMavg7.lpmatrix[,1:20] %*% datDLNM.modeGAMavg7.beta[1:20]
          datDLNM.modeGAMavg7.eta_other.est <- datDLNM.modeGAMavg7.lpmatrix[,21:39] %*% datDLNM.modeGAMavg7.beta[21:39]


          datDLNM.modeGAMavg7.eta_E.sd <- sqrt(diag(datDLNM.modeGAMavg7.lpmatrix[,1:20] %*% vcov(datDLNM.modeGAMavg7)[1:20,1:20] %*% t(datDLNM.modeGAMavg7.lpmatrix[,1:20])))
          datDLNM.modeGAMavg7.eta_other.sd <- sqrt(diag(datDLNM.modeGAMavg7.lpmatrix[,21:39] %*% vcov(datDLNM.modeGAMavg7)[21:39, 21:39] %*% t(datDLNM.modeGAMavg7.lpmatrix[,21:39])))

          datDLNM.modeGAMavg7.eta_E.ul <- datDLNM.modeGAMavg7.eta_E.est + 1.96*datDLNM.modeGAMavg7.eta_E.sd
          datDLNM.modeGAMavg7.eta_E.ll <- datDLNM.modeGAMavg7.eta_E.est - 1.96*datDLNM.modeGAMavg7.eta_E.sd

          datDLNM.modeGAMavg7.eta_other.ul <- datDLNM.modeGAMavg7.eta_other.est + 1.96*datDLNM.modeGAMavg7.eta_other.sd
          datDLNM.modeGAMavg7.eta_other.ll <- datDLNM.modeGAMavg7.eta_other.est - 1.96*datDLNM.modeGAMavg7.eta_other.sd



          ## 2.5 GAM without lag
          dat4 <- dat.DLNM
          dat4 <- dat4[-(1:maxL),]


          datDLNM.modeGAM0 <- mgcv::gam(y ~ s(x, bs = "bs", k = 20) + s(t, bs = "bs", k = 20),
                                        data = dat4,
                                        family = nb())


          datDLNM.modeGAM0.lpmatrix <- predict(datDLNM.modeGAM0, type = "lpmatrix")
          datDLNM.modeGAM0.beta <- datDLNM.modeGAM0$coefficients


          datDLNM.modeGAM0.eta_E.est <- datDLNM.modeGAM0.lpmatrix[,1:20] %*% datDLNM.modeGAM0.beta[1:20]
          datDLNM.modeGAM0.eta_other.est <- datDLNM.modeGAM0.lpmatrix[,21:39] %*% datDLNM.modeGAM0.beta[21:39]


          datDLNM.modeGAM0.eta_E.sd <- sqrt(diag(datDLNM.modeGAM0.lpmatrix[,1:20] %*% vcov(datDLNM.modeGAM0)[1:20,1:20] %*% t(datDLNM.modeGAM0.lpmatrix[,1:20])))
          datDLNM.modeGAM0.eta_other.sd <- sqrt(diag(datDLNM.modeGAM0.lpmatrix[,21:39] %*% vcov(datDLNM.modeGAM0)[21:39, 21:39] %*% t(datDLNM.modeGAM0.lpmatrix[,21:39])))

          datDLNM.modeGAM0.eta_E.ul <- datDLNM.modeGAM0.eta_E.est + 1.96*datDLNM.modeGAM0.eta_E.sd
          datDLNM.modeGAM0.eta_E.ll <- datDLNM.modeGAM0.eta_E.est - 1.96*datDLNM.modeGAM0.eta_E.sd

          datDLNM.modeGAM0.eta_other.ul <- datDLNM.modeGAM0.eta_other.est + 1.96*datDLNM.modeGAM0.eta_other.sd
          datDLNM.modeGAM0.eta_other.ll <- datDLNM.modeGAM0.eta_other.est - 1.96*datDLNM.modeGAM0.eta_other.sd



          ## 2.6 DRF-DLNM (Gasparini et.al. 2017)

          # MATRIX Q OF EXPOSURE HISTORIES
          Q <- Lag(x[1:(Nt+maxL)],0:maxL)
          cb <- crossbasis(Q[,1],lag=c(0,maxL),argvar=list(fun='ps',df=19),
                           arglag=list(fun='ps', df=19))
          cbgamPen <- cbPen(cb)


          datDLNM.modepDLNM <- gam(dat.DLNM$y~cb+s(dat.DLNM$t, bs = "bs", k = 20),family=nb(),paraPen=list(cb=cbgamPen),method='REML')


          pred.lpmatrix <- predict(datDLNM.modepDLNM, type = "lpmatrix")

          datDLNM.modepDLNM.eta_E.est <- as.vector(pred.lpmatrix[, 1:362] %*% datDLNM.modepDLNM$coefficients[1:362])
          datDLNM.modepDLNM.eta_other.est <- as.vector(pred.lpmatrix[, -(1:362)] %*% datDLNM.modepDLNM$coefficients[-(1:362)])

          datDLNM.modepDLNM.eta_E.sd <- sqrt(diag(pred.lpmatrix[, 1:362] %*% vcov(datDLNM.modepDLNM)[1:362, 1:362] %*% t(pred.lpmatrix[, 1:362])))
          datDLNM.modepDLNM.eta_other.sd <- sqrt(diag(pred.lpmatrix[, -(1:362)] %*% vcov(datDLNM.modepDLNM)[-(1:362), -(1:362)] %*% t(pred.lpmatrix[, -(1:362)])))

          datDLNM.modepDLNM.eta_E.ul <- as.vector(datDLNM.modepDLNM.eta_E.est + 1.96*datDLNM.modepDLNM.eta_E.sd)
          datDLNM.modepDLNM.eta_E.ll <- as.vector(datDLNM.modepDLNM.eta_E.est - 1.96*datDLNM.modepDLNM.eta_E.sd)

          datDLNM.modepDLNM.eta_other.ul <- as.vector(datDLNM.modepDLNM.eta_other.est + 1.96*datDLNM.modepDLNM.eta_other.sd)
          datDLNM.modepDLNM.eta_other.ll <- as.vector(datDLNM.modepDLNM.eta_other.est - 1.96*datDLNM.modepDLNM.eta_other.sd)

        }

        cat("finish: ", lst, "\n")
        #### RETURN ###########
        out = list(para = c(lst),
                   datGAM = list(E = list(true.eta = eta_E.gam,
                                          modDLNM = list(funs = datGAM.modDLNM.summ$est,
                                                         est = datGAM.modDLNM.eta_E.est,
                                                         ll = datGAM.modDLNM.eta_E.ll,
                                                         ul = datGAM.modDLNM.eta_E.ul),
                                          modGAM01 = list(est = datGAM.modeGAM01.eta_E.est,
                                                          ll = datGAM.modeGAM01.eta_E.ll,
                                                          ul = datGAM.modeGAM01.eta_E.ul),
                                          modGAMavgmaxL = list(est = datGAM.modeGAMavgmaxL.eta_E.est,
                                                               ll = datGAM.modeGAMavgmaxL.eta_E.ll,
                                                               ul = datGAM.modeGAMavgmaxL.eta_E.ul),
                                          modGAMavg7 = list(est = datGAM.modeGAMavg7.eta_E.est,
                                                            ll = datGAM.modeGAMavg7.eta_E.ll,
                                                            ul = datGAM.modeGAMavg7.eta_E.ul),
                                          modGAM0 = list(est = datGAM.modeGAM0.eta_E.est,
                                                         ll = datGAM.modeGAM0.eta_E.ll,
                                                         ul = datGAM.modeGAM0.eta_E.ul),
                                          modpDLNM = list(est = datGAM.modepDLNM.eta_E.est,
                                                          ll = datGAM.modepDLNM.eta_E.ll,
                                                          ul = datGAM.modepDLNM.eta_E.ul)),
                                 other = list(true.eta = eta_other.gam,
                                              modDLNM = list(funs = datGAM.modDLNM.summ$est,
                                                             est = datGAM.modDLNM.eta_other.est,
                                                             ll = datGAM.modDLNM.eta_other.ll,
                                                             ul = datGAM.modDLNM.eta_other.ul),
                                              modGAM01 = list(est = datGAM.modeGAM01.eta_other.est,
                                                              ll = datGAM.modeGAM01.eta_other.ll,
                                                              ul = datGAM.modeGAM01.eta_other.ul),
                                              modGAMavgmaxL = list(est = datGAM.modeGAMavgmaxL.eta_other.est,
                                                                   ll = datGAM.modeGAMavgmaxL.eta_other.ll,
                                                                   ul = datGAM.modeGAMavgmaxL.eta_other.ul),
                                              modGAMavg7 = list(est = datGAM.modeGAMavg7.eta_other.est,
                                                                ll = datGAM.modeGAMavg7.eta_other.ll,
                                                                ul = datGAM.modeGAMavg7.eta_other.ul),
                                              modGAM0 = list(est = datGAM.modeGAM0.eta_other.est,
                                                             ll = datGAM.modeGAM0.eta_other.ll,
                                                             ul = datGAM.modeGAM0.eta_other.ul),
                                              modpDLNM = list(est = datGAM.modepDLNM.eta_other.est,
                                                              ll = datGAM.modepDLNM.eta_other.ll,
                                                              ul = datGAM.modepDLNM.eta_other.ul))

                   ))


        if(avglag == 0) {
          ## scenario (v)
          out$datDLNM = list(E = list(true.eta = eta_E.DLNM,
                                      modDLNM = list(funs = datDLNM.modDLNM.summ$est,
                                                     est = datDLNM.modDLNM.eta_E.est,
                                                     ll = datDLNM.modDLNM.eta_E.ll,
                                                     ul = datDLNM.modDLNM.eta_E.ul),
                                      modGAM01 = list(est = datDLNM.modeGAM01.eta_E.est,
                                                      ll = datDLNM.modeGAM01.eta_E.ll,
                                                      ul = datDLNM.modeGAM01.eta_E.ul),
                                      modGAMavgmaxL = list(est = datDLNM.modeGAMavgmaxL.eta_E.est,
                                                           ll = datDLNM.modeGAMavgmaxL.eta_E.ll,
                                                           ul = datDLNM.modeGAMavgmaxL.eta_E.ul),
                                      modGAMavg7 = list(est = datDLNM.modeGAMavg7.eta_E.est,
                                                        ll = datDLNM.modeGAMavg7.eta_E.ll,
                                                        ul = datDLNM.modeGAMavg7.eta_E.ul),
                                      modGAM0 = list(est = datDLNM.modeGAM0.eta_E.est,
                                                     ll = datDLNM.modeGAM0.eta_E.ll,
                                                     ul = datDLNM.modeGAM0.eta_E.ul),
                                      modpDLNM = list(est = datDLNM.modepDLNM.eta_E.est,
                                                      ll = datDLNM.modepDLNM.eta_E.ll,
                                                      ul = datDLNM.modepDLNM.eta_E.ul)),
                             other = list(true.eta = eta_other.DLNM,
                                          modDLNM = list(funs = datDLNM.modDLNM.summ$est,
                                                         est = datDLNM.modDLNM.eta_other.est,
                                                         ll = datDLNM.modDLNM.eta_other.ll,
                                                         ul = datDLNM.modDLNM.eta_other.ul),
                                          modGAM01 = list(est = datDLNM.modeGAM01.eta_other.est,
                                                          ll = datDLNM.modeGAM01.eta_other.ll,
                                                          ul = datDLNM.modeGAM01.eta_other.ul),
                                          modGAMavgmaxL = list(est = datDLNM.modeGAMavgmaxL.eta_other.est,
                                                               ll = datDLNM.modeGAMavgmaxL.eta_other.ll,
                                                               ul = datDLNM.modeGAMavgmaxL.eta_other.ul),
                                          modGAMavg7 = list(est = datDLNM.modeGAMavg7.eta_other.est,
                                                            ll = datDLNM.modeGAMavg7.eta_other.ll,
                                                            ul = datDLNM.modeGAMavg7.eta_other.ul),
                                          modGAM0 = list(est = datDLNM.modeGAM0.eta_other.est,
                                                         ll = datDLNM.modeGAM0.eta_other.ll,
                                                         ul = datDLNM.modeGAM0.eta_other.ul),
                                          modpDLNM = list(est = datDLNM.modepDLNM.eta_other.est,
                                                          ll = datDLNM.modepDLNM.eta_other.ll,
                                                          ul = datDLNM.modepDLNM.eta_other.ul))

          )
        }
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
  save.image(file = file.path(resultpath, paste0("02-simulation-", as.character(avglag), "-results-1.RData")))
  rm(simresults1)

  simresults2 <- mclapply(para_run_list2, run_lst, mc.preschedule = FALSE)
  save.image(file = file.path(resultpath, paste0("02-simulation-", as.character(avglag), "-results-2.RData")))
  rm(simresults2)

  simresults3 <- mclapply(para_run_list3, run_lst, mc.preschedule = FALSE)
  save.image(file = file.path(resultpath, paste0("02-simulation-", as.character(avglag), "-results-3.RData")))
  rm(simresults3)

  simresults4 <- mclapply(para_run_list4, run_lst, mc.preschedule = FALSE)
  save.image(file = file.path(resultpath, paste0("02-simulation-", as.character(avglag), "-results-4.RData")))
  rm(simresults4)
 
}


cat("DONE!")
