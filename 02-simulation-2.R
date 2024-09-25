#!/usr/bin/env Rscript
library(aceDLNM)
library(parallel)
library(mgcv)

library(dlnm); library(splines) ; library(tsModel) # packages for dlnm (Gasparrini et.al. 2017)

resultpath <- "results"
if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)


Nt.list <- c(2000)
wltype.list <- c("type1")
fEtype.list <- c("quadratic")
R <- 5000

para <- expand.grid(Nt = Nt.list,
                    wltype = wltype.list,
                    fEtype = fEtype.list,
                    iter = 1:R)
para_run_list <- as.list(data.frame(t(para))) # convert matrix to list by row

## obtain the exposure process with continuous weight
dat.init <- GenerateData(fEtype = "quadratic", wltype = "type1",
                        Nt = 2*14+2000+100,
                        theta = 8,
                        maxL = 14,
                        interpolate = TRUE)

avglag.list <- c(0, 1, 7, 14)
for(avglag in avglag.list){
  cat("Start avglag: ", avglag, ". \n")
  # Xt <- function(tnew) as.numeric(aceDLNM:::Bsplinevec2(tnew, dat.init$knots_x, 4) %*% dat.init$alpha_x)
  Xt.deBoor <- function(tnew) as.numeric(aceDLNM:::deBoor(tnew, dat.init$knots_x, dat.init$alpha_x, 4))
  x.gam.true <- sapply((1:2114)+14, function(s) return(integrate(Vectorize(Xt.deBoor), lower = s-avglag-0.5, upper = s+0.5)$value / (avglag+1)))  
  

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
        
        
        ttmp <- 1:(2*maxL+Nt+100)
        ttmp <- ttmp[-(1:maxL)]
        meangt <- mean(gt(ttmp))
        otherterm <- data.frame(trend = gt(ttmp) - meangt,
                                intercept = 0)
        
        
        ### GENERATE DATA ########
        ## 1. generate data from DLNM
        dat.gen.DLNM <- GenerateData(fEtype = fEtype, wltype = wltype, 
                                     Nt = Nt+100+maxL,
                                     theta = theta,
                                     maxL = maxL,
                                     other = otherterm,
                                     interpolate = TRUE)
        x <- dat.gen.DLNM$x # PM2.5
        t <- dat.gen.DLNM$t # time
        
        ## all data
        dat.DLNM.all <- data.frame(x = dat.gen.DLNM$x, 
                                   t = dat.gen.DLNM$t,
                                   y = dat.gen.DLNM$y)
        
        ## data for test
        dat.DLNM.test <- dat.DLNM.all[-c(1:(Nt+maxL)), ]
        dat.DLNM.test <- subset(dat.DLNM.test, select = -c(y)) # delete outcome variable for testing
        
        ## data for train
        dat.DLNM <- dat.DLNM.all[1:(Nt+maxL),]
        
        
        eta.DLNM.all <- dat.gen.DLNM$eta.sim
        eta.DLNM <- eta.DLNM.all[1:Nt]
        eta.DLNM.test <- eta.DLNM.all[(length(eta.DLNM.all)-99):length(eta.DLNM.all)]
        
        
        ## 2. generate data from GAM with chosen avglag
        # avglag: 0: lag0; 1: lag0-1; 7: lago-7; 14: lag0-14
        
        x.gam.true.min <- min(x.gam.true, na.rm = TRUE)
        x.gam.true.max <- max(x.gam.true, na.rm = TRUE)
        
        fEtmp.gam <- function(x){25*(dnorm(((2*((x-x.gam.true.min)/(x.gam.true.max-x.gam.true.min)+0.18) - 1.1))))}
        fE.gam <- function(x) -fEtmp.gam(x)+11
        
        
        eta.gam.all <- fE.gam(x.gam.true) + gt(ttmp) - meangt
        eta.gam <- eta.gam.all[1:Nt]
        eta.gam.test <- eta.gam.all[(length(eta.gam.all)-99):length(eta.gam.all)]
        
        
        
        y.gam <- sapply(eta.gam, function(eta.) rnbinom(n = 1, size = theta, mu = exp(eta.)))
        
        # plot(density(y.gam))
        # lines(density(dat.DLNM$y), col= "red")
        dat.GAM <- data.frame(x = dat.DLNM$x, 
                              t = dat.DLNM$t,
                              y = c(rep(0,maxL), y.gam))
        dat.GAM.test <- data.frame(x = dat.DLNM.test$x,
                                   t = dat.DLNM.test$t)
        
        
        ##### MODEL FITTING #######
        ##### 1. data from GAM ##########
        # true: eta.gam
        ## 1.1 DLNM with delta = FALSE
        if(Nt == 1000) {
          ## 10 knots for s(t)
          datGAM.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                       smooth = ~s(t, bs = "bs", k = 10),
                                       dat = dat.GAM,
                                       pc = NULL,
                                       kw = 30, 
                                       kE = 20,
                                       maxL = maxL,
                                       eta = TRUE,
                                       verbose = FALSE)
        } else if (Nt == 2000) {
          ## 20 knots for s(t)
          datGAM.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                       smooth = ~s(t, bs = "bs", k = 20),
                                       dat = dat.GAM,
                                       pc = NULL,
                                       kw = 30, 
                                       kE = 20,
                                       maxL = maxL,
                                       eta = TRUE,
                                       verbose = FALSE)
        }
        
        datGAM.modDLNM.summ <- summary(datGAM.modDLNM, plot = FALSE)
        
        datGAM.modDLNM.eta.est <- datGAM.modDLNM$eta$est
        datGAM.modDLNM.eta.ul <- datGAM.modDLNM$eta$ul
        datGAM.modDLNM.eta.ll <- datGAM.modDLNM$eta$ll
        
        datGAM.modDLNM.test <- predict(datGAM.modDLNM, dat.GAM.test)
        datGAM.modDLNM.eta.test.est <- datGAM.modDLNM.test$est
        datGAM.modDLNM.eta.test.ul <- datGAM.modDLNM.test$ul
        datGAM.modDLNM.eta.test.ll <- datGAM.modDLNM.test$ll
        
        # mean((datGAM.modDLNM.eta.est - eta.gam)^2)
        # mean((datGAM.modDLNM.eta.test.est - eta.gam.test)^2)
        # mean((datGAM.modDLNM.eta.ul >= eta.gam) & (datGAM.modDLNM.eta.ll <= eta.gam))
        # mean((datGAM.modDLNM.eta.test.ul >= eta.gam.test) & (datGAM.modDLNM.eta.test.ll <= eta.gam.test))
        
        
        
        rm("datGAM.modDLNM")
        
        # true: eta.gam and eta.gam.test
        
        ## 1.2 GAM with mean of lag-0 and lag-1
        dat1 <- dat.GAM
        x.lagged <- (x + dplyr::lag(x))/2
        dat1$x.mean01 <- x.lagged[1:(Nt+maxL)]
        dat1 <- dat1[-(1:maxL),]
        
        dat1.test <- data.frame(t = dat.GAM.test$t,
                                x.mean01 = x.lagged[-(1:(Nt+maxL))])
        dat1.test <- dat1.test[-(1:maxL),]
        if(Nt == 1000) {
          datGAM.modeGAM01 <- mgcv::gam(y ~ s(x.mean01, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                        data = dat1, 
                                        family = nb())
        } else if (Nt == 2000) {
          datGAM.modeGAM01 <- mgcv::gam(y ~ s(x.mean01, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                        data = dat1, 
                                        family = nb())
        }
        
        datGAM.modeGAM01.pred <- predict(datGAM.modeGAM01, type = "link", se.fit = TRUE)
        
        datGAM.modeGAM01.eta.est <- as.vector(datGAM.modeGAM01.pred$fit)
        datGAM.modeGAM01.eta.sd <- as.vector(datGAM.modeGAM01.pred$se.fit)
        datGAM.modeGAM01.eta.ul <- datGAM.modeGAM01.eta.est + 1.96*datGAM.modeGAM01.eta.sd
        datGAM.modeGAM01.eta.ll <- datGAM.modeGAM01.eta.est - 1.96*datGAM.modeGAM01.eta.sd
        
        datGAM.modeGAM01.test <- predict(datGAM.modeGAM01, newdata = dat1.test, type = "link", se.fit = TRUE)
        datGAM.modeGAM01.eta.test.est <- as.vector(datGAM.modeGAM01.test$fit)
        datGAM.modeGAM01.eta.test.sd <- as.vector(datGAM.modeGAM01.test$se.fit)
        datGAM.modeGAM01.eta.test.ul <- datGAM.modeGAM01.eta.test.est + 1.96*datGAM.modeGAM01.eta.test.sd
        datGAM.modeGAM01.eta.test.ll <- datGAM.modeGAM01.eta.test.est - 1.96*datGAM.modeGAM01.eta.test.sd
        
        # mean((datGAM.modeGAM01.eta.est - eta.gam)^2)
        # mean((datGAM.modeGAM01.eta.test.est - eta.gam.test)^2)
        # mean((datGAM.modeGAM01.eta.ul >= eta.gam) & (datGAM.modeGAM01.eta.ll <= eta.gam))
        # mean((datGAM.modeGAM01.eta.test.ul >= eta.gam.test) & (datGAM.modeGAM01.eta.test.ll <= eta.gam.test))
        
        
        ## 1.3 GAM with mean over 0:14
        x.avgmaxL <- sapply((1:(Nt+100+maxL))+maxL, function(i) mean(x[i - c(0:maxL)]) )
        dat2 <- dat.GAM
        dat2 <- dat2[-(1:maxL),]
        dat2$x.avgmaxL <- x.avgmaxL[1:Nt]
        
        
        
        dat2.test <- data.frame(t = dat.GAM.test$t,
                                x.avgmaxL = x.avgmaxL[-(1:Nt)])
        dat2.test <- dat2.test[-(1:maxL),]
        
        if (Nt == 1000) {
          datGAM.modeGAMavgmaxL <- mgcv::gam(y ~ s(x.avgmaxL, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                             data = dat2, 
                                             family = nb())
        } else if (Nt == 2000) {
          datGAM.modeGAMavgmaxL <- mgcv::gam(y ~ s(x.avgmaxL, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                             data = dat2, 
                                             family = nb())
        }
        
        datGAM.modeGAMavgmaxL.pred <- predict(datGAM.modeGAMavgmaxL, type = "link", se.fit = TRUE)
        
        datGAM.modeGAMavgmaxL.eta.est <- as.vector(datGAM.modeGAMavgmaxL.pred$fit)
        datGAM.modeGAMavgmaxL.eta.sd <- as.vector(datGAM.modeGAMavgmaxL.pred$se.fit)
        datGAM.modeGAMavgmaxL.eta.ul <- datGAM.modeGAMavgmaxL.eta.est + 1.96*datGAM.modeGAMavgmaxL.eta.sd
        datGAM.modeGAMavgmaxL.eta.ll <- datGAM.modeGAMavgmaxL.eta.est - 1.96*datGAM.modeGAMavgmaxL.eta.sd
        
        datGAM.modeGAMavgmaxL.test <- predict(datGAM.modeGAMavgmaxL, newdata = dat2.test, type = "link", se.fit = TRUE)
        datGAM.modeGAMavgmaxL.eta.test.est <- as.vector(datGAM.modeGAMavgmaxL.test$fit)
        datGAM.modeGAMavgmaxL.eta.test.sd <- as.vector(datGAM.modeGAMavgmaxL.test$se.fit)
        datGAM.modeGAMavgmaxL.eta.test.ul <- datGAM.modeGAMavgmaxL.eta.test.est + 1.96*datGAM.modeGAMavgmaxL.eta.test.sd
        datGAM.modeGAMavgmaxL.eta.test.ll <- datGAM.modeGAMavgmaxL.eta.test.est - 1.96*datGAM.modeGAMavgmaxL.eta.test.sd
        
        # mean((datGAM.modeGAMavgmaxL.eta.ul >= eta.gam) & (datGAM.modeGAMavgmaxL.eta.ll <= eta.gam))
        # mean((datGAM.modeGAMavgmaxL.eta.test.ul >= eta.gam.test) & (datGAM.modeGAMavgmaxL.eta.test.ll <= eta.gam.test))
        
        ## 1.4 GAM with mean over 0:7
        x.avg7 <- sapply((1:(Nt+100+maxL))+maxL, function(i) mean(x[i - c(0:7)]) )
        dat3 <- dat.GAM
        dat3 <- dat3[-(1:maxL),]
        dat3$x.avg7 <- x.avg7[1:Nt]
        
        
        
        dat3.test <- data.frame(t = dat.GAM.test$t,
                                x.avg7 = x.avg7[-(1:Nt)])
        dat3.test <- dat3.test[-(1:maxL),]
        
        if (Nt == 1000) {
          datGAM.modeGAMavg7 <- mgcv::gam(y ~ s(x.avg7, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                          data = dat3, 
                                          family = nb())
        } else if (Nt == 2000) {
          datGAM.modeGAMavg7 <- mgcv::gam(y ~ s(x.avg7, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                          data = dat3, 
                                          family = nb())
        }
        
        datGAM.modeGAMavg7.pred <- predict(datGAM.modeGAMavg7, type = "link", se.fit = TRUE)
        
        datGAM.modeGAMavg7.eta.est <- as.vector(datGAM.modeGAMavg7.pred$fit)
        datGAM.modeGAMavg7.eta.sd <- as.vector(datGAM.modeGAMavg7.pred$se.fit)
        datGAM.modeGAMavg7.eta.ul <- datGAM.modeGAMavg7.eta.est + 1.96*datGAM.modeGAMavg7.eta.sd
        datGAM.modeGAMavg7.eta.ll <- datGAM.modeGAMavg7.eta.est - 1.96*datGAM.modeGAMavg7.eta.sd
        
        datGAM.modeGAMavg7.test <- predict(datGAM.modeGAMavg7, newdata = dat3.test, type = "link", se.fit = TRUE)
        datGAM.modeGAMavg7.eta.test.est <- as.vector(datGAM.modeGAMavg7.test$fit)
        datGAM.modeGAMavg7.eta.test.sd <- as.vector(datGAM.modeGAMavg7.test$se.fit)
        datGAM.modeGAMavg7.eta.test.ul <- datGAM.modeGAMavg7.eta.test.est + 1.96*datGAM.modeGAMavg7.eta.test.sd
        datGAM.modeGAMavg7.eta.test.ll <- datGAM.modeGAMavg7.eta.test.est - 1.96*datGAM.modeGAMavg7.eta.test.sd
        
        # mean((datGAM.modeGAMavg7.eta.ul >= eta.gam) & (datGAM.modeGAMavg7.eta.ll <= eta.gam))
        # mean((datGAM.modeGAMavg7.eta.test.ul >= eta.gam.test) & (datGAM.modeGAMavg7.eta.test.ll <= eta.gam.test))
        
        
        ## 1.5 GAM without lag
        dat4 <- dat.GAM
        dat4 <- dat4[-(1:maxL),]
        
        
        dat4.test <- dat.GAM.test
        dat4.test <- dat4.test[-(1:maxL),]
        
        if (Nt == 1000) {
          datGAM.modeGAM0 <- mgcv::gam(y ~ s(x, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                       data = dat3, 
                                       family = nb())
        } else if (Nt == 2000) {
          datGAM.modeGAM0 <- mgcv::gam(y ~ s(x, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                       data = dat3, 
                                       family = nb())
        }
        
        datGAM.modeGAM0.pred <- predict(datGAM.modeGAM0, type = "link", se.fit = TRUE)
        
        datGAM.modeGAM0.eta.est <- as.vector(datGAM.modeGAM0.pred$fit)
        datGAM.modeGAM0.eta.sd <- as.vector(datGAM.modeGAM0.pred$se.fit)
        datGAM.modeGAM0.eta.ul <- datGAM.modeGAM0.eta.est + 1.96*datGAM.modeGAM0.eta.sd
        datGAM.modeGAM0.eta.ll <- datGAM.modeGAM0.eta.est - 1.96*datGAM.modeGAM0.eta.sd
        
        datGAM.modeGAM0.test <- predict(datGAM.modeGAM0, newdata = dat4.test, type = "link", se.fit = TRUE)
        datGAM.modeGAM0.eta.test.est <- as.vector(datGAM.modeGAM0.test$fit)
        datGAM.modeGAM0.eta.test.sd <- as.vector(datGAM.modeGAM0.test$se.fit)
        datGAM.modeGAM0.eta.test.ul <- datGAM.modeGAM0.eta.test.est + 1.96*datGAM.modeGAM0.eta.test.sd
        datGAM.modeGAM0.eta.test.ll <- datGAM.modeGAM0.eta.test.est - 1.96*datGAM.modeGAM0.eta.test.sd
        
        ## 1.6 Type (a) DLNM (Gasparini)
        x.norm <- (x-min(x))/diff(range(x))*10
        
        # MATRIX Q OF EXPOSURE HISTORIES
        Q <- Lag(x.norm[1:(Nt+maxL)],0:maxL)
        cb <- crossbasis(Q[,1],lag=c(0,maxL),argvar=list(fun='ps',df=19),
                         arglag=list(fun='ps', df=29))
        cbgamPen <- cbPen(cb)
        if (Nt == 1000) {
          datGAM.modepDLNM <- gam(dat.GAM$y~cb+s(dat.GAM$t, bs = "bs", k = 10),family=nb(),paraPen=list(cb=cbgamPen),method='REML')  
        } else if (Nt == 2000) {
          datGAM.modepDLNM <- gam(dat.GAM$y~cb+s(dat.GAM$t, bs = "bs", k = 20),family=nb(),paraPen=list(cb=cbgamPen),method='REML')
        }
        
        
        datGAM.modepDLNM.pred.object <- predict(datGAM.modepDLNM, type = "link", se.fit = TRUE)
        
        datGAM.modepDLNM.eta.est <- as.vector(datGAM.modepDLNM.pred.object$fit)
        datGAM.modepDLNM.eta.ul <- as.vector(datGAM.modepDLNM.eta.est + 1.96*datGAM.modepDLNM.pred.object$se.fit)
        datGAM.modepDLNM.eta.ll <- as.vector(datGAM.modepDLNM.eta.est - 1.96*datGAM.modepDLNM.pred.object$se.fit)
        
        
        # mean((datGAM.modepDLNM.eta.est - eta.gam)^2)
        # mean((datGAM.modepDLNM.eta.ul >= eta.gam) & (datGAM.modepDLNM.eta.ll <= eta.gam))
        
        
        if(avglag == 0) {
          ##### 2. data from DLNM ##########
          ## 2.1 DLNM with delta = FALSE
          if(Nt == 1000) {
            ## 10 knots for s(t)
            datDLNM.modDLNM <- aceDLNM(formula = y~sX(t, x),
                                          smooth = ~s(t, bs = "bs", k = 10),
                                          dat = dat.DLNM,
                                          pc = NULL,
                                          kw = 20, 
                                          kE = 20,
                                          maxL = maxL,
                                          eta = TRUE,
                                          verbose = FALSE)
          } else if (Nt == 2000) {
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
          }
          
          datDLNM.modDLNM.summ <- summary(datDLNM.modDLNM, plot = FALSE)
          
          datDLNM.modDLNM.eta.est <- datDLNM.modDLNM$eta$est
          datDLNM.modDLNM.eta.ul <- datDLNM.modDLNM$eta$ul
          datDLNM.modDLNM.eta.ll <- datDLNM.modDLNM$eta$ll
          
          datDLNM.modDLNM.test <- predict(datDLNM.modDLNM, dat.DLNM.test)
          datDLNM.modDLNM.eta.test.est <- datDLNM.modDLNM.test$est
          datDLNM.modDLNM.eta.test.ul <- datDLNM.modDLNM.test$ul
          datDLNM.modDLNM.eta.test.ll <- datDLNM.modDLNM.test$ll
          
          # mean((datDLNM.modDLNM.eta.est - eta.DLNM)^2)
          # mean((datDLNM.modDLNM.eta.test.est - eta.DLNM.test)^2)
          # mean((datDLNM.modDLNM.eta.ul >= eta.DLNM) & (datDLNM.modDLNM.eta.ll <= eta.DLNM))
          # mean((datDLNM.modDLNM.eta.test.ul >= eta.DLNM.test) & (datDLNM.modDLNM.eta.test.ll <= eta.DLNM.test))
          
          rm("datDLNM.modDLNM")
          
          # true: eta.DLNM and eta.DLNM.test
          
          ## 2.2 GAM with mean of lag-0 and lag-1
          dat1 <- dat.DLNM
          x.lagged <- (x + dplyr::lag(x))/2
          dat1$x.mean01 <- x.lagged[1:(Nt+maxL)]
          dat1 <- dat1[-(1:maxL),]
          
          dat1.test <- data.frame(t = dat.DLNM.test$t,
                                  x.mean01 = x.lagged[-(1:(Nt+maxL))])
          dat1.test <- dat1.test[-(1:maxL),]
          if(Nt == 1000) {
            datDLNM.modeGAM01 <- mgcv::gam(y ~ s(x.mean01, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                           data = dat1, 
                                           family = nb())
          } else if (Nt == 2000) {
            datDLNM.modeGAM01 <- mgcv::gam(y ~ s(x.mean01, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                           data = dat1, 
                                           family = nb())
          }
          
          datDLNM.modeGAM01.pred <- predict(datDLNM.modeGAM01, type = "link", se.fit = TRUE)
          
          datDLNM.modeGAM01.eta.est <- as.vector(datDLNM.modeGAM01.pred$fit)
          datDLNM.modeGAM01.eta.sd <- as.vector(datDLNM.modeGAM01.pred$se.fit)
          datDLNM.modeGAM01.eta.ul <- datDLNM.modeGAM01.eta.est + 1.96*datDLNM.modeGAM01.eta.sd
          datDLNM.modeGAM01.eta.ll <- datDLNM.modeGAM01.eta.est - 1.96*datDLNM.modeGAM01.eta.sd
          
          datDLNM.modeGAM01.test <- predict(datDLNM.modeGAM01, newdata = dat1.test, type = "link", se.fit = TRUE)
          datDLNM.modeGAM01.eta.test.est <- as.vector(datDLNM.modeGAM01.test$fit)
          datDLNM.modeGAM01.eta.test.sd <- as.vector(datDLNM.modeGAM01.test$se.fit)
          datDLNM.modeGAM01.eta.test.ul <- datDLNM.modeGAM01.eta.test.est + 1.96*datDLNM.modeGAM01.eta.test.sd
          datDLNM.modeGAM01.eta.test.ll <- datDLNM.modeGAM01.eta.test.est - 1.96*datDLNM.modeGAM01.eta.test.sd
          
          # mean((datDLNM.modeGAM01.eta.est - eta.DLNM)^2)
          # mean((datDLNM.modeGAM01.eta.test.est - eta.DLNM.test)^2)
          # mean((datDLNM.modeGAM01.eta.ul >= eta.DLNM) & (datDLNM.modeGAM01.eta.ll <= eta.DLNM))
          # mean((datDLNM.modeGAM01.eta.test.ul >= eta.DLNM.test) & (datDLNM.modeGAM01.eta.test.ll <= eta.DLNM.test))
          
          
          ## 2.3 GAM with mean over 0:14
          x.avgmaxL <- sapply((1:(Nt+100+maxL))+maxL, function(i) mean(x[i - c(0:maxL)]) )
          dat2 <- dat.DLNM
          dat2 <- dat2[-(1:maxL),]
          dat2$x.avgmaxL <- x.avgmaxL[1:Nt]
          
          
          
          dat2.test <- data.frame(t = dat.DLNM.test$t,
                                  x.avgmaxL = x.avgmaxL[-(1:Nt)])
          dat2.test <- dat2.test[-(1:maxL),]
          
          if (Nt == 1000) {
            datDLNM.modeGAMavgmaxL <- mgcv::gam(y ~ s(x.avgmaxL, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                                data = dat2, 
                                                family = nb())
          } else if (Nt == 2000) {
            datDLNM.modeGAMavgmaxL <- mgcv::gam(y ~ s(x.avgmaxL, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                                data = dat2, 
                                                family = nb())
          }
          
          datDLNM.modeGAMavgmaxL.pred <- predict(datDLNM.modeGAMavgmaxL, type = "link", se.fit = TRUE)
          
          datDLNM.modeGAMavgmaxL.eta.est <- as.vector(datDLNM.modeGAMavgmaxL.pred$fit)
          datDLNM.modeGAMavgmaxL.eta.sd <- as.vector(datDLNM.modeGAMavgmaxL.pred$se.fit)
          datDLNM.modeGAMavgmaxL.eta.ul <- datDLNM.modeGAMavgmaxL.eta.est + 1.96*datDLNM.modeGAMavgmaxL.eta.sd
          datDLNM.modeGAMavgmaxL.eta.ll <- datDLNM.modeGAMavgmaxL.eta.est - 1.96*datDLNM.modeGAMavgmaxL.eta.sd
          
          datDLNM.modeGAMavgmaxL.test <- predict(datDLNM.modeGAMavgmaxL, newdata = dat2.test, type = "link", se.fit = TRUE)
          datDLNM.modeGAMavgmaxL.eta.test.est <- as.vector(datDLNM.modeGAMavgmaxL.test$fit)
          datDLNM.modeGAMavgmaxL.eta.test.sd <- as.vector(datDLNM.modeGAMavgmaxL.test$se.fit)
          datDLNM.modeGAMavgmaxL.eta.test.ul <- datDLNM.modeGAMavgmaxL.eta.test.est + 1.96*datDLNM.modeGAMavgmaxL.eta.test.sd
          datDLNM.modeGAMavgmaxL.eta.test.ll <- datDLNM.modeGAMavgmaxL.eta.test.est - 1.96*datDLNM.modeGAMavgmaxL.eta.test.sd
          
          # mean((datDLNM.modeGAMavgmaxL.eta.est - eta.DLNM)^2)
          # mean((datDLNM.modeGAMavgmaxL.eta.test.est - eta.DLNM.test)^2)
          # mean((datDLNM.modeGAMavgmaxL.eta.ul >= eta.DLNM) & (datDLNM.modeGAMavgmaxL.eta.ll <= eta.DLNM))
          # mean((datDLNM.modeGAMavgmaxL.eta.test.ul >= eta.DLNM.test) & (datDLNM.modeGAMavgmaxL.eta.test.ll <= eta.DLNM.test))
          
          
          
          ## 2.4 GAM with mean over 0:7
          x.avg7 <- sapply((1:(Nt+100+maxL))+maxL, function(i) mean(x[i - c(0:7)]) )
          dat3 <- dat.DLNM
          dat3 <- dat3[-(1:maxL),]
          dat3$x.avg7 <- x.avg7[1:Nt]
          
          
          
          dat3.test <- data.frame(t = dat.DLNM.test$t,
                                  x.avg7 = x.avg7[-(1:Nt)])
          dat3.test <- dat3.test[-(1:maxL),]
          
          if (Nt == 1000) {
            datDLNM.modeGAMavg7 <- mgcv::gam(y ~ s(x.avg7, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                             data = dat3, 
                                             family = nb())
          } else if (Nt == 2000) {
            datDLNM.modeGAMavg7 <- mgcv::gam(y ~ s(x.avg7, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                             data = dat3, 
                                             family = nb())
          }
          
          datDLNM.modeGAMavg7.pred <- predict(datDLNM.modeGAMavg7, type = "link", se.fit = TRUE)
          
          datDLNM.modeGAMavg7.eta.est <- as.vector(datDLNM.modeGAMavg7.pred$fit)
          datDLNM.modeGAMavg7.eta.sd <- as.vector(datDLNM.modeGAMavg7.pred$se.fit)
          datDLNM.modeGAMavg7.eta.ul <- datDLNM.modeGAMavg7.eta.est + 1.96*datDLNM.modeGAMavg7.eta.sd
          datDLNM.modeGAMavg7.eta.ll <- datDLNM.modeGAMavg7.eta.est - 1.96*datDLNM.modeGAMavg7.eta.sd
          
          datDLNM.modeGAMavg7.test <- predict(datDLNM.modeGAMavg7, newdata = dat3.test, type = "link", se.fit = TRUE)
          datDLNM.modeGAMavg7.eta.test.est <- as.vector(datDLNM.modeGAMavg7.test$fit)
          datDLNM.modeGAMavg7.eta.test.sd <- as.vector(datDLNM.modeGAMavg7.test$se.fit)
          datDLNM.modeGAMavg7.eta.test.ul <- datDLNM.modeGAMavg7.eta.test.est + 1.96*datDLNM.modeGAMavg7.eta.test.sd
          datDLNM.modeGAMavg7.eta.test.ll <- datDLNM.modeGAMavg7.eta.test.est - 1.96*datDLNM.modeGAMavg7.eta.test.sd
          
          # mean((datDLNM.modeGAMavg7.eta.est - eta.DLNM)^2)
          # mean((datDLNM.modeGAMavg7.eta.test.est - eta.DLNM.test)^2)
          # mean((datDLNM.modeGAMavg7.eta.ul >= eta.DLNM) & (datDLNM.modeGAMavg7.eta.ll <= eta.DLNM))
          # mean((datDLNM.modeGAMavg7.eta.test.ul >= eta.DLNM.test) & (datDLNM.modeGAMavg7.eta.test.ll <= eta.DLNM.test))
          
          
          
          
          ## 2.5 GAM without lag
          dat4 <- dat.DLNM
          dat4 <- dat4[-(1:maxL),]
          
          
          
          dat4.test <- dat.DLNM.test
          dat4.test <- dat4.test[-(1:maxL),]
          
          if (Nt == 1000) {
            datDLNM.modeGAM0 <- mgcv::gam(y ~ s(x, bs = "bs", k = 20) + s(t, bs = "bs", k = 10), 
                                          data = dat4, 
                                          family = nb())
          } else if (Nt == 2000) {
            datDLNM.modeGAM0 <- mgcv::gam(y ~ s(x, bs = "bs", k = 20) + s(t, bs = "bs", k = 20), 
                                          data = dat4, 
                                          family = nb())
          }
          
          datDLNM.modeGAM0.pred <- predict(datDLNM.modeGAM0, type = "link", se.fit = TRUE)
          
          datDLNM.modeGAM0.eta.est <- as.vector(datDLNM.modeGAM0.pred$fit)
          datDLNM.modeGAM0.eta.sd <- as.vector(datDLNM.modeGAM0.pred$se.fit)
          datDLNM.modeGAM0.eta.ul <- datDLNM.modeGAM0.eta.est + 1.96*datDLNM.modeGAM0.eta.sd
          datDLNM.modeGAM0.eta.ll <- datDLNM.modeGAM0.eta.est - 1.96*datDLNM.modeGAM0.eta.sd
          
          datDLNM.modeGAM0.test <- predict(datDLNM.modeGAM0, newdata = dat4.test, type = "link", se.fit = TRUE)
          datDLNM.modeGAM0.eta.test.est <- as.vector(datDLNM.modeGAM0.test$fit)
          datDLNM.modeGAM0.eta.test.sd <- as.vector(datDLNM.modeGAM0.test$se.fit)
          datDLNM.modeGAM0.eta.test.ul <- datDLNM.modeGAM0.eta.test.est + 1.96*datDLNM.modeGAM0.eta.test.sd
          datDLNM.modeGAM0.eta.test.ll <- datDLNM.modeGAM0.eta.test.est - 1.96*datDLNM.modeGAM0.eta.test.sd
          
          ## 2.6 Type (a) DLNM (Gasparini)
          x.norm <- (x-min(x))/diff(range(x))*10
          
          # MATRIX Q OF EXPOSURE HISTORIES
          Q <- Lag(x.norm[1:(Nt+maxL)],0:maxL)
          cb <- crossbasis(Q[,1],lag=c(0,maxL),argvar=list(fun='ps',df=19),
                           arglag=list(fun='ps', df=19))
          cbgamPen <- cbPen(cb)
          if (Nt == 1000) {
            datDLNM.modepDLNM <- gam(dat.DLNM$y~cb+s(dat.DLNM$t, bs = "bs", k = 10),family=nb(),paraPen=list(cb=cbgamPen),method='REML')  
          } else if (Nt == 2000) {
            datDLNM.modepDLNM <- gam(dat.DLNM$y~cb+s(dat.DLNM$t, bs = "bs", k = 20),family=nb(),paraPen=list(cb=cbgamPen),method='REML')
          }
          
          
          datDLNM.modepDLNM.pred.object <- predict(datDLNM.modepDLNM, type = "link", se.fit = TRUE)
          
          datDLNM.modepDLNM.eta.est <- as.vector(datDLNM.modepDLNM.pred.object$fit)
          datDLNM.modepDLNM.eta.ul <- as.vector(datDLNM.modepDLNM.eta.est + 1.96*datDLNM.modepDLNM.pred.object$se.fit)
          datDLNM.modepDLNM.eta.ll <- as.vector(datDLNM.modepDLNM.eta.est - 1.96*datDLNM.modepDLNM.pred.object$se.fit)
          
        }
        
        cat("finish: ", lst, "\n")
        #### RETURN ###########
        out = list(para = c(lst),
                    datGAM = list(insample = list(true.eta = eta.gam,
                                                  modDLNM = list(funs = datGAM.modDLNM.summ$est,
                                                                 est = datGAM.modDLNM.eta.est,
                                                                 ll = datGAM.modDLNM.eta.ll,
                                                                 ul = datGAM.modDLNM.eta.ul),
                                                  modGAM01 = list(est = datGAM.modeGAM01.eta.est,
                                                                  ll = datGAM.modeGAM01.eta.ll,
                                                                  ul = datGAM.modeGAM01.eta.ul),
                                                  modGAMavgmaxL = list(est = datGAM.modeGAMavgmaxL.eta.est,
                                                                       ll = datGAM.modeGAMavgmaxL.eta.ll,
                                                                       ul = datGAM.modeGAMavgmaxL.eta.ul),
                                                  modGAMavg7 = list(est = datGAM.modeGAMavg7.eta.est,
                                                                    ll = datGAM.modeGAMavg7.eta.ll,
                                                                    ul = datGAM.modeGAMavg7.eta.ul),
                                                  modGAM0 = list(est = datGAM.modeGAM0.eta.est,
                                                                 ll = datGAM.modeGAM0.eta.ll,
                                                                 ul = datGAM.modeGAM0.eta.ul),
                                                  modpDLNM = list(est = datGAM.modepDLNM.eta.est,
                                                                  ll = datGAM.modepDLNM.eta.ll,
                                                                  ul = datGAM.modepDLNM.eta.ul)),
                                  outsample = list(true.eta = eta.gam.test,
                                                   modDLNM = list(est = datGAM.modDLNM.eta.test.est,
                                                                  ll = datGAM.modDLNM.eta.test.ll,
                                                                  ul = datGAM.modDLNM.eta.test.ul),
                                                   modGAM01 = list(est = datGAM.modeGAM01.eta.test.est,
                                                                   ll = datGAM.modeGAM01.eta.test.ll,
                                                                   ul = datGAM.modeGAM01.eta.test.ul),
                                                   modGAMavgmaxL = list(est = datGAM.modeGAMavgmaxL.eta.test.est,
                                                                        ll = datGAM.modeGAMavgmaxL.eta.test.ll,
                                                                        ul = datGAM.modeGAMavgmaxL.eta.test.ul),
                                                   modGAMavg7 = list(est = datGAM.modeGAMavg7.eta.test.est,
                                                                     ll = datGAM.modeGAMavg7.eta.test.ll,
                                                                     ul = datGAM.modeGAMavg7.eta.test.ul),
                                                   modGAM0 = list(est = datGAM.modeGAM0.eta.test.est,
                                                                  ll = datGAM.modeGAM0.eta.test.ll,
                                                                  ul = datGAM.modeGAM0.eta.test.ul))))
                                  
                                  
          if(avglag == 0) {
            out$datDLNM = list(insample = list(true.eta = eta.DLNM,
                                               modDLNM = list(funs = datDLNM.modDLNM.summ$est,
                                                              est = datDLNM.modDLNM.eta.est,
                                                              ll = datDLNM.modDLNM.eta.ll,
                                                              ul = datDLNM.modDLNM.eta.ul),
                                               modGAM01 = list(est = datDLNM.modeGAM01.eta.est,
                                                               ll = datDLNM.modeGAM01.eta.ll,
                                                               ul = datDLNM.modeGAM01.eta.ul),
                                               modGAMavgmaxL = list(est = datDLNM.modeGAMavgmaxL.eta.est,
                                                                    ll = datDLNM.modeGAMavgmaxL.eta.ll,
                                                                    ul = datDLNM.modeGAMavgmaxL.eta.ul),
                                               modGAMavg7 = list(est = datDLNM.modeGAMavg7.eta.est,
                                                                 ll = datDLNM.modeGAMavg7.eta.ll,
                                                                 ul = datDLNM.modeGAMavg7.eta.ul),
                                               modGAM0 = list(est = datDLNM.modeGAM0.eta.est,
                                                              ll = datDLNM.modeGAM0.eta.ll,
                                                              ul = datDLNM.modeGAM0.eta.ul),
                                               modpDLNM = list(est = datDLNM.modepDLNM.eta.est,
                                                               ll = datDLNM.modepDLNM.eta.ll,
                                                               ul = datDLNM.modepDLNM.eta.ul)),
                               outsample = list(true.eta = eta.DLNM.test,
                                                modDLNM = list(est = datDLNM.modDLNM.eta.test.est,
                                                               ll = datDLNM.modDLNM.eta.test.ll,
                                                               ul = datDLNM.modDLNM.eta.test.ul),
                                                modGAM01 = list(est = datDLNM.modeGAM01.eta.test.est,
                                                                ll = datDLNM.modeGAM01.eta.test.ll,
                                                                ul = datDLNM.modeGAM01.eta.test.ul),
                                                modGAMavgmaxL = list(est = datDLNM.modeGAMavgmaxL.eta.test.est,
                                                                     ll = datDLNM.modeGAMavgmaxL.eta.test.ll,
                                                                     ul = datDLNM.modeGAMavgmaxL.eta.test.ul),
                                                modGAMavg7 = list(est = datDLNM.modeGAMavg7.eta.test.est,
                                                                  ll = datDLNM.modeGAMavg7.eta.test.ll,
                                                                  ul = datDLNM.modeGAMavg7.eta.test.ul),
                                                modGAM0 = list(est = datDLNM.modeGAM0.eta.test.est,
                                                               ll = datDLNM.modeGAM0.eta.test.ll,
                                                               ul = datDLNM.modeGAM0.eta.test.ul)))
          }
          return(out)
      },
      error = function(e){
        cat("error in ", lst, ".\n")
        warning(e)
      }
    )
  }
  
  
  
  options(mc.cores = parallel::detectCores()-1)
  
  simresults <- mclapply(para_run_list, run_lst)
  save.image(file = file.path(resultpath, paste0("02-simulation-", as.character(avglag), "-results-July30.RData")))
  rm("simresults")
  cat("Finish avglag: ", avglag, ". \n")
  cat("Saved in", file.path(resultpath, paste0("02-simulation-avglag-", as.character(avglag), "-results-July30.RData")), ".\n")
}  
cat("DONE!")
