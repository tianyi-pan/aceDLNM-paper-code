###### Simulation A ###############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Section 4.1 Simulation A: Absolute Performance
## This code is used to generate the tables/figures in the paper. The results are saved in the "results" folder by running the script "01-simulation-A.R"
## Tianyi Pan
## 2025
###################################

## load packages
library(dplyr)
library(aceDLNM)
library(ggplot2)
library(ggpubr)
library(knitr)

## paths
resultpath <- "results"
figurepath <- "figures"
if (!dir.exists(figurepath)) dir.create(figurepath, recursive = TRUE)


### load data #####
load(file = file.path(resultpath, "01-simulation-results-1.RData"))
load(file = file.path(resultpath, "01-simulation-results-2.RData"))
load(file = file.path(resultpath, "01-simulation-results-3.RData"))
load(file = file.path(resultpath, "01-simulation-results-4.RData"))
load(file = file.path(resultpath, "01-simulation-results-5.RData"))
load(file = file.path(resultpath, "01-simulation-results-6.RData"))
load(file = file.path(resultpath, "01-simulation-results-7.RData"))
load(file = file.path(resultpath, "01-simulation-results-8.RData"))
load(file = file.path(resultpath, "01-simulation-results-9.RData"))
load(file = file.path(resultpath, "01-simulation-results-10.RData"))
load(file = file.path(resultpath, "01-simulation-results-11.RData"))
load(file = file.path(resultpath, "01-simulation-results-12.RData"))
load(file = file.path(resultpath, "01-simulation-results-13.RData"))
load(file = file.path(resultpath, "01-simulation-results-14.RData"))
load(file = file.path(resultpath, "01-simulation-results-15.RData"))
load(file = file.path(resultpath, "01-simulation-results-16.RData"))


simresults <- c(simresults1, simresults2,simresults3,simresults4,
                simresults5, simresults6,simresults7,simresults8,
                simresults9, simresults10,simresults11,simresults12,
                simresults13, simresults14,simresults15,simresults16)


### TABLES ######
# report.delta <- TRUE ## uncomment this line for Delta method
report.delta <- FALSE ## uncomment this line for sampling method

## read data
sumresults <- lapply(simresults, function(results){
  wlresults <- data.frame()
  fEresults <- data.frame()
  gtresults <- data.frame()
  thetaresults <- data.frame()
  if(length(results) > 1){
    if(report.delta) result <- results$model.delta
    else result <- results$model
    para <- results$para
    convergence <- results$model$convergence == 0
    if(!is.null(result$wl)){
      result$wl$Nt = as.numeric(para[1])
      result$wl$wltype = para[2]
      result$wl$fEtype = para[3]
      result$wl$iter = as.numeric(para[4])
      result$wl$time = result$opttime
      result$wl$convergence = convergence
      result$wl$par = results$model$opt$par[3]
      wlresults <- rbind(wlresults, result$wl)
      
      result$fE$Nt = as.numeric(para[1])
      result$fE$wltype = para[2]
      result$fE$fEtype = para[3]
      result$fE$iter = as.numeric(para[4])
      result$fE$time = result$opttime
      result$fE$convergence = convergence
      result$fE$par = results$model$opt$par[2]
      fEresults <- rbind(fEresults, result$fE)

      result$gt$Nt = as.numeric(para[1])
      result$gt$wltype = para[2]
      result$gt$fEtype = para[3]
      result$gt$iter = as.numeric(para[4])
      result$gt$time = result$opttime
      result$gt$true <- result$gt$true - mean(result$gt$true) # center
      result$gt$mode <- result$gt$mode - mean(result$gt$mode) # center
      result$gt$ll <- result$gt$ll - mean(result$gt$mode) # center
      result$gt$ul <- result$gt$ul - mean(result$gt$mode) # center
      result$gt$par = results$model$opt$par[4]
      result$gt$convergence = convergence
      gtresults <- rbind(gtresults, result$gt)
      
      thetaresults <- rbind(thetaresults,
                            data.frame(Nt = as.numeric(para[1]),
                                       wltype = para[2],
                                       fEtype = para[3],
                                       iter = as.numeric(para[4]),
                                       theta = exp(results$model$opt$par[1])))
    }
  }
  return(list(wlresults, fEresults, gtresults, thetaresults))
}
)

wlresults <- lapply(sumresults, "[[", 1)
wlresults <- data.table::rbindlist(wlresults)

fEresults <- lapply(sumresults, "[[", 2)
fEresults <- data.table::rbindlist(fEresults)

gtresults <- lapply(sumresults, "[[", 3)
gtresults <- data.table::rbindlist(gtresults)

thetaresults <- lapply(sumresults, "[[", 4)
thetaresults <- data.table::rbindlist(thetaresults)


wlsumm <- mutate(wlresults, cvg = as.numeric((ul >= true) & (ll <= true))) %>% 
  group_by(Nt, wltype, fEtype, iter) %>%
  summarise(RMSE = sqrt(mean((mode-true)^2)),
            cvg = mean(cvg),
            length = mean(ul-ll)) %>% 
  group_by(Nt, wltype, fEtype) %>% 
  summarise(RMSE.mean = mean(RMSE),
            RMSE.MCSE = sd(RMSE)/sqrt(n()),
            cvg.mean = mean(cvg),
            cvg.MCSE = sd(cvg)/sqrt(n()),
            length.mean = mean(length),
            length.MCSE = sd(length)/sqrt(n()),
            n = n())

fEsumm <- mutate(fEresults, cvg = as.numeric((ul >= true) & (ll <= true))) %>% 
  group_by(Nt, wltype, fEtype, iter) %>%
  summarise(RMSE = sqrt(mean((mode-true)^2)),
            cvg = mean(cvg),
            length = mean(ul-ll)) %>% 
  group_by(Nt, wltype, fEtype) %>% 
  summarise(RMSE.mean = mean(RMSE),
            RMSE.MCSE = sd(RMSE)/sqrt(n()),
            cvg.mean = mean(cvg),
            cvg.MCSE = sd(cvg)/sqrt(n()),
            length.mean = mean(length),
            length.MCSE = sd(length)/sqrt(n()))


thetasumm <- thetaresults %>%
  group_by(Nt, wltype, fEtype) %>%
  summarise(bias.mean = mean(theta-8),
            bias.MCSE = sd(theta)/sqrt(n()),
            RMSE.mean = sqrt(mean((theta-8)^2)),
            RMSE.MCSE = sd((theta-8)^2)/(2*sqrt(mean((theta-8)^2))*sqrt(n())),
            n = n())

gtsumm <- mutate(gtresults, cvg = as.numeric((ul >= true) & (ll <= true))) %>%
  group_by(Nt, wltype, fEtype, iter) %>%
  summarise(RMSE = sqrt(mean((mode-true)^2)),
            cvg = mean(cvg),
            length = mean(ul-ll)) %>% 
  group_by(Nt, wltype, fEtype) %>% 
  summarise(RMSE.mean = mean(RMSE),
            RMSE.MCSE = sd(RMSE)/sqrt(n()),
            cvg.mean = mean(cvg),
            cvg.MCSE = sd(cvg)/sqrt(n()),
            length.mean = mean(length),
            length.MCSE = sd(length)/sqrt(n()))


### report TABLES ###### 

## w and f
# choose one of them
wltype. <- "type1"
wltype. <- "type2"
wltype. <- "type3"

kable(filter(wlsumm,
       wltype == wltype., fEtype == "cubic") %>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE, n),
      digits = 3)
kable(filter(fEsumm, 
             wltype == wltype., fEtype == "cubic")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)

kable(filter(wlsumm, 
             wltype == wltype., fEtype == "quadratic")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)
kable(filter(fEsumm, 
             wltype == wltype., fEtype == "quadratic")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)

kable(filter(wlsumm, 
             wltype == wltype., fEtype == "linear")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)
kable(filter(fEsumm, 
             wltype == wltype., fEtype == "linear")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)


## gt
# choose one of them
wltype. <- "type1"
wltype. <- "type2"
wltype. <- "type3"

kable(filter(gtsumm, 
             wltype == wltype., fEtype == "cubic")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)
kable(filter(gtsumm, 
             wltype == wltype., fEtype == "quadratic")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)
kable(filter(gtsumm, 
             wltype == wltype., fEtype == "linear")%>% 
        select(Nt, RMSE.mean, cvg.mean, length.mean, RMSE.MCSE, cvg.MCSE, length.MCSE),
      digits = 3)

## theta
# choose one of them
wltype. <- "type1"
wltype. <- "type2"
wltype. <- "type3"

kable(filter(thetasumm,
             wltype == wltype., fEtype == "cubic") %>% 
        select(Nt, bias.mean, RMSE.mean, bias.MCSE, RMSE.MCSE),
      digits = 3)
kable(filter(thetasumm,
             wltype == wltype., fEtype == "quadratic") %>% 
        select(Nt, bias.mean, RMSE.mean, bias.MCSE, RMSE.MCSE),
      digits = 3)
kable(filter(thetasumm,
             wltype == wltype., fEtype == "linear") %>% 
        select(Nt, bias.mean, RMSE.mean, bias.MCSE, RMSE.MCSE),
      digits = 3)



### plot for point estimates ######
## ggplot settings
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25 * 1.8
GGPLOTSTSIZE <- 24 * 1.8
GGPLOTLEGENDSIZE <- 23 * 1.8
PLOTWIDTH <- 12
PLOTHEIGHT <- 9


for (wltype. in c("type1", "type2", "type3")) {
  for (fEtype. in c("cubic", "quadratic", "linear")){
    for (Nt. in c(1000, 2000)) {
      cat(wltype.)
      wlname <- switch(wltype., 
                      "type1" = "Type (i)",
                      "type2" = "Type (ii)",
                      "type3" = "Type (iii)")
      wllim <- switch(wltype., 
                      "type1" = c(-0.1, 0.45),
                      "type2" = c(-0.1, 0.45),
                      "type3" = c(-0.1, 1))
      fEname <- switch(fEtype.,
                       "cubic" = "Type (i)",
                       "quadratic" = "Type (ii)",
                       "linear" = "Type (iii)")
      fElim <- switch(fEtype., 
                      "cubic" = c(-2, 4),
                      "quadratic" = c(-2, 4.5),
                      "linear" = c(-2, 2))
                      
      wlresult_draw1 <- filter(wlresults,
                               Nt == Nt., 
                               wltype == wltype., 
                               fEtype == fEtype.)
      
      iter.sample <- sort(unique(wlresult_draw1$iter))[1:100]

      wlresult_draw1 <- filter(wlresult_draw1,
                               iter %in% iter.sample)

      wlresult_draw1_mean <- filter(wlresults,
                               Nt == Nt., 
                               wltype == wltype., 
                               fEtype == fEtype.) %>% 
                              group_by(l) %>% 
                              summarise(mode = mean(mode),
                                        true = mean(true))
        
      pw1 <- ggplot()
      pw1 <- pw1 + geom_line(data = wlresult_draw1, aes(x = l, y = mode, group = iter), linewidth = 1.5, alpha = 0.08, color = "#2c7bb6") + 
        geom_line(data = wlresult_draw1_mean, aes(x = l, y = true), color = "black", linewidth = 2, alpha = 0.5)
      pw1 <- pw1 + geom_line(data = wlresult_draw1_mean, aes(x = l, y = mode), linetype = "dashed", linewidth = 2, color = "blue")
      pw1 <- pw1 + xlab("Lag") + ylab("Weight")
      pw1 <- pw1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
                         legend.text=element_text(size=GGPLOTLEGENDSIZE),
                         plot.subtitle=element_text(size=GGPLOTSTSIZE),
                         legend.spacing.y = unit(0.3, 'cm') ) + 
             scale_x_continuous(breaks=seq(0,14,1)) + 
             ylim(wllim)

      pw1 <- pw1 + ggtitle(paste0(wlname," w")) + geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7)
      
      fEresult_draw1 <- filter(fEresults,
                               Nt == Nt., 
                               wltype == wltype., 
                               fEtype == fEtype.,
                               iter %in% iter.sample)
      fEresult_draw1_mean <- filter(fEresults,
                               Nt == Nt., 
                               wltype == wltype., 
                               fEtype == fEtype.) %>% 
                              group_by(E) %>% 
                              summarise(mode = mean(mode),
                                        true = mean(true))
      pE1 <- ggplot()
      pE1 <- pE1 + geom_line(data = fEresult_draw1, aes(x = E, y = mode, group = iter), linewidth = 1.5, alpha = 0.08, color = "#2c7bb6") + 
        geom_line(data = fEresult_draw1_mean, aes(x = E, y = true), color = "black", linewidth = 2, alpha = 0.5)
      pE1 <- pE1 + geom_line(data = fEresult_draw1_mean, aes(x = E, y = mode), linetype = "dashed", linewidth = 2, color = "blue")
      pE1 <- pE1 + xlab("Adaptive Cumulative Exposure") + ylab("ACERF")
      pE1 <- pE1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
                         legend.text=element_text(size=GGPLOTLEGENDSIZE),
                         plot.subtitle=element_text(size=GGPLOTSTSIZE),
                         legend.spacing.y = unit(0.3, 'cm') ) + 
              ylim(fElim)                         
      pE1 <- pE1 + ggtitle(paste0(fEname," f")) 
      
      gtresult_draw1 <- filter(gtresults,
                               Nt == Nt., 
                               wltype == wltype., 
                               fEtype == fEtype.,
                               iter %in% iter.sample)
      gtresult_draw1_mean <- filter(gtresults,
                               Nt == Nt., 
                               wltype == wltype., 
                               fEtype == fEtype.) %>% 
                              group_by(x) %>% 
                              summarise(mode = mean(mode),
                                        true = mean(true))
      pg1 <- ggplot()
      pg1 <- pg1 + geom_line(data = gtresult_draw1, aes(x = x, y = mode, group = iter), linewidth = 1.5, alpha = 0.08, color = "#2c7bb6") + 
        geom_line(data = gtresult_draw1_mean, aes(x = x, y = true), color = "black", linewidth = 2, alpha = 0.5)
      pg1 <- pg1 + geom_line(data = gtresult_draw1_mean, aes(x = x, y = mode), linetype = "dashed", linewidth = 2, color = "blue")
      pg1 <- pg1 + xlab("t") + ylab("Time Trend")
      pg1 <- pg1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
                         legend.text=element_text(size=GGPLOTLEGENDSIZE),
                         plot.subtitle=element_text(size=GGPLOTSTSIZE),
                         legend.spacing.y = unit(0.3, 'cm') )
      pg1 <- pg1 + ggtitle("h") + ylim(c(-0.6, 0.6))

      
      p <- ggarrange(pw1, 
                     pE1, 
                     pg1, 
                     ncol=3, nrow=1)
      ggsave(p, file = file.path(figurepath, paste0("sim1-wl-", wltype., "-fE-", fEtype., "-Nt-", Nt., ".pdf")), width = 3*PLOTWIDTH, height = PLOTHEIGHT)
    }
  }
}  

## for Figure 1 in the main text 
for(Nt. in c(1000, 2000)) {
  pw.list <- list()
  pf.list <- list()
  for (wltype. in c("type1", "type2", "type3")) {
        cat(Nt.)      
        wlname <- switch(wltype., 
                        "type1" = "Type (i)",
                        "type2" = "Type (ii)",
                        "type3" = "Type (iii)")
        wllim <- c(-0.1, 1)
        fEtype. <- switch(wltype., 
                        "type1" = "cubic",
                        "type2" = "quadratic",
                        "type3" = "linear")
        fEname <- switch(fEtype.,
                        "cubic" = "Type (i)",
                        "quadratic" = "Type (ii)",
                        "linear" = "Type (iii)")
        fElim <- c(-2, 4.5)
                        
        wlresult_draw1 <- filter(wlresults,
                                Nt == Nt., 
                                wltype == wltype., 
                                fEtype == fEtype.)
        
        iter.sample <- sort(unique(wlresult_draw1$iter))[1:100]

        wlresult_draw1 <- filter(wlresult_draw1,
                                iter %in% iter.sample)

        wlresult_draw1_mean <- filter(wlresults,
                                Nt == Nt., 
                                wltype == wltype., 
                                fEtype == fEtype.) %>% 
                                group_by(l) %>% 
                                summarise(mode = mean(mode),
                                          true = mean(true))
          
        pw1 <- ggplot()
        pw1 <- pw1 + geom_line(data = wlresult_draw1, aes(x = l, y = mode, group = iter), linewidth = 1.5, alpha = 0.08, color = "#2c7bb6") + 
          geom_line(data = wlresult_draw1_mean, aes(x = l, y = true), color = "black", linewidth = 2, alpha = 0.5)
        pw1 <- pw1 + geom_line(data = wlresult_draw1_mean, aes(x = l, y = mode), linetype = "dashed", linewidth = 2, color = "blue")
        pw1 <- pw1 + xlab("Lag") + ylab("Weight")
        pw1 <- pw1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
                          legend.text=element_text(size=GGPLOTLEGENDSIZE),
                          plot.subtitle=element_text(size=GGPLOTSTSIZE),
                          legend.spacing.y = unit(0.3, 'cm') ) + 
              scale_x_continuous(breaks=seq(0,14,1)) + 
              ylim(wllim)

        pw1 <- pw1 + ggtitle(paste0(wlname," w")) + geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7)
        pw.list <- c(pw.list, list(pw1))

        fEresult_draw1 <- filter(fEresults,
                                Nt == Nt., 
                                wltype == wltype., 
                                fEtype == fEtype.,
                                iter %in% iter.sample)
        fEresult_draw1_mean <- filter(fEresults,
                                Nt == Nt., 
                                wltype == wltype., 
                                fEtype == fEtype.) %>% 
                                group_by(E) %>% 
                                summarise(mode = mean(mode),
                                          true = mean(true))
        pE1 <- ggplot()
        pE1 <- pE1 + geom_line(data = fEresult_draw1, aes(x = E, y = mode, group = iter), linewidth = 1.5, alpha = 0.08, color = "#2c7bb6") + 
          geom_line(data = fEresult_draw1_mean, aes(x = E, y = true), color = "black", linewidth = 2, alpha = 0.5)
        pE1 <- pE1 + geom_line(data = fEresult_draw1_mean, aes(x = E, y = mode), linetype = "dashed", linewidth = 2, color = "blue")
        pE1 <- pE1 + xlab("Adaptive Cumulative Exposure") + ylab("ACERF")
        pE1 <- pE1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
                          legend.text=element_text(size=GGPLOTLEGENDSIZE),
                          plot.subtitle=element_text(size=GGPLOTSTSIZE),
                          legend.spacing.y = unit(0.3, 'cm') ) + 
                ylim(fElim)                         
        pE1 <- pE1 + ggtitle(paste0(fEname," f")) # + geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7)
        pf.list <- c(pf.list, list(pE1))

  }  
  p <- ggarrange(pw.list[[1]], pw.list[[2]] + ylab(" "), pw.list[[3]] + ylab(" "), 
                pf.list[[1]], pf.list[[2]] + ylab(" "), pf.list[[3]] + ylab(" "),
                ncol=3, nrow=2)
  ggsave(p, file = file.path(figurepath, paste0("sim1-main-", Nt., ".pdf")), width = 3*PLOTWIDTH, height = 2*PLOTHEIGHT)
}



### Web Appendix H.1 ##############
## Web Figure 2 a #######
library(forecast); library(fpp3)

data("PM25Waterloo")
data <- PM25.waterloo[1:(1014), ]

tsdata <- tsibble(
  date = data$date,
  PM25 = data$PM25,
  index = date
)

tsdata |>
  mutate(date = ymd(date)) |>
  as_tsibble(index = date)

tsdata |>
  features(PM25, feat_stl)


p <- tsdata |>
  autoplot() + theme_bw() + xlab("date")

ggsave(p, file = "figures/WaterlooPlot.pdf", width = 5, height = 3)

## Web Figure 2 b #######
library(lubridate)
library(tidyverse)
library(mgcv)

data("PM25Waterloo")
data <- PM25.waterloo[1:(1014), ]

starting.date <- data$date[1]

data <- mutate(data, date.num = as.numeric(difftime(date, starting.date, units = "days"))) 
# add Day of Week (DOW)
data <- mutate(data, DOW = as.factor(wday(date)))
# add month
data <- mutate(data, month = month(date))
data$logPM25 <- log(data$PM25)

mod <- gam(logPM25 ~ s(date.num, bs = "bs", k = 10) + 
             s(month, bs = "cc", k = 10),
           data = data)

pdf(file = "figures/WaterlooGAM.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(mod, select = 1, xlab = "date", ylab = "h1(date)", xaxt='n')
x1 <- mod$model$date
axis(1, at = pretty(x1), labels = format(starting.date + pretty(x1), "%Y-%m-%d"))

plot(mod, select = 2, xlab = "month", ylab = "h2(month)", xaxt='n')
x2 <- mod$model$month
axis(1, at = 1:12, labels = month.abb)

dev.off()