library(dplyr)
library(aceDLNM)

resultpath <- "results"
figurepath <- "figures"
if (!dir.exists(figurepath)) dir.create(figurepath, recursive = TRUE)

library(ggplot2)
library(ggpubr)

theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25
GGPLOTSTSIZE <- 24
GGPLOTLEGENDSIZE <- 23

PLOTWIDTH <- 12
PLOTHEIGHT <- 9


### load data #####
load(file = file.path(resultpath, "01-simulation-results-1-July30.RData"))
load(file = file.path(resultpath, "01-simulation-results-2-July30.RData"))
load(file = file.path(resultpath, "01-simulation-results-3-July30.RData"))
load(file = file.path(resultpath, "01-simulation-results-4-July30.RData"))
load(file = file.path(resultpath, "01-simulation-results-5-July30.RData"))
load(file = file.path(resultpath, "01-simulation-results-6-July30.RData"))
load(file = file.path(resultpath, "01-simulation-results-7-July30.RData"))
load(file = file.path(resultpath, "01-simulation-results-8-July30.RData"))


simresults <- c(simresults1, simresults2,simresults3,simresults4,
                simresults5, simresults6,simresults7,simresults8)


### TABLES ######
# report.delta <- TRUE
report.delta <- FALSE

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
  group_by(Nt, wltype, fEtype) %>%
  summarise(bias.mean = mean(mode-true),
            RMSE = sqrt(mean((mode-true)^2)),
            cvg.mean = mean(cvg),
            cvg.MCSE = sqrt(mean(cvg)*(1-mean(cvg))/(n()/500)),
            length = mean(ul-ll),
            time = mean(time),
            convergence = mean(convergence),
            par = mean(par),
            n = n())




fEsumm <- mutate(fEresults, cvg = as.numeric((ul >= true) & (ll <= true))) %>% 
  # filter(convergence == 1) %>%
  group_by(Nt, wltype, fEtype) %>%
  summarise(bias.mean = mean(mode-true), 
            RMSE = sqrt(mean((mode-true)^2)),
            cvg.mean = mean(cvg),
            cvg.MCSE = sqrt(mean(cvg)*(1-mean(cvg))/(n()/500)),
            length = mean(ul-ll),
            convergence = mean(convergence),
            par = mean(par),
            time = mean(time))


thetasumm <- thetaresults %>%
  group_by(Nt, wltype, fEtype) %>%
  summarise(bias.mean = mean(theta-8),
            RMSE = sqrt(mean((theta-8)^2)),
            n = n())


gtsumm <- mutate(gtresults, cvg = as.numeric((ul >= true) & (ll <= true))) %>%
  group_by(Nt, wltype, fEtype) %>%
  summarise(bias.mean = mean(mode-true),
            RMSE = sqrt(mean((mode-true)^2)),
            cvg.mean = mean(cvg),
            cvg.MCSE = sqrt(mean(cvg)*(1-mean(cvg))/n()),
            convergence = mean(convergence),
            par = mean(par),
            length = mean(ul-ll))

### TABLES
library(knitr)
wltype. <- "type3"
kable(filter(wlsumm,
       wltype == wltype., fEtype == "cubic") %>% 
        select(Nt, RMSE, cvg.mean, length),
      digits = 3, "latex")
kable(filter(fEsumm, 
             wltype == wltype., fEtype == "cubic")%>% 
        select(Nt, RMSE, cvg.mean, length, time) %>% 
        mutate(time = round(time)),
      digits = 3, "latex")

kable(filter(wlsumm, 
             wltype == wltype., fEtype == "quadratic")%>% 
        select(Nt, RMSE, cvg.mean, length),
      digits = 3, "latex")
kable(filter(fEsumm, 
             wltype == wltype., fEtype == "quadratic")%>% 
        select(Nt, RMSE, cvg.mean, length, time) %>% 
        mutate(time = round(time)),
      digits = 3, "latex")

kable(filter(wlsumm, 
             wltype == wltype., fEtype == "linear")%>% 
        select(Nt, RMSE, cvg.mean, length),
      digits = 3, "latex")
kable(filter(fEsumm, 
             wltype == wltype., fEtype == "linear")%>% 
        select(Nt, RMSE, cvg.mean, length, time) %>% 
        mutate(time = round(time)),
      digits = 3, "latex")


## gt
wltype. <- "type3"
kable(filter(gtsumm, 
             wltype == wltype., fEtype == "cubic")%>% 
        select(Nt, RMSE, cvg.mean, length),
      digits = 3, "latex")
kable(filter(gtsumm, 
             wltype == wltype., fEtype == "quadratic")%>% 
        select(Nt, RMSE, cvg.mean, length),
      digits = 3, "latex")
kable(filter(gtsumm, 
             wltype == wltype., fEtype == "linear")%>% 
        select(Nt, RMSE, cvg.mean, length),
      digits = 3, "latex")

## theta
kable(filter(thetasumm,
             wltype == wltype., fEtype == "cubic") %>% 
        select(Nt, bias.mean, RMSE),
      digits = 3, "latex")
kable(filter(thetasumm,
             wltype == wltype., fEtype == "quadratic") %>% 
        select(Nt, bias.mean, RMSE),
      digits = 3, "latex")
kable(filter(thetasumm,
             wltype == wltype., fEtype == "linear") %>% 
        select(Nt, bias.mean, RMSE),
      digits = 3, "latex")



### plot for point estimates ######
library(ggplot2)
library(ggpubr)

theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25 * 1.8
GGPLOTSTSIZE <- 24 * 1.8
GGPLOTLEGENDSIZE <- 23 * 1.8

PLOTWIDTH <- 12
PLOTHEIGHT <- 9


iter.sample <- 1:100

for (wltype. in c("type1", "type2", "type3")) {
  for (fEtype. in c("cubic", "quadratic", "linear")){
    for (Nt. in c(1000, 2000)) {
      cat(wltype.)
      wlname <- switch(wltype., 
                      "type1" = "Type 1",
                      "type2" = "Type 2",
                      "type3" = "Type 3")
      wllim <- switch(wltype., 
                      "type1" = c(-0.1, 0.45),
                      "type2" = c(-0.1, 0.45),
                      "type3" = c(-0.1, 1))
      fEname <- switch(fEtype.,
                       "cubic" = "Type 1",
                       "quadratic" = "Type 2",
                       "linear" = "Type 3")
      fElim <- switch(fEtype., 
                      "cubic" = c(-2, 4),
                      "quadratic" = c(-2, 4.5),
                      "linear" = c(-2, 2))
      wlresult_draw1 <- filter(wlresults,
                               Nt == Nt., 
                               wltype == wltype., 
                               fEtype == fEtype.,
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
      # ggsave(pw1, file = file.path(figurepath, paste0("sim1-w-wl-", wltype., "-fE-", fEtype., "-Nt-", Nt., ".pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)
      
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
      # ggsave(pE1, file = file.path(figurepath, paste0("sim1-f-wl-", wltype., "-fE-", fEtype., "-Nt-", Nt., ".pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)
      
      
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
      # ggsave(pg1, file = file.path(figurepath, paste0("sim1-gt-wl-", wltype., "-fE-", fEtype., "-Nt-", Nt., ".pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)
      
      p <- ggarrange(pw1, 
                     pE1, 
                     pg1, 
                     ncol=3, nrow=1)
      ggsave(p, file = file.path(figurepath, paste0("sim1-wl-", wltype., "-fE-", fEtype., "-Nt-", Nt., ".pdf")), width = 3*PLOTWIDTH, height = PLOTHEIGHT)
    }
  }
}  




### other plots ###########

## 1. Weight functions
maxL <- 14
maxLreal <- 15
wll <- function(l, type){
 switch (type,
         type1 = {
           wlc <- function(l) dnorm(l, mean = 3, sd = 3.5)^2
           wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
           wl <- function(l) dnorm(l, mean = 3, sd = 3.5)/wl_de
         },
         type2 = {
           wlc <- function(l) (1/(1+exp(l-8)))^2
           wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
           wl <- function(l) 1/(1+exp(l-8))/wl_de
         },
         type3 = {
           wlc <- function(l) (1/(1+exp(l-1)))^2
           wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
           wl <- function(l) 1/(1+exp(l-1))/wl_de
         }
      )
  return(wl(l))
}

l <- seq(0, 15, length.out = 1000)
wltrue1 <- wll(l, type = "type1")
wltrue2 <- wll(l, type = "type2")
wltrue3 <- wll(l, type = "type3")

wl_df <- data.frame(l = rep(l, 3),
                    wl = c(wltrue1,wltrue2,wltrue3),
                    Type = rep(c("Type1", "Type2", "Type3"), each = 1000))

p <- ggplot(data = wl_df, aes(x = l, y = wl)) +
  geom_line(aes(color = Type), linewidth = 3) + 
  xlab("Lag") + ylab("Weight w(lag)")
p <- p + theme(text = element_text(size = GGPLOTTEXTSIZE),
             legend.text=element_text(size=GGPLOTLEGENDSIZE),
             plot.subtitle=element_text(size=GGPLOTSTSIZE),
             legend.spacing.y = unit(0.9, 'cm') ) +
  guides(color = guide_legend(title = NULL)) + 
  scale_x_continuous(breaks=seq(0,14,1))

p <- p + scale_color_manual(values = c("Type1" = "#d7191c", "Type2" = "#fdae61", "Type3" = "#2c7bb6"))
ggsave(p, file = file.path(figurepath, "true-wl.pdf"), width = PLOTWIDTH, height = PLOTHEIGHT)



## 2. Effect of E
dat <- GenerateData(fEtype = "cubic", wltype = "type1", Nt = 1000,
                    theta = 8, maxL = 14)
Emin <- dat$true.f$Emin
Emax <- dat$true.f$Emax
fEE <- function(E, type){
  switch (type,
    cubic = {
      fEtmp <- function(E) (E-(Emin + (Emax-Emin)*0.3))*(E-(Emin + (Emax-Emin)*0.2))*(E-(Emin + (Emax-Emin)*0.9))
      fE <- function(E) fEtmp(E) / (fEtmp(Emax)/2.5) + 2.5
    },
    linear = {
      fEtmp <- function(E) 0
      fE <- function(E) (E-(Emax+Emin)/2) / ((Emax - Emin)/3.5) + 3
    },
    quadratic = {
      fEtmp <- function(x){25*(dnorm(((2*((x-Emin)/(Emax-Emin)+0.18) - 1.1))))}
      fE <- function(x) -fEtmp(x)+11
    }
  )
  return(fE(E))
}

E <- seq(Emin, Emax, length.out = 1000)
fEtrue1 <- fEE(E, type = "cubic")
fEtrue2 <- fEE(E, type = "quadratic")
fEtrue3 <- fEE(E, type = "linear")


fE_df <- data.frame(E = rep(E, 3),
                    fE = c(fEtrue1,fEtrue2,fEtrue3),
                    Type = rep(c("Type1", "Type2", "Type3"), each = 1000))
fE_df$Type <- factor(fE_df$Type, levels = c("Type1", "Type2", "Type3"))
p <- ggplot(data = fE_df, aes(x = E, y = fE)) +
  geom_line(aes(color = Type), linewidth = 3) + 
  xlab("ACE") + ylab("ACERF f(E)")
p <- p + theme(text = element_text(size = GGPLOTTEXTSIZE),
               legend.text=element_text(size=GGPLOTLEGENDSIZE),
               plot.subtitle=element_text(size=GGPLOTSTSIZE),
               legend.spacing.y = unit(0.9, 'cm') ) +
  guides(color = guide_legend(title = NULL))

p <- p + scale_color_manual(values = c("Type1" = "#d7191c", "Type2" = "#fdae61", "Type3" = "#2c7bb6"))


ggsave(p, file = file.path(figurepath, "true-fE.pdf"), width = PLOTWIDTH, height = PLOTHEIGHT)

## h(t)

gt <- function(x) 0.5*sin(x/150)
gt_df <- data.frame(t = 1:2014,
                    gt = gt(1:2014))

p <- ggplot(data = gt_df, aes(x = t, y = gt)) +
  geom_line(linewidth = 3) + 
  xlab("t") + ylab("Time Trend h(t)")
p <- p + theme(text = element_text(size = GGPLOTTEXTSIZE),
               legend.text=element_text(size=GGPLOTLEGENDSIZE),
               plot.subtitle=element_text(size=GGPLOTSTSIZE),
               legend.spacing.y = unit(0.9, 'cm') ) +
  guides(color = guide_legend(title = NULL))

p <- p + geom_vline(xintercept = 1014, linetype = "dashed", linewidth = 1.5, color = "grey30")

ggsave(p, file = file.path(figurepath, "true-gt.pdf"), width = PLOTWIDTH, height = PLOTHEIGHT)

