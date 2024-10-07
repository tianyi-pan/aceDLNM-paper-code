library(dplyr)
library(aceDLNM)
library(ggplot2)

theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 25
GGPLOTSTSIZE <- 24
GGPLOTLEGENDSIZE <- 23

PLOTWIDTH <- 12
PLOTHEIGHT <- 9


resultpath <- "results"
figurepath <- "figures"
if (!dir.exists(figurepath)) dir.create(figurepath, recursive = TRUE)

## load data
# choose one of them
avglag <- 0
avglag <- 1
avglag <- 7
avglag <- 14

load(file = file.path(resultpath, paste0("02-simulation-", 
                                         as.character(avglag), 
                                         "-results-July30.RData")))


### TABLES ##########
datGAM_results <- lapply(simresults, function(results){

  out <- list(insample = list(),
              outsample = list())

  if(length(results) > 1){
    Nt <- as.numeric(results$para[1])
    iter <- as.numeric(results$para[4])
    
    ## insample
    result_insample <- results$datGAM$insample
    true.eta_insample <- result_insample$true.eta
    
    
    out$insample$DLNM.insample.RMSE <- sqrt(mean((result_insample$modDLNM$est - true.eta_insample)^2))
    out$insample$DLNM.insample.Cvg <- mean((result_insample$modDLNM$ul >= true.eta_insample ) & (result_insample$modDLNM$ll <= true.eta_insample ))
    
    out$insample$GAM01.insample.RMSE <- sqrt(mean((result_insample$modGAM01$est - true.eta_insample)^2))
    out$insample$GAM01.insample.Cvg <- mean((result_insample$modGAM01$ul >= true.eta_insample ) & (result_insample$modGAM01$ll <= true.eta_insample ))
    
    out$insample$GAM0.insample.RMSE <- sqrt(mean((result_insample$modGAM0$est - true.eta_insample)^2))
    out$insample$GAM0.insample.Cvg <- mean((result_insample$modGAM0$ul >= true.eta_insample ) & (result_insample$modGAM0$ll <= true.eta_insample ))
    
    
    out$insample$GAMavgmaxL.insample.RMSE <- sqrt(mean((result_insample$modGAMavgmaxL$est - true.eta_insample)^2))
    out$insample$GAMavgmaxL.insample.Cvg <- mean((result_insample$modGAMavgmaxL$ul >= true.eta_insample ) & (result_insample$modGAMavgmaxL$ll <= true.eta_insample ))
    
    out$insample$GAMavg7.insample.RMSE <- sqrt(mean((result_insample$modGAMavg7$est - true.eta_insample)^2))
    out$insample$GAMavg7.insample.Cvg <- mean((result_insample$modGAMavg7$ul >= true.eta_insample ) & (result_insample$modGAMavg7$ll <= true.eta_insample ))
    
    out$insample$pDLNM.insample.RMSE <- sqrt(mean((result_insample$modpDLNM$est - true.eta_insample)^2))
    out$insample$pDLNM.insample.Cvg <- mean((result_insample$modpDLNM$ul >= true.eta_insample ) & (result_insample$modpDLNM$ll <= true.eta_insample ))
    
    
    out$insample <- as.data.frame(out$insample)
    out$insample$Nt <- Nt
    out$insample$iter <- iter
    
    ## outsample
    result_outsample <- results$datGAM$outsample
    true.eta_outsample <- result_outsample$true.eta
    
    out$outsample$DLNM.outsample.RMSE <- sqrt(mean((result_outsample$modDLNM$est - true.eta_outsample)^2))
    out$outsample$DLNM.outsample.Cvg <- mean((result_outsample$modDLNM$ul >= true.eta_outsample) & (result_outsample$modDLNM$ll <= true.eta_outsample))
    
    out$outsample$GAM01.outsample.RMSE <- sqrt(mean((result_outsample$modGAM01$est - true.eta_outsample)^2))
    out$outsample$GAM01.outsample.Cvg <- mean((result_outsample$modGAM01$ul >= true.eta_outsample) & (result_outsample$modGAM01$ll <= true.eta_outsample))
    
    
    out$outsample$GAM0.outsample.RMSE <- sqrt(mean((result_outsample$modGAM0$est - true.eta_outsample)^2))
    out$outsample$GAM0.outsample.Cvg <- mean((result_outsample$modGAM0$ul >= true.eta_outsample) & (result_outsample$modGAM0$ll <= true.eta_outsample))
    
    
    out$outsample$GAMavgmaxL.outsample.RMSE <- sqrt(mean((result_outsample$modGAMavgmaxL$est - true.eta_outsample)^2))
    out$outsample$GAMavgmaxL.outsample.Cvg <- mean((result_outsample$modGAMavgmaxL$ul >= true.eta_outsample) & (result_outsample$modGAMavgmaxL$ll <= true.eta_outsample))
    
    out$outsample$GAMavg7.outsample.RMSE <- sqrt(mean((result_outsample$modGAMavg7$est - true.eta_outsample)^2))
    out$outsample$GAMavg7.outsample.Cvg <- mean((result_outsample$modGAMavg7$ul >= true.eta_outsample) & (result_outsample$modGAMavg7$ll <= true.eta_outsample))
    
    out$outsample <- as.data.frame(out$outsample)
    out$outsample$Nt <- Nt
    out$outsample$iter <- iter
    
  }
  return(out)
}
)


datGAM_insample <- lapply(datGAM_results, "[[", "insample")
datGAM_insample <- data.table::rbindlist(datGAM_insample)

datGAM_outsample <- lapply(datGAM_results, "[[", "outsample")
datGAM_outsample <- data.table::rbindlist(datGAM_outsample)



if(avglag == 0) {
datDLNM_results <- lapply(simresults, function(results){
  
  out <- list(insample = list(),
              outsample = list())
  
  if(length(results) > 1){
    Nt <- as.numeric(results$para[1])
    iter <- as.numeric(results$para[4])
    
    ## insample
    result_insample <- results$datDLNM$insample
    true.eta_insample <- result_insample$true.eta
    
    
    out$insample$DLNM.insample.RMSE <- sqrt(mean((result_insample$modDLNM$est - true.eta_insample)^2))
    out$insample$DLNM.insample.Cvg <- mean((result_insample$modDLNM$ul >= true.eta_insample) & (result_insample$modDLNM$ll <= true.eta_insample))
    
    out$insample$GAM01.insample.RMSE <- sqrt(mean((result_insample$modGAM01$est - true.eta_insample)^2))
    out$insample$GAM01.insample.Cvg <- mean((result_insample$modGAM01$ul >= true.eta_insample) & (result_insample$modGAM01$ll <= true.eta_insample))
    
    out$insample$GAM0.insample.RMSE <- sqrt(mean((result_insample$modGAM0$est - true.eta_insample)^2))
    out$insample$GAM0.insample.Cvg <- mean((result_insample$modGAM0$ul >= true.eta_insample) & (result_insample$modGAM0$ll <= true.eta_insample))
    
    
    out$insample$GAMavgmaxL.insample.RMSE <- sqrt(mean((result_insample$modGAMavgmaxL$est - true.eta_insample)^2))
    out$insample$GAMavgmaxL.insample.Cvg <- mean((result_insample$modGAMavgmaxL$ul >= true.eta_insample) & (result_insample$modGAMavgmaxL$ll <= true.eta_insample))
    
    out$insample$GAMavg7.insample.RMSE <- sqrt(mean((result_insample$modGAMavg7$est - true.eta_insample)^2))
    out$insample$GAMavg7.insample.Cvg <- mean((result_insample$modGAMavg7$ul >= true.eta_insample) & (result_insample$modGAMavg7$ll <= true.eta_insample))
    
    out$insample$pDLNM.insample.RMSE <- sqrt(mean((result_insample$modpDLNM$est - true.eta_insample)^2))
    out$insample$pDLNM.insample.Cvg <- mean((result_insample$modpDLNM$ul >= true.eta_insample) & (result_insample$modpDLNM$ll <= true.eta_insample))
    
    
    out$insample <- as.data.frame(out$insample)
    out$insample$Nt <- Nt
    out$insample$iter <- iter
    
    ## outsample
    result_outsample <- results$datDLNM$outsample
    true.eta_outsample <- result_outsample$true.eta
    
    out$outsample$DLNM.outsample.RMSE <- sqrt(mean((result_outsample$modDLNM$est - true.eta_outsample)^2))
    out$outsample$DLNM.outsample.Cvg <- mean((result_outsample$modDLNM$ul >= true.eta_outsample) & (result_outsample$modDLNM$ll <= true.eta_outsample))
    

    out$outsample$GAM01.outsample.RMSE <- sqrt(mean((result_outsample$modGAM01$est - true.eta_outsample)^2))
    out$outsample$GAM01.outsample.Cvg <- mean((result_outsample$modGAM01$ul >= true.eta_outsample) & (result_outsample$modGAM01$ll <= true.eta_outsample))
    
    out$outsample$GAM0.outsample.RMSE <- sqrt(mean((result_outsample$modGAM0$est - true.eta_outsample)^2))
    out$outsample$GAM0.outsample.Cvg <- mean((result_outsample$modGAM0$ul >= true.eta_outsample) & (result_outsample$modGAM0$ll <= true.eta_outsample))
    
    
    out$outsample$GAMavgmaxL.outsample.RMSE <- sqrt(mean((result_outsample$modGAMavgmaxL$est - true.eta_outsample)^2))
    out$outsample$GAMavgmaxL.outsample.Cvg <- mean((result_outsample$modGAMavgmaxL$ul >= true.eta_outsample) & (result_outsample$modGAMavgmaxL$ll <= true.eta_outsample))
    
    out$outsample$GAMavg7.outsample.RMSE <- sqrt(mean((result_outsample$modGAMavg7$est - true.eta_outsample)^2))
    out$outsample$GAMavg7.outsample.Cvg <- mean((result_outsample$modGAMavg7$ul >= true.eta_outsample) & (result_outsample$modGAMavg7$ll <= true.eta_outsample))
    
    out$outsample <- as.data.frame(out$outsample)
    out$outsample$Nt <- Nt
    out$outsample$iter <- iter
    
  }
  return(out)
}
)

datDLNM_insample <- lapply(datDLNM_results, "[[", "insample")
datDLNM_insample <- data.table::rbindlist(datDLNM_insample)

datDLNM_outsample <- lapply(datDLNM_results, "[[", "outsample")
datDLNM_outsample <- data.table::rbindlist(datDLNM_outsample)


}

## table 
## datGAM01
n <- nrow(datGAM_insample)
cvg.MCSE <- function(cvg) sqrt(mean(cvg)*(1-mean(cvg))/(n*500))
cvg.MCSE <- function(cvg) sqrt(mean(cvg)*(1-mean(cvg))/(n))

knitr::kable(apply(filter(datGAM_insample, Nt == 2000)%>%select(-iter), 2, mean), digits = 3)
knitr::kable(apply(filter(datGAM_insample, Nt == 2000)%>%select(-iter), 2, cvg.MCSE), digits = 3)


if(avglag == 0) {
  ## datDLNM
  knitr::kable(apply(filter(datDLNM_insample, Nt == 2000)%>%select(-iter), 2, mean), digits = 3)
}




### FIGURES ##########
## figures
datGAM_DLNM <- lapply(simresults, function(results){
  
  wlresults <- data.frame()
  fEresults <- data.frame()
  gtresults <- data.frame()
  
  if(length(results) > 1){
    Nt <- as.numeric(results$para[1])
    iter <- as.numeric(results$para[4])
    
    wlresult <- results$datGAM$insample$modDLNM$funs$wl
    fEresult <- results$datGAM$insample$modDLNM$funs$fE
    gtresult <- results$datGAM$insample$modDLNM$funs$smooth
    
    wlresult$Nt <- Nt
    wlresult$iter <- iter
    
    fEresult$Nt <- Nt
    fEresult$iter <- iter
    
    gtresult$Nt <- Nt
    gtresult$iter <- iter
    
    wlresults <- rbind(wlresults, wlresult)
    fEresults <- rbind(fEresults, fEresult)
    gtresults <- rbind(gtresults, gtresult)
  }
  return(list(wlresults = wlresults,
              fEresults = fEresults,
              gtresults = gtresults))
}
)


wlresults <- lapply(datGAM_DLNM, "[[", 1)
wlresults <- data.table::rbindlist(wlresults)


iter.sample <- sort(unique(wlresults$iter))[1:100] # the first 100 iters
wlresult_draw <- filter(wlresults,
                        Nt == 2000,
                        iter %in% iter.sample)
wlresult_draw_mean <-  group_by(wlresults, l) %>%
  summarise(mode = mean(mode),
            ul = mean(ul),
            ll = mean(ll))
 



wlresult_draw <- mutate(wlresult_draw, draw = ifelse(l < (avglag + 1), "inside", "outside"))
wlresult_draw_mean <- mutate(wlresult_draw_mean, draw = ifelse(l < (avglag + 1), "inside", "outside"))




pw1 <- ggplot(data = wlresult_draw, aes(x = l, y = mode, group = iter, color = draw))
pw1 <- pw1 + geom_line(aes(color = draw), linewidth = 1.5, alpha = 0.08)
pw1 <- pw1 + geom_line(data = wlresult_draw_mean, aes(x = l, y = mode, color = draw), inherit.aes = FALSE,
                       linewidth = 2, linetype = "longdash")

# pw1 <- pw1 + geom_line(aes(x = l, y = ul, color = draw), linetype = "dashed", linewidth = 1.5, alpha = 0.7) +
  # geom_line(aes(x = l, y = ll, color = draw), linetype = "dashed", linewidth = 1.5, alpha = 0.7)
pw1 <- pw1 + xlab("Lag") + ylab("Weight w(l)")
pw1 <- pw1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
                   legend.text=element_text(size=GGPLOTLEGENDSIZE),
                   plot.subtitle=element_text(size=GGPLOTSTSIZE),
                   legend.spacing.y = unit(0.3, 'cm') )

pw1 <- pw1 + geom_vline(xintercept = avglag+1, linewidth = 1,color = "grey",alpha = 0.9) + geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7) +
  # geom_vline(xintercept = 0, linewidth = 1, alpha = 0.7) + geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7) +
  scale_x_continuous(breaks=seq(0,14,1)) +
  annotate('segment', x = 0, xend = avglag+1, y = sqrt(1/(avglag+1)), yend = sqrt(1/(avglag+1)), linewidth = 2, color = "blue") +
  theme(axis.text.x =
          element_text(face = ifelse((0:14) <= avglag, "bold", "plain"),
                        color = ifelse((0:14) <= avglag, "blue", "black"))) +
  scale_color_manual(values = c("inside" = "#2c7bb6", "outside" = "black")) +
  scale_fill_manual(values = c("inside" = "#2c7bb6", "outside" = "black")) +
  theme(legend.position="none")


# pw1
ggsave(pw1, file = file.path(figurepath, paste0("sim2-DLNM-avglag-", avglag, ".pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)

