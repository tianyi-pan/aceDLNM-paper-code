###### Supplementary Simulation: Misspecification Scenario A ############
## Paper: Estimating Associations Between Cumulative Exposure and Health via Generalized Distributed Lag Non-Linear Models using Penalized Splines
## Web Appendix H.4.1 Supplementary Simulation: Misspecification Scenario A
## This code is used to generate the tables/figures in the paper. The results are saved in the "results" folder by running the script "05-simulation-appendix.R"
## Tianyi Pan
## 2025
###################################

## load packages
library(dplyr)
# devtools::install_github("tianyi-pan/aceDLNM")
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
load(file = file.path(resultpath, "05-simulation-appendix-results-1.RData"))
load(file = file.path(resultpath, "05-simulation-appendix-results-2.RData"))
load(file = file.path(resultpath, "05-simulation-appendix-results-3.RData"))
load(file = file.path(resultpath, "05-simulation-appendix-results-4.RData"))

simresults <- c(simresults1, simresults2, simresults3, simresults4)


####### TABLES #############
datC1_results <- lapply(simresults, function(results){

  out <- list(E = list(),
              other = list())

  if(length(results) > 1){
    Nt <- as.numeric(results$para[1])
    wltype <- results$para[2]
    fxtype <- results$para[3]
    iter <- as.numeric(results$para[4])

    ## distributed lag terms
    result_E <- results$datC1$E
    true.eta_E <- result_E$true.eta
    
    out$E$DLNM.E.RMSE <- sqrt(mean((result_E$modDLNM$est - true.eta_E)^2))
    out$E$DLNM.E.Cvg <- mean((result_E$modDLNM$ul >= true.eta_E ) & (result_E$modDLNM$ll <= true.eta_E ))
    

    out$E$pDLNM.E.RMSE <- sqrt(mean((result_E$modpDLNM$est - true.eta_E)^2))
    out$E$pDLNM.E.Cvg <- mean((result_E$modpDLNM$ul >= true.eta_E ) & (result_E$modpDLNM$ll <= true.eta_E ))
    

    out$E <- as.data.frame(out$E)
    out$E$Nt <- Nt
    out$E$wltype <- wltype
    out$E$fxtype <- fxtype
    out$E$iter <- iter
    
    ## other terms
    result_other <- results$datC1$other
    true.eta_other <- result_other$true.eta
    
    out$other$DLNM.other.RMSE <- sqrt(mean((result_other$modDLNM$est - true.eta_other)^2))
    out$other$DLNM.other.Cvg <- mean((result_other$modDLNM$ul >= true.eta_other) & (result_other$modDLNM$ll <= true.eta_other))
    
    out$other$pDLNM.other.RMSE <- sqrt(mean((result_other$modpDLNM$est - true.eta_other)^2))
    out$other$pDLNM.other.Cvg <- mean((result_other$modpDLNM$ul >= true.eta_other) & (result_other$modpDLNM$ll <= true.eta_other))
    

    out$other <- as.data.frame(out$other)
    out$other$Nt <- Nt
    out$other$wltype <- wltype
    out$other$fxtype <- fxtype
    out$other$iter <- iter
    
  }
  return(out)
}
)


datC1_E <- lapply(datC1_results, "[[", "E")
datC1_E <- data.table::rbindlist(datC1_E)

datC1_other <- lapply(datC1_results, "[[", "other")
datC1_other <- data.table::rbindlist(datC1_other)



#### report tables
datC1_E.DLNM.summ <- group_by(datC1_E, Nt, fxtype, wltype) %>% summarise( 
  RMSE.mean = mean(DLNM.E.RMSE),
  RMSE.MCSE = sd(DLNM.E.RMSE) / sqrt(n()),
  Cvg.mean = mean(DLNM.E.Cvg),
  Cvg.MCSE = sd(DLNM.E.Cvg) / sqrt(n()),
  n = n()
)

datC1_E.DLNM.summ %>% filter(wltype == "type1") %>% knitr::kable(digits = 3)


datC1_E.pDLNM.summ <- group_by(datC1_E, Nt, fxtype, wltype) %>% summarise( 
  RMSE.mean = mean(pDLNM.E.RMSE),
  RMSE.MCSE = sd(pDLNM.E.RMSE) / sqrt(n()),
  Cvg.mean = mean(pDLNM.E.Cvg),
  Cvg.MCSE = sd(pDLNM.E.Cvg) / sqrt(n()),
  n = n()
)

datC1_E.pDLNM.summ %>% filter(wltype == "type1") %>% knitr::kable(digits = 3)


filter(datC1_E.summ, Nt == 2000) %>% 
  select(wltype, RMSE.mean, Cvg.mean) %>%  #fxtype
  knitr::kable(digits = 3)


datC1_other.summ <- group_by(datC1_other, Nt, wltype, fxtype) %>% summarise(
  RMSE.mean = mean(DLNM.other.RMSE),
  RMSE.MCSE = sd(DLNM.other.RMSE) / sqrt(n()),
  Cvg.mean = mean(DLNM.other.Cvg),
  Cvg.MCSE = sd(DLNM.other.Cvg) / sqrt(n())
)

knitr::kable(datC1_other.summ,  digits = 3)


########### Plotting the association functions ############
a.list <- sqrt(seq(0,1, by = 0.1))
prop.curv <- a.list^2
prop.curv*100
x.med <- (x.min + x.max)/2

for(a in a.list) {
  
  fx <- function(x) ((x-(x.max+x.min)/2) / ((x.max - x.min)/3.5) + 3)/2 - a*sqrt(0.5)*(x-x.med)^2/((x.max-x.med)^2) + a*sqrt(0.5)
  curve(fx, from = x.min, to = x.max, add = TRUE, col = "red")
}

x.seq <- seq(x.min, x.max, length.out = 500)
df <- lapply(a.list, function(a) {
  fx <- function(x) ((x-(x.max+x.min)/2) / ((x.max - x.min)/3.5) + 3)/2 - a*sqrt(0.5)*(x-x.med)^2/((x.max-x.med)^2) + a*sqrt(0.5)
  return(data.frame(x = x.seq,
                    f = fx(x.seq),
                    J = rep(a^2, 500)))
})
df <- data.table::rbindlist(df)

df$J <- factor(as.character(paste0(df$J*100, "%")), levels = paste0(as.character(seq(1,0,by = -0.1)*100), "%"))

p <- ggplot(data = df, aes(x = x, y = f, group = J, color = J)) + 
  geom_line(linewidth = 2) + 
  scale_color_manual(values =colorRampPalette(c("#08306b", "#bdd7e7"))(11),
                     name = "Relative J(f)") + 
  theme(text = element_text(size = 1*GGPLOTTEXTSIZE),
          legend.text=element_text(size=1.2*GGPLOTLEGENDSIZE),
          plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
          plot.subtitle=element_text(size=0.9*GGPLOTSTSIZE),
          legend.spacing.y = unit(0.3, 'cm'))
ggsave(p, file = "figures/fcurve.pdf", width = PLOTWIDTH, height = 1.2*PLOTHEIGHT)
