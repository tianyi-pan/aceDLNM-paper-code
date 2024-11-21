## load packages
library(tidyverse)
library(dplyr)
library(mgcv)
library(lubridate) # for Day Of Week (DOW)

library(aceDLNM)



CD.selected <- c("CD3530", "CD3521","CD3525",
                  "CD5915", "CD4806")

for(CD in CD.selected) {
for(maxL in c(14)) {

ylimu <- 0.3
yliml <- -0.26
E0 <- x0 <- 10

CDname <- switch(CD, "CD3530" = "Waterloo",
                     "CD4806" = "Calgary",
                     "CD3525" = "Hamilton",
                     "CD5915" = "Vancouver",
                     "CD3521" = "Peel")

dataname <- paste0(CD, "_DLNM_morbPulm_modeldata.RData")
loaddata <- file.exists(dataname)

modelname <- file.path(paste0(CD, "-", maxL, "-Pulm-modeldata.RData"))
loadmodel <- file.exists(modelname)

if(loaddata){
    cat("Load data")
    load(file = dataname)
}


if(loadmodel) {
    cat("load model")
    load(file = modelname)
} else {
    cat("Fitting Model for ", CD, " with maxL", maxL, "\n")
    results <- aceDLNM(formula = y~sX(t, x),
                                smooth = ~ s(temp, bs = "bs", k = 10) + 
                                        s(t, bs = "bs", k = 10) + 
                                        s(month.num, bs = "cc", k = 10),
                                fe.varying = ~ DOW,
                                maxL = maxL,
                                pc = NULL,
                                kw = 20, kE = 20, 
                                CI.R = 3000,
                                dat = dat.fit,
                                verbose = TRUE)
    cat(CD, " theta: ", results$point$log_theta, "\n")
    save(results, file = modelname)
}





## draw figures
library(ggplot2)
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 35
GGPLOTSTSIZE <- 40
GGPLOTLEGENDSIZE <- 23
PLOTWIDTH <- 12
PLOTHEIGHT <- 9

res <- residuals(results, plot = FALSE)
qq_dlnm <- res$p.qq + 
    theme_bw() + 
    theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
    labs(subtitle = "ACE-DLNM: QQ-Plot")
res_dlnm <- res$p.res + ylim(-5,5) + 
        theme_bw() + 
        scale_x_continuous(breaks = seq(min(dat.fit$t), max(dat.fit$t), by = round((max(dat.fit$t) - min(dat.fit$t))/30))) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        theme(text = element_text(size = GGPLOTTEXTSIZE),
                legend.text=element_text(size=GGPLOTLEGENDSIZE),
                plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
                plot.subtitle=element_text(size=GGPLOTSTSIZE),
                legend.spacing.y = unit(0.3, 'cm'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        labs(subtitle = "ACE-DLNM: Residual Plot")




# results$point$betaF[1:7]
CI <- apply(results$CI.sample$betaF_sample[,1:7], 2, function(col.) quantile(col., probs = c(0.025, 0.975)))
CI <- rbind(results$point$betaF[1:7], CI)
CI <- t(CI)
knitr::kable(CI, digits = 3)
print(results$point$log_theta, digits = 4)

# set the scale of x
E.mode <- c(results$data$B_inner %*% results$point$alpha_w)
# quantile(E.mode, c(0.0025, 0.9975))
# x.quantile <- quantile(dat.fit$x, c(0.0025, 0.9975))
x.ll <- max(min(E.mode), min(dat.fit$x))
x.ul <- min(max(E.mode), max(dat.fit$x))

E.eval <- seq(x.ll, x.ul, length.out = 500)
summ <- summary(results, plot = FALSE, E0 = E0, E.eval = E.eval)


pw1 <- ggplot(data = summ$est$wl, aes(x = l, y = mode))
pw1 <- pw1 + geom_line(linewidth = 1.5, color = "#2c7bb6")
pw1 <- pw1 + geom_ribbon(aes(x = l, ymax = ul, ymin = ll), alpha = 0.3, fill = "#2c7bb6")
pw1 <- pw1 + xlab("Lag") + ylab("Weight")
pw1 <- pw1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
    legend.text=element_text(size=GGPLOTLEGENDSIZE),
    plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
    plot.subtitle=element_text(size=GGPLOTSTSIZE),
    legend.spacing.y = unit(0.3, 'cm') ) + 
    ylim(c(-0.7,1.3))
pw1 <- pw1 + geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7) + 
    scale_x_continuous(breaks=seq(0,maxL,1))
pw1 <- pw1 + labs(# title = paste0(CDname, ": Respiratory"),
                  subtitle = "ACE-DLNM: w")
pw_dlnm <- pw1








# pf1 <- ggplot(data = summ$est$fE0, aes(x = E, y = mode))
draw.fE <- summ$est$fE0
draw.fE$ul <- sapply(draw.fE$ul, function(a) min(a, ylimu))
draw.fE$ll <- sapply(draw.fE$ll, function(a) max(a, yliml))
pf1 <- ggplot(data = draw.fE, aes(x = E, y = mode))
pf1 <- pf1 + geom_line(linewidth = 1.5, color = "#2c7bb6")
pf1 <- pf1 + geom_ribbon(aes(ymax = ul, ymin = ll), alpha = 0.3, fill = "#2c7bb6")
pf1 <- pf1 + ylab("ACERF") + xlab("Adaptive Cumulative Exposure")
pf1 <- pf1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
    legend.text=element_text(size=GGPLOTLEGENDSIZE),
    plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
    plot.subtitle=element_text(size=GGPLOTSTSIZE),
    legend.spacing.y = unit(0.3, 'cm') ) + 
    geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7) +
    ylim(c(yliml, ylimu)) + 
    xlim(c(x.ll, x.ul))
pf1 <- pf1 + labs(# title = "  ",
           subtitle = "ACE-DLNM: f")
pf_dlnm <- pf1





pg1 <- ggplot(data = filter(summ$est$smooth, var == "temp"), aes(x = x, y = mode))
pg1 <- pg1 + geom_line(linewidth = 1.5, color = "#2c7bb6")
pg1 <- pg1 + geom_ribbon(aes(x = x, ymax = ul, ymin = ll), alpha = 0.3, fill = "#2c7bb6")
pg1 <- pg1 + ylab("Association of Temperature") + xlab("Temperature")
pg1 <- pg1 + theme(text = element_text(size = GGPLOTTEXTSIZE),
    legend.text=element_text(size=GGPLOTLEGENDSIZE),
    plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
    plot.subtitle=element_text(size=GGPLOTSTSIZE),
    legend.spacing.y = unit(0.3, 'cm') ) + 
    ylim(c(-0.25, 0.13))
pg1 <- pg1 + labs(subtitle = "ACE-DLNM: h(temp)")
ptemp_dlnm <- pg1


starting.date <- as.Date("2001-01-01")

gt_df <- filter(summ$est$smooth, var == "t")
gt_df$t <- starting.date + gt_df$x
pg2 <- ggplot(data = gt_df, aes(x = t, y = mode))
pg2 <- pg2 + geom_line(linewidth = 1.5, color = "#2c7bb6")
pg2 <- pg2 + geom_ribbon(aes(x = t, ymax = ul, ymin = ll), alpha = 0.3, fill = "#2c7bb6")
pg2 <- pg2 + ylab("Long-Term Trend") + xlab("Time")
pg2 <- pg2 + theme(text = element_text(size = GGPLOTTEXTSIZE),
    legend.text=element_text(size=GGPLOTLEGENDSIZE),
    plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
    plot.subtitle=element_text(size=GGPLOTSTSIZE),
    legend.spacing.y = unit(0.3, 'cm') ) + 
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") + 
    ylim(c(-0.22, 0.3))
pg2 <- pg2 + labs(subtitle = "ACE-DLNM: h(t)")
ptime_dlnm <- pg2



gt_df <- filter(summ$est$smooth, var == "month.num")
pg2 <- ggplot(data = gt_df, aes(x = x, y = mode))
pg2 <- pg2 + geom_line(linewidth = 1.5, color = "#2c7bb6")
pg2 <- pg2 + geom_ribbon(aes(x = x, ymax = ul, ymin = ll), alpha = 0.3, fill = "#2c7bb6")
pg2 <- pg2 + ylab("Seasonality") + xlab("Month")
pg2 <- pg2 + theme(text = element_text(size = GGPLOTTEXTSIZE),
    legend.text=element_text(size=GGPLOTLEGENDSIZE),
    plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
    plot.subtitle=element_text(size=GGPLOTSTSIZE),
    legend.spacing.y = unit(0.3, 'cm') )  + 
    scale_x_continuous(breaks=seq(0,12,1)) + 
    ylim(c(-0.48, 0.4))
pg2 <- pg2 + labs(subtitle = "ACE-DLNM: h(month)")
pmonth_dlnm <- pg2







### mgcv ############
colors <- c("#ca0020", "#f4a582", "#bababa", "#404040")
labels <- c("lag0", "avg lag0-1", "avg lag0-7", "avg lag0-14")


if(maxL == 14){
dat.fit.mgcv <- dat.fit

dat.fit.mgcv$avglag01 <- (dat.fit.mgcv$x+ dplyr::lag(dat.fit.mgcv$x, 1))/sqrt(2)

dat.fit.mgcv$avglag02 <- (dat.fit.mgcv$x + 
                          dplyr::lag(dat.fit.mgcv$x, 1) + 
                          dplyr::lag(dat.fit.mgcv$x, 2) + 
                          dplyr::lag(dat.fit.mgcv$x, 3) + 
                          dplyr::lag(dat.fit.mgcv$x, 4) + 
                          dplyr::lag(dat.fit.mgcv$x, 5) + 
                          dplyr::lag(dat.fit.mgcv$x, 6) + 
                          dplyr::lag(dat.fit.mgcv$x, 7)
                          )/sqrt(8)



dat.fit.mgcv$avglag03 <- (dat.fit.mgcv$x + 
                          dplyr::lag(dat.fit.mgcv$x, 1) + 
                          dplyr::lag(dat.fit.mgcv$x, 2) + 
                          dplyr::lag(dat.fit.mgcv$x, 3) + 
                          dplyr::lag(dat.fit.mgcv$x, 4) + 
                          dplyr::lag(dat.fit.mgcv$x, 5) + 
                          dplyr::lag(dat.fit.mgcv$x, 6) + 
                          dplyr::lag(dat.fit.mgcv$x, 7) + 
                          dplyr::lag(dat.fit.mgcv$x, 8) + 
                          dplyr::lag(dat.fit.mgcv$x, 9) + 
                          dplyr::lag(dat.fit.mgcv$x, 10) + 
                          dplyr::lag(dat.fit.mgcv$x, 11) + 
                          dplyr::lag(dat.fit.mgcv$x, 12) + 
                          dplyr::lag(dat.fit.mgcv$x, 13) + 
                          dplyr::lag(dat.fit.mgcv$x, 14)
                          )/sqrt(15)

mod.fit.mgcv <- bam(y~s(x, bs = "bs", k = 20) + 
                      s(month.num, bs = "cc", k = 10) + 
                      s(t, bs = "bs", k = 10) + 
                      s(temp, bs = "bs", k = 10) + 
                      DOW,
                    data = dat.fit.mgcv, 
                    family = nb(),
                    discrete = TRUE)

mod.fit.mgcv.pois <- bam(y~s(x, bs = "bs", k = 20) + 
                      s(month.num, bs = "cc", k = 10) + 
                      s(t, bs = "bs", k = 10) + 
                      s(temp, bs = "bs", k = 10) + 
                      DOW,
                    data = dat.fit.mgcv, 
                    family = "poisson",
                    discrete = TRUE)

mod.fit.avg01.mgcv <- bam(y~s(avglag01, bs = "bs", k = 20) + 
                            s(month.num, bs = "cc", k = 10) + 
                            s(t, bs = "bs", k = 10) + 
                            s(temp, bs = "bs", k = 10) + 
                            DOW,
                         data = dat.fit.mgcv, 
                         family = nb(),
                         discrete = TRUE)

mod.fit.avg02.mgcv <- bam(y~s(avglag02, bs = "bs", k = 20) + 
                            s(month.num, bs = "cc", k = 10) + 
                            s(t, bs = "bs", k = 10) + 
                            s(temp, bs = "bs", k = 10) + 
                            DOW,
                         data = dat.fit.mgcv, 
                         family = nb(),
                         discrete = TRUE)

mod.fit.avg03.mgcv <- bam(y~s(avglag03, bs = "bs", k = 20) + 
                            s(month.num, bs = "cc", k = 10) + 
                            s(t, bs = "bs", k = 10) + 
                            s(temp, bs = "bs", k = 10) + 
                            DOW,
                         data = dat.fit.mgcv, 
                         family = nb(),
                         discrete = TRUE)

                     
summ.mgcv <- summary(mod.fit.mgcv)
summ.mgcv.table <- as.data.frame(summ.mgcv$p.table)
summ.mgcv.table$ul <- summ.mgcv.table[,1] + 1.96*summ.mgcv.table[,2]
summ.mgcv.table$ll <- summ.mgcv.table[,1] - 1.96*summ.mgcv.table[,2]
knitr::kable(select(summ.mgcv.table, Estimate, ll, ul), digits = 3)
print(family(mod.fit.mgcv)$getTheta(), digits = 4)

summ.mgcv <- summary(mod.fit.avg03.mgcv)
summ.mgcv.table <- as.data.frame(summ.mgcv$p.table)
summ.mgcv.table$ul <- summ.mgcv.table[,1] + 1.96*summ.mgcv.table[,2]
summ.mgcv.table$ll <- summ.mgcv.table[,1] - 1.96*summ.mgcv.table[,2]
knitr::kable(select(summ.mgcv.table, Estimate, ll, ul), digits = 3)
print(family(mod.fit.avg03.mgcv)$getTheta(), digits = 4)

## RQR
mod <- mod.fit.mgcv
mu <- mod$fitted.values # mu
theta <- exp(family(mod)$getTheta()) # theta
y <- dat.fit.mgcv$y

residual <- mapply(function(yi, mui){
        yilim <- max(yi-1,0)
        ai <- pnbinom(yilim, size = theta, mu = mui)
        bi <- pnbinom(yi, size = theta, mu = mui)
        ui <- runif(1, min = ai, max = bi)
        ri <- qnorm(ui)
        return(ri)
    },
    y, mu)

p.qq <- ggplot(data.frame(xx = residual), aes(sample = xx)) + stat_qq() + stat_qq_line() + xlab("Theoretical quantiles") + ylab("Sample quantiles") + 
    theme(text = element_text(size = GGPLOTTEXTSIZE),
    legend.text=element_text(size=GGPLOTLEGENDSIZE),
    plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
    plot.subtitle=element_text(size=GGPLOTSTSIZE),
    legend.spacing.y = unit(0.3, 'cm') ) + 
    labs(subtitle = "GAM lag0: QQ-Plot")
qq_mgcv <- p.qq
# ggsave(p.qq, file = file.path(paste0(CD,"-CP-morb-mgcv-qqplot.pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)

q.res <- ggplot(data.frame(t = dat.fit.mgcv$t, xx = residual), aes(x = t, y = xx)) +
        geom_point() + geom_hline(yintercept = 0) + 
        ylim(-5,5) + 
        scale_x_continuous(breaks = seq(min(dat.fit.mgcv$t), max(dat.fit.mgcv$t), by = round((max(dat.fit.mgcv$t) - min(dat.fit.mgcv$t))/30))) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        xlab("Time") + ylab("Randomized Quantile Residuals") + 
         theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') ) + 
        labs(subtitle = "GAM lag0: Residual Plot")
res_mgcv <- q.res
# ggsave(q.res, file = file.path(paste0(CD, "-CP-morb-mgcv-resplot.pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)


res.mgcv <- function(mod, title){
    mu <- mod$fitted.values # mu
    theta <- exp(family(mod)$getTheta()) # theta
    y <- dat.fit.mgcv$y

    residual <- mapply(function(yi, mui){
            yilim <- max(yi-1,0)
            ai <- pnbinom(yilim, size = theta, mu = mui)
            bi <- pnbinom(yi, size = theta, mu = mui)
            ui <- runif(1, min = ai, max = bi)
            ri <- qnorm(ui)
            return(ri)
        },
        y, mu)

    p.qq <- ggplot(data.frame(xx = residual), aes(sample = xx)) + stat_qq() + stat_qq_line() + xlab("Theoretical quantiles") + ylab("Sample quantiles") + 
        theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') ) + 
        labs(subtitle = "GAM: QQ-Plot")
    ggsave(p.qq, file = file.path("figures", paste0(CD,"-", title, "-qqplot.pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)

    q.res <- ggplot(data.frame(t = dat.fit.mgcv$t, xx = residual), aes(x = t, y = xx)) +
            geom_point() + geom_hline(yintercept = 0) + 
            ylim(-5,5) + 
            scale_x_continuous(breaks = seq(min(dat.fit.mgcv$t), max(dat.fit.mgcv$t), by = round((max(dat.fit.mgcv$t) - min(dat.fit.mgcv$t))/30))) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
            xlab("Time") + ylab("Randomized Quantile Residuals") + 
            theme(text = element_text(size = GGPLOTTEXTSIZE),
            legend.text=element_text(size=GGPLOTLEGENDSIZE),
            plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
            plot.subtitle=element_text(size=GGPLOTSTSIZE),
            legend.spacing.y = unit(0.3, 'cm') ) + 
            labs(subtitle = "GAM: Residual Plot")
    ggsave(q.res, file = file.path("figures", paste0(CD, "-", title, "-resplot.pdf")), width = PLOTWIDTH, height = PLOTHEIGHT)


}

res.mgcv(mod.fit.avg01.mgcv, "avg01")
res.mgcv(mod.fit.avg02.mgcv, "avg07")
res.mgcv(mod.fit.avg03.mgcv, "avg014")

preddat <- data.frame(x = seq(min(dat.fit.mgcv$x, na.rm = T), max(dat.fit.mgcv$x, na.rm = T), length.out=300),
                      avglag01 = seq(min(dat.fit.mgcv$x, na.rm = T), max(dat.fit.mgcv$x, na.rm = T), length.out=300),
                      avglag02 = seq(min(dat.fit.mgcv$x, na.rm = T), max(dat.fit.mgcv$x, na.rm = T), length.out=300),
                      avglag03 = seq(min(dat.fit.mgcv$x, na.rm = T), max(dat.fit.mgcv$x, na.rm = T), length.out=300),
                      temp = 20,
                      month.num = 1,
                      DOW = dat.fit.mgcv$DOW[1],
                      t = 100)
preddat.E0 <- data.frame(x = rep(x0,300),
                      avglag01 = x0,
                      avglag02 = x0,
                      avglag03 = x0,
                      temp = 20,
                      month.num = 1,
                      DOW = dat.fit.mgcv$DOW[1],
                      t = 100)  




## f(x) - f(x0)

get_con <- function(mod, dat0, dat1) {
     pred0 <- predict(mod, newdata = dat0, se.fit = TRUE, type = 'link')
     pred0.plus1 <- predict(mod, newdata = dat1, se.fit = TRUE, type = 'link')
     log_rr <- pred0.plus1$fit - pred0$fit
     model_matrix0 <- model.matrix(mod, newdata = dat0)
     model_matrix0.plus1 <- model.matrix(mod, newdata = dat1)

     beta_var <- vcov(mod)
     # U <- chol(beta_var)
     # log_rr_se <- as.numeric(sqrt(colSums( (U %*% t(model_matrix0.plus1-model_matrix0))^2 )))

     log_rr_se <- sqrt(diag((model_matrix0.plus1-model_matrix0) %*% beta_var %*% t((model_matrix0.plus1-model_matrix0))))
     log_rr_l = log_rr - 1.96*log_rr_se
     log_rr_u = log_rr + 1.96*log_rr_se
     return(data.frame(x = dat0$x, log_rr = log_rr, log_rr_l = log_rr_l, log_rr_u = log_rr_u))
}

E0_mgcv0 <- get_con(mod.fit.mgcv, preddat.E0, preddat)
E0_mgcv01 <- get_con(mod.fit.avg01.mgcv, preddat.E0, preddat)
E0_mgcv02 <- get_con(mod.fit.avg02.mgcv, preddat.E0, preddat)
E0_mgcv03 <- get_con(mod.fit.avg03.mgcv, preddat.E0, preddat)

E0_mgcv0$exposure <- "lag0"
E0_mgcv0$x <- preddat$x
E0_mgcv01$exposure <- "avg lag0-1"
E0_mgcv01$x <- preddat$x
E0_mgcv02$exposure <- "avg lag0-2"
E0_mgcv02$x <- preddat$x
E0_mgcv03$exposure <- "avg lag0-3"
E0_mgcv03$x <- preddat$x

E0_mgcv <- rbind(E0_mgcv0,
                 E0_mgcv01,
                 E0_mgcv02,
                 E0_mgcv03)
E0_mgcv$exposure <- factor(E0_mgcv$exposure, levels = c("lag0", "avg lag0-1", "avg lag0-2", "avg lag0-3"))
E0_mgcv$fitted <- E0_mgcv$log_rr
E0_mgcv$ul <- E0_mgcv$log_rr_u
E0_mgcv$ll <- E0_mgcv$log_rr_l

E0_mgcv$ul_draw <- sapply(E0_mgcv$ul, function(a) min(a, ylimu))
E0_mgcv$ll_draw <- sapply(E0_mgcv$ll, function(a) max(a, yliml))

p <- ggplot(data = E0_mgcv, aes(x = x, y = fitted)) + 
 geom_line(aes(color = exposure),linewidth = 1.5) + 
 geom_ribbon(aes(x = x, ymax = ul_draw, ymin = ll_draw, fill = exposure), alpha = 0.3) + #fill = "#2c7bb6"
 ylab("ERF") + xlab("Cumulative Exposure") + 
 theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=0.8*GGPLOTTEXTSIZE),
        legend.title = element_text(size = 0.8*GGPLOTTEXTSIZE),
        plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') ) + 
 labs(#title = "  ",
      subtitle = "GAM: f") + 
 ylim(c(yliml,ylimu)) + 
 xlim(c(x.ll, x.ul)) + 
 geom_hline(yintercept = 0, linewidth = 1, alpha = 0.7) + 
 scale_color_manual(values = colors, name = "Exposure", labels = labels) + 
 scale_fill_manual(values = colors, name = "Exposure", labels = labels)
pf_mgcv <- p



## plot h(temp)
preddat <- data.frame(x = 20,
                      avglag01 = 20,
                      avglag02 = 20,
                      avglag03 = 20,
                      temp = seq(min(dat.fit.mgcv$temp, na.rm = T), max(dat.fit.mgcv$temp, na.rm = T), length.out=300),
                      DOW = dat.fit.mgcv$DOW[1],
                      month.num = 1,
                      t = 100)
pred0 <- predict(mod.fit.mgcv, newdata = preddat, se.fit = TRUE, type = 'terms')
fitted0 <- pred0$fit[,5]
pred01 <- predict(mod.fit.avg01.mgcv, newdata = preddat, se.fit = TRUE, type = 'terms')
fitted01 <- pred01$fit[,5]
pred02 <- predict(mod.fit.avg02.mgcv, newdata = preddat, se.fit = TRUE, type = 'terms')
fitted02 <- pred02$fit[,5]
pred03 <- predict(mod.fit.avg03.mgcv, newdata = preddat, se.fit = TRUE, type = 'terms')
fitted03 <- pred03$fit[,5]

se0 <- pred0$se.fit[,5]
ul0 <- fitted0 + 1.96*se0
ll0 <- fitted0 - 1.96*se0
se01 <- pred01$se.fit[,5]
ul01 <- fitted01 + 1.96*se01
ll01 <- fitted01 - 1.96*se01
se02 <- pred0$se.fit[,5]
ul02 <- fitted02 + 1.96*se02
ll02 <- fitted02 - 1.96*se02
se03 <- pred0$se.fit[,5]
ul03 <- fitted03 + 1.96*se03
ll03 <- fitted03 - 1.96*se03

gam.df0 <- data.frame(x = preddat$temp,
                     fitted = fitted0,
                     ul = ul0,
                     ll = ll0,
                     exposure = "lag0")
gam.df01 <- data.frame(x = preddat$temp,
                      fitted = fitted01,
                      ul = ul01,
                      ll = ll01,
                      exposure = "avg lag0-1")
gam.df02 <- data.frame(x = preddat$temp,
                      fitted = fitted02,
                      ul = ul02,
                      ll = ll02,
                      exposure = "avg lag0-2")
gam.df03 <- data.frame(x = preddat$temp,
                      fitted = fitted03,
                      ul = ul03,
                      ll = ll03,
                      exposure = "avg lag0-3")
gam.df <- rbind(gam.df0,
                gam.df01,
                gam.df02,
                gam.df03)
gam.df$exposure <- factor(gam.df$exposure, levels = c("lag0", "avg lag0-1", "avg lag0-2", "avg lag0-3"))

gam.df$ul_draw <- sapply(gam.df$ul, function(a) min(a, 0.13))
gam.df$ll_draw <- sapply(gam.df$ll, function(a) max(a, -0.25))

p <- ggplot(data = gam.df, aes(x = x, y = fitted)) + 
 geom_line(linewidth = 1.5, aes(color = exposure)) + 
 geom_ribbon(aes(x = x, ymax = ul_draw, ymin = ll_draw,fill = exposure), alpha = 0.3) + 
 ylab("Association of Temperature") + xlab("Temperature") + 
 theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') ) + 
 ylim(c(-0.25, 0.13)) + 
 labs(subtitle = "GAM: h(temp)") + 
 scale_color_manual(values = colors, name = "Exposure", labels = labels) + 
 scale_fill_manual(values = colors, name = "Exposure", labels = labels)
ptemp_mgcv <- p


## plot h(t)
preddat <- data.frame(x = 20,
                      avglag01 = 20,
                      avglag02 = 20,
                      avglag03 = 20,
                      temp = 20,
                      DOW = dat.fit.mgcv$DOW[1],
                      month.num = 1,
                      t = dat.fit.mgcv$t)
pred0 <- predict(mod.fit.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted0 <- pred0$fit[,4]
pred01 <- predict(mod.fit.avg01.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted01 <- pred01$fit[,4]
pred02 <- predict(mod.fit.avg02.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted02 <- pred02$fit[,4]
pred03 <- predict(mod.fit.avg03.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted03 <- pred03$fit[,4]

se0 <- pred0$se.fit[,4]
ul0 <- fitted0 + 1.96*se0
ll0 <- fitted0 - 1.96*se0
se01 <- pred01$se.fit[,4]
ul01 <- fitted01 + 1.96*se01
ll01 <- fitted01 - 1.96*se01
se02 <- pred02$se.fit[,4]
ul02 <- fitted02 + 1.96*se02
ll02 <- fitted02 - 1.96*se02
se03 <- pred03$se.fit[,4]
ul03 <- fitted03 + 1.96*se03
ll03 <- fitted03 - 1.96*se03



gam.df0 <- data.frame(t = dat.fit.mgcv$t,
                      fitted = fitted0,
                      ul = ul0,
                      ll = ll0,
                      exposure = "lag0")
gam.df01 <- data.frame(t = dat.fit.mgcv$t,
                       fitted = fitted01,
                       ul = ul01,
                       ll = ll01,
                       exposure = "avg lag0-1")
gam.df02 <- data.frame(t = dat.fit.mgcv$t,
                       fitted = fitted02,
                       ul = ul02,
                       ll = ll02,
                       exposure = "avg lag0-2")
gam.df03 <- data.frame(t = dat.fit.mgcv$t,
                       fitted = fitted03,
                       ul = ul03,
                       ll = ll03,
                       exposure = "avg lag0-3")

gam.df <- rbind(gam.df0,
                gam.df01,
                gam.df02,
                gam.df03)  
gam.df$exposure <- factor(gam.df$exposure, levels = c("lag0", "avg lag0-1", "avg lag0-2", "avg lag0-3"))
gam.df$t <- starting.date + gam.df$t
gam.df$ul_draw <- sapply(gam.df$ul, function(a) min(a, 0.3))
gam.df$ll_draw <- sapply(gam.df$ll, function(a) max(a, -0.22))


pg2 <- ggplot(data = gam.df, aes(x = t, y = fitted))
pg2 <- pg2 + geom_line(linewidth = 1.5, aes(color = exposure))
pg2 <- pg2 + geom_ribbon(aes(x = t, ymax = ul_draw, ymin = ll_draw, fill = exposure), alpha = 0.3)
pg2 <- pg2 + ylab("Long-Term Trend") + xlab("Time")
pg2 <- pg2 + theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') ) + 
    labs(subtitle = "GAM: h(t)") + 
    ylim(c(-0.22, 0.3)) + 
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") + 
    scale_color_manual(values = colors, name = "Exposure", labels = labels) + 
    scale_fill_manual(values = colors, name = "Exposure", labels = labels)

ptime_mgcv <- pg2



## plot h(month)
preddat <- data.frame(x = 20,
                      avglag01 = 20,
                      avglag02 = 20,
                      avglag03 = 20,
                      temp = 20,
                      DOW = dat.fit.mgcv$DOW[1],
                      month.num = seq(1,12,length.out = 100),
                      t = 1)
pred0 <- predict(mod.fit.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted0 <- pred0$fit[,3]
pred01 <- predict(mod.fit.avg01.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted01 <- pred01$fit[,3]
pred02 <- predict(mod.fit.avg02.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted02 <- pred02$fit[,3]
pred03 <- predict(mod.fit.avg03.mgcv, newdata = preddat, se.fit = TRUE, unconditional = FALSE, type = 'terms')
fitted03 <- pred03$fit[,3]

se0 <- pred0$se.fit[,3]
ul0 <- fitted0 + 1.96*se0
ll0 <- fitted0 - 1.96*se0
se01 <- pred01$se.fit[,3]
ul01 <- fitted01 + 1.96*se01
ll01 <- fitted01 - 1.96*se01
se02 <- pred02$se.fit[,3]
ul02 <- fitted02 + 1.96*se02
ll02 <- fitted02 - 1.96*se02
se03 <- pred03$se.fit[,3]
ul03 <- fitted03 + 1.96*se03
ll03 <- fitted03 - 1.96*se03

gam.df0 <- data.frame(month = preddat$month.num,
                     fitted = fitted0,
                     ul = ul0,
                     ll = ll0,
                     exposure = "lag0")
gam.df01 <- data.frame(month = preddat$month.num,
                     fitted = fitted01,
                     ul = ul01,
                     ll = ll01,
                     exposure = "avg lag0-1")
gam.df02 <- data.frame(month = preddat$month.num,
                     fitted = fitted02,
                     ul = ul02,
                     ll = ll02,
                     exposure = "avg lag0-2")
gam.df03 <- data.frame(month = preddat$month.num,
                     fitted = fitted03,
                     ul = ul03,
                     ll = ll03,
                     exposure = "avg lag0-3")

gam.df <-  rbind(gam.df0,
                 gam.df01,
                 gam.df02,
                 gam.df03)
gam.df$exposure <- factor(gam.df$exposure, levels = c("lag0", "avg lag0-1", "avg lag0-2", "avg lag0-3"))


gam.df$ul_draw <- sapply(gam.df$ul, function(a) min(a, 0.4))
gam.df$ll_draw <- sapply(gam.df$ll, function(a) max(a, -0.48))

pg2 <- ggplot(data = gam.df, aes(x = month, y = fitted))
pg2 <- pg2 + geom_line(linewidth = 1.5, aes(color = exposure))
pg2 <- pg2 + geom_ribbon(aes(x = month, ymax = ul_draw, ymin = ll_draw, fill = exposure), alpha = 0.3)
pg2 <- pg2 + ylab("Seasonality") + xlab("month")
pg2 <- pg2 + theme(text = element_text(size = GGPLOTTEXTSIZE),
        legend.text=element_text(size=GGPLOTLEGENDSIZE),
        plot.title = element_text(size = 1.2*GGPLOTSTSIZE, face = "bold"),
        plot.subtitle=element_text(size=GGPLOTSTSIZE),
        legend.spacing.y = unit(0.3, 'cm') ) + 
        labs(subtitle = "GAM: h(month)") + 
        ylim(c(-0.48, 0.4)) + 
        scale_x_continuous(breaks=seq(0,12,1)) + 
        scale_color_manual(values = colors, name = "Exposure", labels = labels) + 
        scale_fill_manual(values = colors, name = "Exposure", labels = labels)

pmonth_mgcv <- pg2


## save plots

library(ggpubr)
p1 <- ggarrange(pw_dlnm, pf_dlnm, pf_mgcv, ncol=3, nrow=1, common.legend = TRUE, legend="right")
p1 <- annotate_figure(p1, top = text_grob(paste0(CDname, ": Respiratory"), 
               color = "black", face = "bold", size = 1.2*GGPLOTSTSIZE))
ggsave(p1, file = file.path("figures", paste0(CD, "-Pulm-w-f.pdf")), width = 3*PLOTWIDTH, height = PLOTHEIGHT)

ptemp <- ggarrange(ptemp_dlnm, ptemp_mgcv, nrow=1, common.legend = TRUE, legend="right")
ggsave(ptemp, file = file.path("figures", paste0(CD, "-Pulm-temp.pdf")), width = 2*PLOTWIDTH, height = PLOTHEIGHT)

ptime <- ggarrange(ptime_dlnm, ptime_mgcv, nrow=1, common.legend = TRUE, legend="right")
ggsave(ptime, file = file.path("figures", paste0(CD, "-Pulm-time.pdf")), width = 2*PLOTWIDTH, height = PLOTHEIGHT)

pmonth <- ggarrange(pmonth_dlnm, pmonth_mgcv, nrow=1, common.legend = TRUE, legend="right")
ggsave(pmonth, file = file.path("figures", paste0(CD, "-Pulm-month.pdf")), width = 2*PLOTWIDTH, height = PLOTHEIGHT)

p_qq <- ggarrange(qq_dlnm, qq_mgcv, nrow=1, common.legend = TRUE, legend="right")
ggsave(p_qq, file = file.path("figures", paste0(CD, "-Pulm-qq.pdf")), width = 2*PLOTWIDTH, height = PLOTHEIGHT)

p_res <- ggarrange(res_dlnm, res_mgcv, nrow=1, common.legend = TRUE, legend="right")
ggsave(p_res, file = file.path("figures", paste0(CD, "-Pulm-res.pdf")), width = 2*PLOTWIDTH, height = PLOTHEIGHT)

}
}
}
