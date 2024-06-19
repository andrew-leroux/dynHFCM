## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  # collapse = TRUE,
  class.source = "fold-show",
  echo = TRUE#,
  #comment = "#>"
)
library(tidyverse)
library(ggplot2)
library(mgcv)

## ----load_data, include = TRUE------------------------------------------------
library(dynHFCM)
data("hfcm_sim_data")
str(hfcm_sim_data)

## ----create_landmark_data-----------------------------------------------------
## grid of observation
tmin <- 0 # time of first observation
tmax <- 1 # time of last observation (administrative censoring)

## landmark times and prediciton windows
S <- seq(tmin, tmax, len=51)[-51]       # landmark times for both landmark models
W1 <- diff(c(S, tmax))                  # prediction windows for first landmark model (w=0.04)
W2 <- rep(Inf, length(S))               # prediction windows for second landmark model (w=infinity)


## create landmark datasets for estiamting the landmark historical functional Cox models
data_lm_w1 <- make_lm_data(data=hfcm_sim_data, vars_tv="Z",vars_tf=c(), id="id", event_time="T", 
                           time="tind", S=S, wide=FALSE, censor_var="E", w=W1)
data_lm_w2 <- make_lm_data(data=hfcm_sim_data, vars_tv="Z",vars_tf=c(), id="id", event_time="T", 
                           time="tind", S=S, wide=FALSE, censor_var="E", w=W2)

## ----view_landmark_data_w1----------------------------------------------------
str(data_lm_w1)
sum(data_lm_w1$id == 1)

## ----view_landmark_data_w2----------------------------------------------------
str(data_lm_w2)

## ----fit_landmark_models------------------------------------------------------
# fit landmark model with w=0.02
st_time_lm_w1  <- Sys.time()
fit_lm_w1 <- gam(cbind(event_time_lm, s) ~ s(tind, smat, by=Z_int_cn), family=cox.ph, weights=event_lm, data=data_lm_w1)
(fit_time_lm_w1 <- difftime(Sys.time(), st_time_lm_w1, units="mins"))
        
# fit landmark model with w=infinity
st_time_lm_w2  <- Sys.time()
fit_lm_w2 <- gam(cbind(event_time_lm, s) ~ s(tind, smat, by=Z_int_cn), family=cox.ph, weights=event_lm, data=data_lm_w2)
(fit_time_lm_w2 <- difftime(Sys.time(), st_time_lm_w2, units="mins"))

## ----get_estimated_coefficients-----------------------------------------------
nupred_coef <- 100
nspred_coef <- 100
upred_coef <- seq(tmin, max(S), len=nupred_coef)
spred_coef <- seq(min(S), max(S), len=nspred_coef)
## get estimated coefficient
# to get the estimated coefficients, need to evaluate the coefficient using predict.gam with type="terms"
df_pred_lm_coef <- data.frame(tind = rep(upred_coef, nspred_coef),
                              smat = rep(spred_coef, each=nupred_coef),
                              Z_int_cn = 1)
df_pred_lm_coef <- subset(df_pred_lm_coef, tind < smat)
gamma_hat_lm_w1 <- predict(fit_lm_w1, newdata=df_pred_lm_coef, type="terms")
gamma_hat_lm_w2 <- predict(fit_lm_w2, newdata=df_pred_lm_coef, type="terms")
## combine results into a single data frame for plotting
df_pred_lm_coef$gamma_hat_w1 <- gamma_hat_lm_w1[,"s(tind,smat):Z_int_cn"]
df_pred_lm_coef$gamma_hat_w2 <- gamma_hat_lm_w2[,"s(tind,smat):Z_int_cn"]

## ----plot_estimated_coefficients, fig.height=5, fig.width=15------------------
## add in the true coefficient value
df_pred_lm_coef <- 
  df_pred_lm_coef %>% 
  mutate(gamma_true = cos(2 * pi * (smat-tind) / 0.5)  * 10) 
## plot the results
# transform to long format for faceting then plot
df_pred_lm_coef %>% 
  pivot_longer(cols=c("gamma_hat_w1","gamma_hat_w2","gamma_true"), 
               names_to = "result", values_to="gamma_hat") %>% 
  mutate(result = factor(result, levels=c("gamma_true","gamma_hat_w1","gamma_hat_w2"),
                         labels=c("Truth","Estimate w=0.04","Estimate w=Infinity"))) %>% 
  ggplot() + 
  geom_raster(aes(x=tind,y=smat,fill=gamma_hat)) + facet_grid(~result) + 
  scale_fill_gradientn(colours=fields::tim.colors(50)) + theme_classic() + xlab("Historical Time (u)") + 
  ylab("Landmark Time (s)")

