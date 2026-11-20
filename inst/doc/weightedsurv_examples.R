## ----suggested-packages, include = FALSE--------------------------------------
# Packages in Suggests: a chunk that uses one is evaluated only when it is installed
has <- vapply(c("DiagrammeR", "dplyr", "data.table", "gt", "adjustedCurves", "pammtools"),
              requireNamespace, logical(1), quietly = TRUE)
# numbers in the text that need such a package
if_has <- function(pkg, x) if (all(has[pkg])) x else paste0("(needs ", paste(pkg, collapse = " and "), ")")

## ----echo = FALSE, eval = has[["DiagrammeR"]]---------------------------------
library(DiagrammeR)
# Define the flowchart in DOT syntax with horizontal layout
flowchart <- grViz('digraph wt_rg_S_flowchart {
  rankdir = LR
  node [shape = box, style = filled, fillcolor = "#e3e3e3"]
  A [label = "Start: Call wt.rg.S"]
  B [label = "Input Validation"]
  C [label = "Which scheme?"]
  D [label = "FH: S_left^rho * (1 - S_left)^gamma"]
  E [label = "Schemper: S_left / Scensor_left"]
  F [label = "XO: S_left / Ybar"]
  G [label = "MB: 1 / max(S_left, Shat_tzero)"]
  H [label = "custom_time: w0.tau or w1.tau by t.tau"]
  I [label = "fh_exp1: exp(S_left^0.5 * (1 - S_left)^0.5)"]
  J [label = "fh_exp2: pmin(exp(wt01), wmax)"]
  K [label = "Return weights"]
  L [label = "If details=TRUE, return list; else, return weights"]
  M [label = "End"]

  A -> B -> C
  C -> D [label = "fh"]
  C -> E [label = "schemper"]
  C -> F [label = "XO"]
  C -> G [label = "MB"]
  C -> H [label = "custom_time"]
  C -> I [label = "fh_exp1"]
  C -> J [label = "fh_exp2"]
  D -> K
  E -> K
  F -> K
  G -> K
  H -> K
  I -> K
  J -> K
  K -> L -> M
}
', width = 800, height =400)
# Display the flowchart
flowchart

## ----message = FALSE, warning = FALSE-----------------------------------------
oldpar <- par(no.readonly = TRUE)
library(survival)
# To install weightedsurv: the released version from CRAN
# install.packages("weightedsurv")
# or the development version from GitHub
# devtools::install_github("larry-leon/weightedsurv")
library(weightedsurv)
if (has[["dplyr"]]) library(dplyr)
# formatting of the numbers quoted in the text
fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

## -----------------------------------------------------------------------------

library(survival)
df_gbsg <- gbsg
df_gbsg$tte <- df_gbsg$rfstime / 30.4375
df_gbsg$event <- df_gbsg$status
df_gbsg$treat <- df_gbsg$hormon
df_gbsg$grade3 <- ifelse(df_gbsg$grade == "3", 1, 0)


# Note that these are default names
# For alternative names, for example if the 
# time-to-event outcome (tte) is time_months
# then tte.name <- "time_months"
tte.name <- "tte"
event.name <- "event"
treat.name <- "treat"
arms <- c("treat", "control")


## -----------------------------------------------------------------------------
# Returns standard log-rank (FH(0,0) via scheme="fh" rho,gamma = 0 in scheme_params)

dfcount_gbsg <- df_counting(df=df_gbsg, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = arms, 
by.risk = 12, scheme = "fh", scheme_params = list(rho = 0, gamma =0), lr.digits = 4, cox.digits =3)

# test default names
#test <- df_counting(df=df_gbsg, by.risk = 12, scheme = "fh", scheme_params = list(rho = 0, gamma =0), lr.digits = 4, cox.digits =3)


## ----eval=FALSE---------------------------------------------------------------
# # MB weighting
# #test <- df_counting(df=df_gbsg, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, scheme = "MB", scheme_params = list(mb_tstar = 12))
# # 0/1 weighting (0 for time <= t.tau, 1 for time > t.tau)
# #test <- df_counting(df=df_gbsg, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, scheme = "custom_time", scheme_params = list(t.tau = 6, w0.tau = 0, w1.tau =1))

## ----fig.width = 10, fig.height=6, warning = FALSE, message = FALSE-----------
g <- plot_weight_schemes(dfcount_gbsg)
print(g)


## ----fig.width = 10, fig.height=6---------------------------------------------
par(mfrow=c(1,2))
# Mine
# Set ymax a little above 1.0 to allow for logrank placement in topleft
# topleft is default; xmed.fraction = 0.65 positions the display of the median estimates below the HR estimates ("topright" legend [default])
# For details, please see summary of xmed.fraction in the Appendix below.

plot_weighted_km(dfcount=dfcount_gbsg, conf.int=FALSE, show.logrank = TRUE, 
                 put.legend.lr = "topleft", ymax = 1.05, xmed.fraction = 0.65)
title(main="GBSG (trial) data un-weighted")


# Compare with survfit
plot_km(df=df_gbsg, tte.name=tte.name, event.name=event.name, treat.name=treat.name)
title(main="GBSG (trial) data un-weighted
      via survfit")

## ----fig.width = 8, fig.height=6----------------------------------------------
plot_weighted_km(dfcount=dfcount_gbsg, conf.int=TRUE, show.logrank = TRUE, ymax = 1.05)
title(main="GBSG (trial) data un-weighted")


## -----------------------------------------------------------------------------
gbsg_validate <- within(rotterdam, {
  rfstime <- ifelse(recur == 1, rtime, dtime)
  t_months <- rfstime / 30.4375
  time_months <- t_months
  status <- pmax(recur, death)
  ignore <- (recur == 0 & death == 1 & rtime < dtime)
  status2 <- ifelse(recur == 1 | ignore, recur, death)
  rfstime2 <- ifelse(recur == 1 | ignore, rtime, dtime)
  time_months2 <- rfstime2 / 30.4375
  grade3 <- ifelse(grade == "3", 1, 0)
  treat <- hormon
  event <- status2
  tte <- time_months
  id <- as.numeric(1:nrow(rotterdam))
  SG0 <- ifelse(er <= 0, 0, 1)
})

## -----------------------------------------------------------------------------
df_rotterdam_np <- subset(gbsg_validate, nodes > 0)

## ----eval = has[["dplyr"]]----------------------------------------------------
library(dplyr)
rotterdam2 <- df_rotterdam_np %>%
  mutate(
    meno = factor(meno, levels = c(0, 1), labels = c("Pre", "Post")),
    grade3 = factor(grade3, levels = c(0, 1), labels = c("Not Grade 3", "Grade 3")),
    chemo = factor(chemo, levels = c(0, 1), labels = c("No", "Yes")),
    er_negative = factor(ifelse(er == 0, 1, 0), labels = c("ER negative","ER positive"))
  )


baseline_table <- create_baseline_table(
  data = rotterdam2,
  treat_var = "treat",
  vars_continuous = c("age", "pgr", "er", "nodes"),
  vars_categorical = c("size","meno","chemo","er_negative"),
  show_pvalue = TRUE,
  show_smd = TRUE
)



## ----results = 'asis', eval = has[["dplyr"]]----------------------------------
baseline_table

## ----ess-node-positive--------------------------------------------------------
pw_np <- ps_weights(treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er,
                    data = df_rotterdam_np, type = "stabilized")
df_rotterdam_np$sw.weights <- pw_np$weights

dfcount_np_wtd <- get_dfcounting(df=df_rotterdam_np, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, 
                                  weight.name="sw.weights")
knitr::kable(dfcount_np_wtd$ess, digits = c(0, 0, 1, 3))

## ----echo = FALSE-------------------------------------------------------------
# numbers quoted in the text
ess_np <- dfcount_np_wtd$ess
n_treat_np <- sum(df_rotterdam_np$treat == 1)
n_treat_restricted <- sum(df_rotterdam_np$treat == 1 & df_rotterdam_np$meno == 1 & df_rotterdam_np$chemo == 0)

## ----ps-weights-restricted----------------------------------------------------
df_rotterdam <- subset(df_rotterdam_np, meno == 1 & chemo == 0)

pw <- ps_weights(treat ~ age + size + grade3 + nodes + pgr + er,
                 data = df_rotterdam, type = "stabilized")
pw
df_rotterdam$sw.weights <- pw$weights

## ----balance-restricted-------------------------------------------------------
std_diff <- function(x, a, w) {
  (weighted.mean(x[a == 1], w[a == 1]) - weighted.mean(x[a == 0], w[a == 0])) /
    sqrt((var(x[a == 1]) + var(x[a == 0])) / 2)
}
X_ps <- model.matrix(~ age + size + grade3 + nodes + pgr + er, data = df_rotterdam)[, -1]
balance <- data.frame(
  covariate = colnames(X_ps),
  unweighted = apply(X_ps, 2, std_diff, a = df_rotterdam$treat, w = rep(1, nrow(df_rotterdam))),
  weighted = apply(X_ps, 2, std_diff, a = df_rotterdam$treat, w = df_rotterdam$sw.weights),
  row.names = NULL
)
knitr::kable(balance, digits = 3)

## -----------------------------------------------------------------------------
dfcount_rotterdam_unwtd <- get_dfcounting(df=df_rotterdam, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24)

dfcount_rotterdam_wtd <- get_dfcounting(df=df_rotterdam, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, 
                                  weight.name="sw.weights", weight_variance = "estimated", weight_model = pw)

knitr::kable(dfcount_rotterdam_wtd$ess, digits = c(0, 0, 1, 3))

## ----echo = FALSE-------------------------------------------------------------
# numbers quoted in the text: hazard ratios (FH(0,0), asymptotic intervals)
hr_text <- function(dfcount) {
  ci <- cox_rhogamma(dfcount = dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0))$hr_ci_asy
  paste0(fmt(ci$hr, 2), " (", fmt(ci$lower, 2), ", ", fmt(ci$upper, 2), ")")
}

## ----fig.width = 10, fig.height=6---------------------------------------------
par(mfrow=c(1,2))
plot_weighted_km(dfcount=dfcount_rotterdam_unwtd, xmed.fraction = 0.65)
title(main="Rotterdam un-weighted K-M curves")
plot_weighted_km(dfcount=dfcount_rotterdam_wtd, xmed.fraction = 0.65)
title(main="Rotterdam weighted K-M curves")


## ----fig.width = 10, fig.height=6---------------------------------------------
par(mfrow=c(1,2))
plot_weighted_km(dfcount=dfcount_gbsg, conf.int = TRUE, xmed.fraction = 0.65, show.logrank = TRUE, ymax = 1.05)
title(main="GBSG (trial) data un-weighted K-M curves")
plot_weighted_km(dfcount=dfcount_rotterdam_wtd, conf.int = TRUE, xmed.fraction = 0.65, show.logrank = TRUE, ymax = 1.05)
title(main="Rotterdam weighted K-M curves")

## ----rotterdam-band-resamples, fig.width = 10, fig.height = 6-----------------

par(mfrow=c(1,2))
 
 temp <- plotKM.band_subgroups(
    df=df_rotterdam,
    xlabel="Months", ylabel="Survival difference",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    tau_add=42, by.risk = 12, risk_delta=0.075, risk.pad=0.03, ymax.pad = 0.125,
    tte.name = tte.name, treat.name = treat.name, event.name = event.name, weight.name = "sw.weights",
    draws.band = 1000, show_resamples = TRUE,
    weight_variance = "estimated", weight_model = pw
  )
  legend("topleft", c("Difference estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.7)
  title(main="Rotterdam data
        propensity score weighted")


## ----rotterdam-gbsg-bands, fig.width = 10, fig.height = 10--------------------

par(mfrow=c(2,2))
 
 temp <- plotKM.band_subgroups(
    df=df_rotterdam,
    xlabel="Months", ylabel="Survival difference",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    tau_add=42, by.risk = 12, risk_delta=0.075, risk.pad=0.03, ymax.pad = 0.125,
    tte.name = tte.name, treat.name = treat.name, event.name = event.name, weight.name = "sw.weights",
    draws.band = 1000, show_resamples = FALSE,
    weight_variance = "estimated", weight_model = pw
  )
  legend("topleft", c("Difference estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.7)
  title(main="Rotterdam data
        propensity score weighted")

get_bands <- cumulative_rmst_bands(df = df_rotterdam, fit = temp$fit_itt, 
tte.name = tte.name, event.name = event.name, treat.name = treat.name , weight.name = "sw.weights", 
draws_sb = 1000, xlab="months", rmst_max_cex = 0.7,
weight_variance = "estimated", weight_model = pw)

legend("topleft", c("Cumulative RMST estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.7)

# kept for the comparison of the two weight_variance options below
band_rotterdam <- temp
rmst_rotterdam <- get_bands
  


temp <- plotKM.band_subgroups(
    df=df_gbsg, 
    Maxtau = NULL,
    xlabel="Months", ylabel="Survival difference",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    by.risk = 6, risk_delta=0.075, risk.pad=0.03,
    tte.name = tte.name, treat.name = treat.name, event.name = event.name, 
    draws.band = 1000, show_resamples = FALSE
  )
  legend("topleft", c("Difference estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.75)
  title(main="GBSG data")

  
get_bands <- cumulative_rmst_bands(df = df_gbsg, fit = temp$fit_itt, 
tte.name = tte.name, event.name = event.name, treat.name = treat.name, 
draws_sb = 1000, xlab="months", rmst_max_cex = 0.75)
  legend("topleft", c("Cumulative RMST estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.75)
  
  

## ----known-vs-estimated-------------------------------------------------------
dfcount_rotterdam_known <- get_dfcounting(df=df_rotterdam, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, 
                                  weight.name="sw.weights")

# survival difference and RMST with the weights treated as known, on the time grid of the band plot
fit_known <- KM_diff(df = df_rotterdam, tte.name = tte.name, event.name = event.name, treat.name = treat.name,
                     weight.name = "sw.weights", at_points = band_rotterdam$fit_itt$at_points, draws.band = 1000)
rmst_known <- cumulative_rmst_bands(df = df_rotterdam, fit = fit_known,
                                    tte.name = tte.name, event.name = event.name, treat.name = treat.name,
                                    weight.name = "sw.weights", draws_sb = 1000, plot = FALSE)

option_row <- function(dfcount, rmst) {
  cox <- cox_rhogamma(dfcount = dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0))
  ci <- cox$hr_ci_asy
  r <- rmst$rmst_maxtau_ci
  data.frame(
    weights = cox$weight_variance,
    z = fmt(dfcount$z.score, 2),
    HR = paste0(fmt(ci$hr, 2), " (", fmt(ci$lower, 2), ", ", fmt(ci$upper, 2), ")"),
    RMST = paste0(fmt(r[1], 1), " (", fmt(r[2], 1), ", ", fmt(r[3], 1), ")")
  )
}
known_vs_estimated <- rbind(option_row(dfcount_rotterdam_known, rmst_known),
                            option_row(dfcount_rotterdam_wtd, rmst_rotterdam))
knitr::kable(known_vs_estimated, align = "lrrr",
             col.names = c("Weights treated as", "Log-rank z", "Hazard ratio (95% CI)",
                           paste0("RMST difference at ", fmt(max(rmst_rotterdam$at_points), 1), " months (95% CI)")))

## ----eval = has[["dplyr"]]----------------------------------------------------

df_gbsg_treat <- subset(df_gbsg, treat == 1)
df_rotterdam_control <- subset(df_rotterdam_np, treat == 0)

df_gbsg_treat <- df_gbsg_treat %>%
  mutate(
    size = ifelse(size <= 20, "<=20", ifelse(size > 20 & size <=50, "20-50", ">50")),
  )


# note: GBSG does not have "chemo" variable

df_gbsg_treat <- df_gbsg_treat[,c("tte","event","treat","age","meno","size","grade3","nodes","pgr","er")]
df_rotterdam_control <- df_rotterdam_control[,c("tte","event","treat","age","meno","size","grade3","nodes","pgr","er")]

# cross-trial-comparison (CTC)

df_CTC <- rbind(df_gbsg_treat, df_rotterdam_control)

fit.ps <- glm(treat ~ age + meno + size + grade3 + nodes + pgr + er, data = df_CTC, family="binomial")
pihat <- fit.ps$fitted
pihat.null <- glm(treat ~ 1, family="binomial", data = df_CTC)$fitted
wt.1 <- pihat.null / pihat
wt.0 <- (1 - pihat.null) / (1 - pihat)

df_CTC$sw.weights <- ifelse(df_CTC$treat == 1, wt.1, wt.0)

## ----eval = has[["dplyr"]]----------------------------------------------------

df_CTC$er_negative <- ifelse(df_CTC$er == 0, 1, 0) 

baseline_table <- create_baseline_table(
  data = df_CTC,
  treat_var = "treat",
  vars_continuous = c("age", "pgr", "er", "nodes"),
  vars_categorical = c("size","meno","er_negative"),
  show_pvalue = TRUE,
  show_smd = TRUE
)



## ----results = 'asis', eval = has[["dplyr"]]----------------------------------
baseline_table

## ----eval = has[["dplyr"]]----------------------------------------------------
dfcount_CTC_unwtd <- get_dfcounting(df = df_CTC, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24)

dfcount_CTC_wtd <- get_dfcounting(df=df_CTC, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, 
                                  weight.name="sw.weights", weight_variance = "known")

knitr::kable(dfcount_CTC_wtd$ess, digits = c(0, 0, 1, 3))

## ----echo = FALSE, eval = has[["dplyr"]]--------------------------------------
# numbers quoted in the text
ess_CTC <- dfcount_CTC_wtd$ess
gbsg_treated <- subset(gbsg, hormon == 1)
rott_control <- subset(rotterdam, nodes > 0 & hormon == 0)
pct <- function(x) paste0(fmt(100 * mean(x), 0), "%")
treated_at_risk <- function(t) sum(df_CTC$tte[df_CTC$treat == 1] >= t)
horizon <- 60

## ----fig.width = 10, fig.height=6, eval = has[["dplyr"]]----------------------
par(mfrow=c(2,2))
plot_weighted_km(dfcount=dfcount_gbsg, conf.int = TRUE, xmed.fraction = 0.65, show.logrank = TRUE, ymax = 1.2, show.med = FALSE)
title(main="GBSG (trial) data un-weighted K-M curves")
plot_weighted_km(dfcount=dfcount_CTC_wtd, conf.int = TRUE, xmed.fraction = 0.65, show.logrank = TRUE, ymax = 1.2, show.med = FALSE)
title(main="CTC weighted K-M curves")
plot_weighted_km(dfcount=dfcount_CTC_unwtd, conf.int = TRUE, xmed.fraction = 0.65, show.logrank = TRUE, ymax = 1.2, show.med = FALSE)
title(main="CTC un-weighted K-M curves")

## ----echo = FALSE, eval = has[["dplyr"]]--------------------------------------
# cumulative log-rank z; the curve is drawn from the 10th event in each arm
cum_lr <- wlr_cumulative(dfcount_CTC_wtd, scheme_params = list(rho = 0, gamma = 0), scheme = "fh", return_cumulative = TRUE)
event_times <- function(a) sort(df_CTC$tte[df_CTC$treat == a & df_CTC$event == 1])
t_start <- max(event_times(1)[10], event_times(0)[10])
z_at <- function(t) cum_lr$z.score[findInterval(t, cum_lr$time)]
last_treated <- max(df_CTC$tte[df_CTC$treat == 1])

## ----cumulative-logrank, fig.width = 8, fig.height=6, fig.cap = if (has[["dplyr"]]) paste0("Cumulative weighted log-rank z-statistic (cross-trial comparison), drawn from the 10th event in each arm (", fmt(t_start, 1), " months). Dashed line: ", horizon, " months."), eval = has[["dplyr"]]----
par(mfrow=c(1,1))
temp <- wlr_cumulative(dfcount_CTC_wtd, scheme_params = list(rho = 0, gamma = 0), scheme = "fh", return_cumulative = TRUE)
shown <- temp$time >= t_start
# y-axis: the range of the plotted curve and of the reference line
ylim <- range(c(temp$z.score[shown], qnorm(0.975)))
ylim <- ylim + c(-0.05, 0.05) * diff(ylim)
plot(temp$time[shown], temp$z.score[shown], xlab="time", ylab="cumulative log-rank Z", type="s", ylim=ylim)
abline(h = qnorm(0.975), lwd=2, col="red")
abline(v = horizon, lwd=1, col="blue", lty=2)


## ----ctc-horizon, eval = has[["dplyr"]]---------------------------------------
# follow-up cut at the horizon
df_CTC_cut <- df_CTC
df_CTC_cut$event <- ifelse(df_CTC_cut$tte > horizon, 0, df_CTC_cut$event)
df_CTC_cut$tte <- pmin(df_CTC_cut$tte, horizon)
dfcount_CTC_cut <- get_dfcounting(df=df_CTC_cut, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, 
                                  weight.name="sw.weights", weight_variance = "known")

ctc_row <- function(dfcount, followup) {
  cox <- cox_rhogamma(dfcount = dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0))
  ci <- cox$hr_ci_asy
  data.frame(followup = followup, weights = cox$weight_variance, z = fmt(dfcount$z.score, 2),
             HR = paste0(fmt(ci$hr, 2), " (", fmt(ci$lower, 2), ", ", fmt(ci$upper, 2), ")"))
}
knitr::kable(rbind(ctc_row(dfcount_CTC_wtd, "all"), ctc_row(dfcount_CTC_cut, paste0("cut at ", horizon, " months"))),
             align = "llrr", col.names = c("Follow-up", "Weights treated as", "Log-rank z", "Hazard ratio (95% CI)"))

## ----ctc-bands, fig.width = 10, fig.height = 10, eval = has[["dplyr"]]--------
draws <- 1000
par(mfrow=c(2,2))
 
temp <- plotKM.band_subgroups(
    df=df_gbsg, 
    Maxtau = NULL,
    xlabel="Months", ylabel="Survival difference",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    by.risk = 6, risk_delta=0.075, risk.pad=0.03,
    tte.name = tte.name, treat.name = treat.name, event.name = event.name, 
    draws.band = draws, show_resamples = FALSE
  )
  legend("topleft", c("Difference estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.75)
  title(main="GBSG data")

  
get_bands <- cumulative_rmst_bands(df = df_gbsg, fit = temp$fit_itt, 
tte.name = tte.name, event.name = event.name, treat.name = treat.name, 
draws_sb = draws, xlab="months", rmst_max_cex = 0.75)
  legend("topleft", c("Cumulative RMST estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.75)



 temp <- plotKM.band_subgroups(
    df=df_CTC,
    xlabel="Months", ylabel="Survival difference",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    tau_add=42, by.risk = 12, risk_delta=0.075, risk.pad=0.03, ymax.pad = 0.125,
    tte.name = tte.name, treat.name = treat.name, event.name = event.name, weight.name ="sw.weights",
    draws.band = draws, show_resamples = FALSE, weight_variance = "known"
  )
  legend("topleft", c("Difference estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.7)
  title(main="CTC data")

get_bands <- cumulative_rmst_bands(df = df_CTC, fit = temp$fit_itt, 
tte.name = tte.name, event.name = event.name, treat.name = treat.name , weight.name = "sw.weights", 
draws_sb = draws, xlab="months", rmst_max_cex = 0.7, weight_variance = "known")

legend("topleft", c("Cumulative RMST estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.7)

# kept for the reading at the horizon
band_CTC <- temp
rmst_CTC <- get_bands


## ----echo = FALSE, eval = has[["dplyr"]]--------------------------------------
# cross-trial estimates at the horizon (pointwise 95% limits), from the fits plotted above
at_h <- function(x, at) x[findInterval(horizon, at)]
f <- band_CTC$fit_itt
diff_h <- paste0(fmt(at_h(f$dhat, f$at_points), 2), " (", fmt(at_h(f$lower, f$at_points), 2), ", ", fmt(at_h(f$upper, f$at_points), 2), ")")
rmst_h <- paste0(fmt(at_h(rmst_CTC$rmst_time, rmst_CTC$at_points), 1), " (", fmt(at_h(rmst_CTC$rmst_time_lower, rmst_CTC$at_points), 1), ", ",
                 fmt(at_h(rmst_CTC$rmst_time_upper, rmst_CTC$at_points), 1), ")")

## ----fig.width = 8, fig.height=6----------------------------------------------
par(mfrow=c(1,1))
sg_labs <- c("er == 0","er > 0", "meno == 0", "meno ==1")
sg_cols <- c("blue", "brown", "green", "turquoise")

# Randomly sample colors from the built-in palette
# n_sg <- length(sg_labs)
# set.seed(123) # for reproducibility
# sg_cols <- sample(colors(), n_sg)

  temp <- plotKM.band_subgroups(
    df=df_gbsg,
    sg_labels = sg_labs,
    sg_colors = sg_cols,
    xlabel="Months", ylabel="Survival difference",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    tau_add=42, by.risk=6, risk_delta=0.05, risk.pad=0.03, draws.band = 1000,
    tte.name = tte.name, treat.name = treat.name, event.name = event.name 
  )
  legend("topleft", c("ITT", sg_labs),
         lwd=2, col=c("black", sg_cols), bty="n", cex=0.75)
  title(main="ITT and subgroups")


## ----fig.width = 12, fig.height=6---------------------------------------------
par(mfrow=c(1,2))
  temp <- plotKM.band_subgroups(
    df = subset(df_gbsg, er > 0),
    xlabel="Months", ylabel="Survival difference",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    by.risk=6, risk_delta=0.075, risk.pad=0.03,
    tte.name = tte.name, treat.name = treat.name, event.name = event.name, draws.band = 1000, show_resamples = TRUE
  )
  legend("topleft", c("Difference estimate", "95% pointwise CIs"),
         lwd=c(2,1), col=c("black"), lty=c(1,2), bty="n", cex=0.75)
  title(main="Estrogen receptor positive (er>0) sub-population")
  
  

## ----message = FALSE, warning = FALSE, eval = has[["data.table"]] && has[["gt"]]----
library(data.table)
library(gt)

# Restricted Rotterdam population, propensity-score weighted:
#   dfcount_rotterdam_wtd    weights treated as estimated (weight_variance = "estimated")
#   dfcount_rotterdam_known  weights treated as known (the default)
# cox_rhogamma() takes the option from the df_counting() result it is given.
#
# Note: with subject weights (here propensity-score weights) that are not all equal, the asymptotic
# SE of the log hazard ratio is the score-residual (robust) SE, and the resampling SE (draws > 0) is
# built from the same per-subject score contributions, so the two estimate the same quantity.
# hr_ci_asy is centred at the estimate; hr_ci_star is centred at the de-biased estimate and uses
# the resampling SE.

cox_row <- function(analysis, dfcount, interval = c("asymptotic", "resampling"), ...) {
  interval <- match.arg(interval)
  get_rg <- cox_rhogamma(dfcount = dfcount, ..., draws = if (interval == "resampling") 1000 else 0, verbose = FALSE)
  ci <- if (interval == "resampling") get_rg$hr_ci_star else get_rg$hr_ci_asy
  data.table(analysis = analysis, weights = get_rg$weight_variance, interval = interval,
             z = get_rg$z.score, hr = get_rg$hr_ci_asy$hr, lower = ci$lower, upper = ci$upper)
}

# compare with coxph() 
#cat("Comparison with coxph with weighting, coxph:", dfcount_rotterdam_wtd$cox_results$cox_text, "\n")

res_update <- rbind(
  cox_row("FH(0,0)", dfcount_rotterdam_wtd, "asymptotic", scheme = "fh", scheme_params = list(rho = 0, gamma = 0)),
  cox_row("FH(0,0)", dfcount_rotterdam_known, "resampling", scheme = "fh", scheme_params = list(rho = 0, gamma = 0)),
  cox_row("FH(0,1)", dfcount_rotterdam_known, "resampling", scheme = "fh", scheme_params = list(rho = 0, gamma = 1)),
  cox_row("FH_exp2", dfcount_rotterdam_known, "resampling", scheme = "fh_exp2"),
  cox_row("MB(12)", dfcount_rotterdam_known, "resampling", scheme = "MB", scheme_params = list(mb_tstar = 12)),
  cox_row("zero_one(6)", dfcount_rotterdam_known, "resampling", scheme = "custom_time", scheme_params = list(t.tau = 6, w0.tau = 0, w1.tau = 1))
)

res_update %>% gt()  %>% fmt_number(decimals = 3) 


## -----------------------------------------------------------------------------
dfcount_rotterdam_known <- get_dfcounting(df=df_rotterdam, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, 
                                  weight.name="sw.weights", check.seKM = TRUE, draws = 0)

## ----echo = FALSE-------------------------------------------------------------
# largest absolute difference from survfit's weighted standard error, at survfit's event times
se_diff_survfit <- function(a) {
  sf <- summary(survfit(Surv(tte, event) ~ 1, data = df_rotterdam[df_rotterdam$treat == a, ], weights = sw.weights))
  v <- if (a == 1) dfcount_rotterdam_known$sig2_surv1 else dfcount_rotterdam_known$sig2_surv0
  max(abs(sqrt(v[findInterval(sf$time, dfcount_rotterdam_known$at_points)]) - sf$std.err))
}

## -----------------------------------------------------------------------------
dfcount_rotterdam_draws <- get_dfcounting(df=df_rotterdam, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, 
                                  weight.name="sw.weights", check.seKM = TRUE, draws = 5000)

## ----eval = has[["adjustedCurves"]] && has[["pammtools"]]---------------------
# Compare with adjustedCurves
library(adjustedCurves)
library(pammtools)

df_rotterdam$hormon2 <- with(df_rotterdam, factor(hormon, labels=c("No","Yes")))

ps_model <- glm(hormon2 ~ age+size+grade3+nodes+pgr+er, data = df_rotterdam, family="binomial")

iptw <- adjustedsurv(data=df_rotterdam, variable="hormon2", ev_time="tte", event = "event", method = "iptw_km", treatment_model = ps_model, conf_int = TRUE)


## ----fig.width = 8, fig.height = 6--------------------------------------------
plot_weighted_km(dfcount=dfcount_rotterdam_known, conf.int = TRUE)
 title(main="Rotterdam weighted K-M curves")

## ----fig.width = 8, fig.height = 6, eval = has[["adjustedCurves"]] && has[["pammtools"]]----
 plot(iptw, conf_int = TRUE, legend.title = "Hormonal therapy")

## ----fig.width = 12, fig.height = 6, eval = has[["adjustedCurves"]] && has[["pammtools"]]----
# Check SEs
par(mfrow=c(1,2))

df_mine <- dfcount_rotterdam_known

df_check <- subset(iptw$adj, group == "No")
yymax <- max(c(sqrt(df_mine$sig2_surv0[df_mine$idv0.check]), df_check$se))
toget <- df_mine$idv0.check
tt <- df_mine$at_points[toget]
yy <- sqrt(df_mine$sig2_surv0)[toget]
plot(tt,yy, type="s", lty=1, col="lightgrey", lwd=4, ylim=c(0,yymax), xlab="time", ylab="SEs")
with(df_check, lines(time, se, type="s", lty=2, lwd=1, col="red"))
title(main="Control SEs")


df_check <- subset(iptw$adj, group == "Yes")
yymax <- max(c(sqrt(df_mine$sig2_surv1[df_mine$idv1.check]), df_check$se))
toget <- df_mine$idv1.check
tt <- df_mine$at_points[toget]
yy <- sqrt(df_mine$sig2_surv1)[toget]
plot(tt,yy, type="s", lty=1, col="lightgrey", lwd=4, ylim=c(0,yymax), xlab="time", ylab="SEs")
with(df_check, lines(time, se, type="s", lty=2, lwd=1, col="red"))
title(main = "Treat SEs")


## ----echo = FALSE, eval = has[["adjustedCurves"]] && has[["pammtools"]]-------
# weightedsurv (grey) against adjustedCurves iptw_km (red) at the times adjustedCurves reports
vs_adjustedCurves <- function(group, surv, sig2) {
  ad <- subset(iptw$adj, group == group & time > 0)
  ad <- ad[ad$group == group, ]
  i <- findInterval(ad$time, dfcount_rotterdam_known$at_points)
  ok <- ad$se > 0 & sig2[i] > 0
  list(surv_diff = max(abs(surv[i] - ad$surv)), se_ratio = range(sqrt(sig2[i][ok]) / ad$se[ok]))
}
ac0 <- vs_adjustedCurves("No", dfcount_rotterdam_known$surv0, dfcount_rotterdam_known$sig2_surv0)
ac1 <- vs_adjustedCurves("Yes", dfcount_rotterdam_known$surv1, dfcount_rotterdam_known$sig2_surv1)

## -----------------------------------------------------------------------------
# dfcount_gbsg (un-weighted) contains logrank (chisq) stats from specified scheme (default is standard logrank rho=0 and gamma=0) from 3 sources:
# (1) lr^2 / sig2_lr: the weighted log-rank statistic and its variance (wlr_cumulative)
# (2) logrank_results$chisq: without subject weights (as here), survdiff, which only provides rho options
#     so gamma != 0 is not available; with subject weights that are not all equal it is instead the
#     robust score test of the weighted Cox model (coxph, robust = TRUE), see logrank_results$method
# (3) z.score^2: the standardized weighted log-rank statistic (Cox score-test perspective)
check_results(dfcount_gbsg)
# Check wlr_dhat_estimates which calculates from any scheme based on counting dataset 
# In addition, calculates dhat's at tzero as well as correlation with log-rank
temp <- wlr_dhat_estimates(dfcounting = dfcount_gbsg, scheme_params = list(rho = 0, gamma = 0), scheme = "fh", tzero = 12)
temp$z.score^2


## -----------------------------------------------------------------------------
res <- df_counting(df=df_gbsg, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, strata.name="grade")

## -----------------------------------------------------------------------------
with(res,lr_stratified^2/sig2_lr_stratified)
# logrank_results: survdiff, stratified (un-weighted data; with subject weights that are not all equal
# this element is the stratified robust score test of the weighted Cox model)
res$logrank_results$chisq
 # Cox score stratified
with(res,z.score_stratified^2)

## -----------------------------------------------------------------------------
with(res,lr^2/sig2_lr)
with(res,z.score^2)

# Restore when done
par(oldpar)

