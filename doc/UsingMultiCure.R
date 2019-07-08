## ---- echo=FALSE, fig.cap="Figure 1: Multistate Cure Model Structure", out.width = '80%'----
knitr::include_graphics("imagesUsingMultiCure/ModelDiagram.png")
library(MultiCure)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
NONE = SimulateMultiCure(type = 'NoMissingness') 
COV = SimulateMultiCure(type = 'CovariateMissingness')
CENS = SimulateMultiCure(type = 'UnequalCensoring')

## ---- echo = FALSE, eval = TRUE,  fig.width = 7, fig.height= 4-----------
library(survival)
par(mfrow = c(1,2))
plot(survfit(Surv(NONE$Y_R, NONE$delta_R)~1), mark.time = T, main = 'Time to Recurrence', xlab = 'Years', ylab = 'Event-Free Probability', cex.main = 0.8)
plot(survfit(Surv(NONE$Y_D, NONE$delta_D)~1), mark.time = T, main = 'Overall Survival', xlab = 'Years', ylab = 'Event-Free Probability', cex.main = 0.8)

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  ### Prepare Data
#  Cov = data.frame(X1 = NONE$X1,X2 = NONE$X2)
#  VARS = names(Cov)
#  TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#  datWIDE = data.frame( Y_R = NONE$Y_R, Y_D = NONE$Y_D, delta_R = NONE$delta_R ,
#            delta_D = NONE$delta_D, G = NONE$G)
#  
#  ### Fit Model
#  fit = MultiCure(iternum = 50, datWIDE, Cov, ASSUME = 'SameHazard', TransCov=TransCov,
#            BASELINE = 'weib')
#  OUT = VarianceEM(fit,iternum=20, bootnum=50, datWIDE, Cov, ASSUME = 'SameHazard', TransCov=TransCov,
#            BASELINE = 'weib')

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  ### Prepare Data
#  Cov = data.frame(X1 = COV$X1, X2 = COV$X2)
#  VARS = names(Cov)
#  TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#  datWIDE = data.frame( Y_R = COV$Y_R, Y_D = COV$Y_D, delta_R = COV$delta_R ,
#            delta_D = COV$delta_D, G = COV$G)
#  
#  ### Obtain Point Estimates
#  fit = MultiCure(iternum = 200, datWIDE, Cov, COVIMPUTEFUNCTION_Example, COVIMPUTEINITIALIZE_Example,
#            IMPNUM = 10,ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'weib')
#  beta = apply(fit[[5]][,190:200], 1, mean)
#  
#  ### Variance Estimation
#  OUT = VarianceMCEM(fit,var_method = 'default', datWIDE = datWIDE,  ASSUME = 'SameHazard',
#            TransCov = TransCov, BASELINE = 'weib', COVIMPUTEFUNCTION = COVIMPUTEFUNCTION_Example,
#            COVIMPUTEINITIALIZE = COVIMPUTEINITIALIZE_Example, POSTITER = 5)	

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  ### Prepare Data
#  Cov = data.frame(X1 = CENS$X1, X2 = CENS$X2)
#  VARS = names(Cov)
#  TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#  datWIDE = data.frame( Y_R = CENS$Y_R, Y_D = CENS$Y_D, delta_R = CENS$delta_R,
#            delta_D = CENS$delta_D, G = CENS$G)
#  
#  ### Obtain Point Estimates
#  fit = MultiCure(iternum = 200, datWIDE  = datWIDE, Cov = Cov, IMPNUM = 10,
#            ASSUME = 'SameHazard', TransCov = TransCov, BASELINE = 'weib',
#            UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBREJECTION)
#  
#  ### Variance Estimation
#  OUT = VarianceMCEM(fit,var_method = 'default', datWIDE = datWIDE,  ASSUME = 'SameHazard',
#            TransCov = TransCov, BASELINE = 'weib', UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBREJECTION,
#            POSTITER = 5)

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  STATEOCCUPANCYWEIB(times = seq(0,max(datWIDE$Y_D),1), TransCov,
#            newCov = data.frame(X1 = c(0,0.5), X2 = c(0,0.5)), beta = fit[[1]], alpha = fit[[2]],
#            scale = fit[[3]], shape = fit[[4]])

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  Haz = BaselineHazard_NOIMP(datWIDE, Cov, beta = fit[[1]], alpha = fit[[2]], TransCov, ASSUME = 'SameHazard', p = fit[[5]])
#  STATEOCCUPANCYCOX_NOIMP(times = seq(0,max(datWIDE$Y_D),1), TransCov,
#            newCov = data.frame(X1 = c(0,0.5), X2 = c(0,0.5)), beta = fit[[1]], alpha = fit[[2]],
#            Haz_13 = Haz[[1]], Haz_24 = Haz[[2]], Haz_14 = Haz[[3]], Haz_34 = Haz[[4]])

