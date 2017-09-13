


#' VarianceEM
#' @description The function VarianceEM performs variance estimation for the Multistate Cure Model when the model is fit using the EM Algorithm. Variances are estimated by fitting the multistate cure model (via EM) to bootstrap samples of the data.
#' @param fit Multistate cure model fit from MultiCure
#' @param iternum Number of iterations of the EM algorithm to apply to each bootstrap sample
#' @param bootnum Number of bootstrap samples
#' @param datWIDE A data frame with the following columns: 
#' \itemize{
#' \item Y_R, the recurrence event/censoring time
#' \item delta_R, the recurrence event/censoring indicator
#'\item Y_D, the death event/censoring time 
#' \item delta_D, the death event/censoring indicator
#' \item G, the cure status variable. This takes value 1 for known non-cured, 0 for "known" cured and NA for unknown cure status
#'}
#' @param Cov  A matrix containing the covariates
#' @param ASSUME This variables indicates what equality assumptions we are making regarding the 24 and 14 transitions. The possible options are:
#' \itemize{
#' \item 'SameHazard': Lambda_14(t) = Lambda_24(t)
#' \item 'AllSeparate': No restrictions on Lambda_14(t) and Lambda_24(t)
#' \item  'ProportionalHazard': Lambda_14(t) = Lambda_24(t) exp(Beta0)
#' \item  'SameBaseHaz': Lambda^0_14(t) = Lambda^0_24(t), No restrictions on beta_14 and beta_24
#' }
#' @param TransCov a list with elements: Trans13, Trans24, Trans14, Trans34, PNonCure. Each list element is a vector containing the names of the variables in Cov to be used in the model for the corresponding transition. 13 is NonCured -> Recurrence, 24 is Cured -> Death, 14 is NonCured -> Death, 34 is Recurrence -> Death. PNonCure contains the names of the covariates for the logistic regression for P(NonCure). 
#' @param BASELINE This variable indicates the assumptions about the baseline hazard form. This can take values 'weib' and 'cox'
#' @param PENALTY This variable indicates whether we are using any variable selection in the model fitting. Right now, the options are 'None' (no variable selection), 'Ridge' (ridge regression for all covariates in all models) and 'Lasso' (lasso for all covariates in all models, only implemented for Cox baseline hazards)
#'
#' @return OUT a matrix containing the following:
#' \itemize{
#' \item Estimate an estimate of the multistate cure model parameter
#' \item Variance an estimate of variance of the the multistate cure model parameter
#'}
#' @examples
#' attach(SimulateMultiCure(type = "NoMissingness"))
#' Cov = data.frame(X1,X2)
#' VARS = names(Cov)
#' TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#' datWIDE = data.frame( Y_R, Y_D, delta_R , delta_D, G)
#' fit = MultiCure(iternum=100, datWIDE, Cov, ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "weib") 
#' OUT = VarianceEM(fit,iternum=20, bootnum=50, datWIDE, Cov, ASSUME = "SameHazard", TransCov, BASELINE = "weib")
#' @export

VarianceEM = function(fit,iternum, bootnum, datWIDE, Cov, ASSUME, TransCov, BASELINE, PENALTY = 'None'){				
	if(BASELINE == 'weib'){
		PARAMINIT = list(beta = fit[[1]], alpha = fit[[2]], scale = fit[[3]], shape = fit[[4]])		
	}else{
		PARAMINIT = list(beta = fit[[1]], alpha = fit[[2]])		
	}
	if(BASELINE == 'weib'){
		WRAPPERFUNC2 = function(data, indices){
			datWIDETEMP = as.data.frame(data[indices,1:5])
			CovImpTEMP = as.data.frame(data[indices,6:length(data[1,])])
			fitTEMP = MultiCure(iternum,datWIDE = datWIDETEMP,Cov=CovImpTEMP, ASSUME = ASSUME, 
				TransCov= TransCov, BASELINE = BASELINE, PENALTY = PENALTY, PARAMINIT = PARAMINIT)
			return(c(fitTEMP[[1]], fitTEMP[[2]], fitTEMP[[3]], fitTEMP[[4]]))
		}	
		WRAPPERFUNC = tryCatch(WRAPPERFUNC2,error = function(cond){return(rep(NA, P))},warning = function(cond){return(rep(NA, P))})	
	}else{
		WRAPPERFUNC2 = function(data, indices){
			datWIDETEMP = as.data.frame(data[indices,1:5])
			CovImpTEMP = as.data.frame(data[indices,6:length(data[1,])])
			fitTEMP = MultiCure(iternum,datWIDE = datWIDETEMP,Cov=CovImpTEMP, ASSUME = ASSUME,
				TransCov= TransCov, BASELINE = BASELINE, PENALTY = PENALTY, PARAMINIT = PARAMINIT)
			return(c(fitTEMP[[1]], fitTEMP[[2]]))
		}				
		WRAPPERFUNC = tryCatch(WRAPPERFUNC2,error = function(cond){return(rep(NA, P))},warning = function(cond){return(rep(NA, P))})	
	}
	datWIDETEMP = subset(datWIDE, select = c(Y_R, Y_D, delta_R, delta_D, G))
	#sum(datWIDETEMP$Y_R<datWIDETEMP$Y_D & datWIDETEMP$delta_R == 0 & is.na(datWIDETEMP$G))
	DAT = data.frame(datWIDETEMP, Cov)
	BOOT = boot::boot(data = DAT, statistic = WRAPPERFUNC, R = bootnum)
	SDs = apply(BOOT$t,2,sd, na.rm=T)
	SAVE_VAR = (SDs)^2

	BOOT$t_sub = ifelse(abs(BOOT$t)>5 | is.na(BOOT$t), matrix(NA,ncol = ncol(BOOT$t), nrow = nrow(BOOT$t)), BOOT$t)
	SDs_sub = apply(BOOT$t_sub,2,sd, na.rm=T)
	SAVE_VAR_sub = (SDs_sub)^2
	if(ASSUME == 'ProportionalHazard'){
		TransCov$Trans14 = c(TransCov$Trans14, 'INT')
	}
	if(BASELINE == 'weib')
	{
		PARAM = c(fit[[1]], fit[[2]], fit[[3]], fit[[4]])	
		VARNAMES = c(paste(TransCov$Trans13, '_13', sep = ''),paste(TransCov$Trans24, '_24', sep = ''),
					paste(TransCov$Trans14, '_14', sep = ''), paste(TransCov$Trans34, '_34', sep = ''),
					paste(c('Intercept',TransCov$PNonCure), '_NonCure', sep = ''),
					paste('Scale', c('_13', '_24', '_14', '_34'), sep = ''), 
					paste('Shape', c('_13', '_24', '_14', '_34'), sep = ''))		
	}else if(BASELINE == 'cox'){
		PARAM = c(fit[[1]], fit[[2]])	
		VARNAMES = c(paste(TransCov$Trans13, '_13', sep = ''),paste(TransCov$Trans24, '_24', sep = ''),
					paste(TransCov$Trans14, '_14', sep = ''), paste(TransCov$Trans34, '_34', sep = ''),
					paste(c('Intercept',TransCov$PNonCure), '_NonCure', sep = ''))		
	}

	return(data.frame(VARNAMES,Estimate = PARAM, Variance = SAVE_VAR_sub))
}
