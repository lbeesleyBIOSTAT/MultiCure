
#' VarianceMCEM_BOOT
#' @description The function VarianceMCEM_BOOT performs variance estimation for the Multistate Cure Model when the model is fit using a Monte Carlo EM Algorithm. For each imputed dataset, this function estimates parameter variances by fitting the multistate cure model to bootstrap samples of the imputed dataset. Then, Rubin's rules are used to obtain a single set of parameter estimates and standard errors across imputed datasets. Important!!!: The function ProperDraws_MC should be used to obtain the imputed datasets to be used in this function. This ensures that the imputed datasets are (roughly) proper imputations and that Rubin's rules can then be applied. 
#'
#' @param fit Multistate cure model fit from MultiCure
#' @param bootnum Number of bootstrap samples used for each imputed dataset
#' @param datWIDE A data frame with the following columns: 
#' \itemize{
#' \item Y_R, the recurrence event/censoring time
#' \item delta_R, the recurrence event/censoring indicator
#'\item Y_D, the death event/censoring time 
#' \item delta_D, the death event/censoring indicator
#' \item G, the cure status variable. This takes value 1 for known non-cured, 0 for "known" cured and NA for unknown cure status
#'}
#' @param CovImp  A list with IMPNUM elements containing the imputations of Cov output from the function ProperDraws_MC
#' @param GImp  A matrix with IMPNUM elements containing the imputations of G output from the function ProperDraws_MC
#' @param YRImp  A matrix with IMPNUM elements containing the imputations of Y_R output from the function ProperDraws_MC
#' @param deltaRImp  A matrix with IMPNUM elements containing the imputations of delta_R output from the function ProperDraws_MC
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
#' \item Estimate an estimate of the multistate cure model parameter from Rubin's Rules
#' \item Variance an estimate of variance of the the multistate cure model parameter from Rubin's Rules
#' \item v The estimated degrees of freedom of the t-distribution of the parameter estimate from Rubin's Rules
#'}
#' @details This function provides parameter estimates and estimated variances. The parameter estimates are obtained using Rubin's rules, but an alternative estimate of the multistate cure model parameter can be obtained by averaging the parameter estimates from the last few iterations of the model fitting algorithm. In our experience, we found that the approach that averages across the last few iterations (rather than estimated using Rubin's rules) provides a better estimate of the parameter of interest.
#'
#' @examples
#' attach(SimulateMultiCure(type = "UnequalCensoring"))
#' Cov = data.frame(X1,X2)
#' VARS = names(Cov)
#' TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#' datWIDE = data.frame( Y_R, Y_D, delta_R , delta_D, G)
#' fit = MultiCure(iternum = 100, datWIDE, Cov, ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "weib", IMPNUM = 10) ### Note: This will take a moment
#' Proper = ProperDraws_MC(datWIDE,Cov, CovImp = fit[[9]], GImp = fit[[10]], YRImp = fit[[11]], deltaRImp = fit[[12]], ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "weib") ### Note: This will take a moment
#' OUT = VarianceMCEM_BOOT(fit, bootnum = 50, datWIDE, CovImp = Proper[[1]], GImp = Proper[[2]], YRImp = Proper[[3]], deltaRImp = Proper[[4]],  ASSUME = "SameHazard", TransCov, BASELINE = "weib")
#' @export

VarianceMCEM_BOOT = function(fit,bootnum, datWIDE, CovImp, GImp, YRImp, deltaRImp, ASSUME, TransCov, BASELINE, PENALTY = 'None'){			
	
	########################################
	### Initializing, Defining Functions ###
	########################################
	
	ASSUME = match.arg(ASSUME, choices = c('SameHazard', 'AllSeparate', 'SameBaseHaz', 'ProportionalHazard'))
	BASELINE = match.arg(BASELINE, choices = c('weib','cox'))
	PENALTY = match.arg(PENALTY, choices = c('None', 'Ridge', 'Lasso'))

	if(BASELINE == 'weib'){
		PARAMINIT = list(beta = fit[[1]], alpha = fit[[2]], scale = fit[[3]], shape = fit[[4]])		
	}else{
		PARAMINIT = list(beta = fit[[1]], alpha = fit[[2]])		
	}
	if(BASELINE == 'weib'){
		WRAPPERFUNC2 = function(data, indices){
			datWIDETEMP = as.data.frame(data[indices,1:5])
			CovImpTEMP = as.data.frame(data[indices,6:length(data[1,])])
			fitTEMP = MultiCure(iternum=1,datWIDE = datWIDETEMP,Cov=CovImpTEMP, ASSUME = ASSUME, 
				TransCov= TransCov, BASELINE = BASELINE, PENALTY = PENALTY, PARAMINIT = PARAMINIT)
			return(c(fitTEMP[[1]], fitTEMP[[2]], fitTEMP[[3]], fitTEMP[[4]]))
		}		
		WRAPPERFUNC = tryCatch(WRAPPERFUNC2,error = function(cond){return(rep(NA,P))},warning = function(cond){return(rep(NA,P))})	
	}else{
		WRAPPERFUNC2 = function(data, indices){
			datWIDETEMP = as.data.frame(data[indices,1:5])
			CovImpTEMP = as.data.frame(data[indices,6:length(data[1,])])
			fitTEMP = MultiCure(iternum=1,datWIDE = datWIDETEMP,Cov=CovImpTEMP, ASSUME = ASSUME,
				TransCov= TransCov, BASELINE = BASELINE, PENALTY = PENALTY, PARAMINIT = PARAMINIT)
			return(c(fitTEMP[[1]], fitTEMP[[2]]))
		}	
		WRAPPERFUNC = tryCatch(WRAPPERFUNC2,error = function(cond){return(rep(NA,P))},warning = function(cond){return(rep(NA,P))})				
	}
	IMPNUM = length(CovImp)
	SAVE_VAR = c()
	SAVE_PARAM = c()
	
	##########################################
	### Fitting Model to Bootstrap Samples ###
	##########################################
	
	for(i in 1:IMPNUM){
		CovImpTEMP = CovImp[[i]]
		datWIDETEMP = datWIDE
		datWIDETEMP$G = GImp[,i]
		datWIDETEMP$delta_R = deltaRImp[,i]
		datWIDETEMP$Y_R = YRImp[,i]
		datWIDETEMP = subset(datWIDETEMP, select = c(Y_R, Y_D, delta_R, delta_D, G))
		DAT = data.frame(datWIDETEMP, CovImpTEMP)
		BOOT = boot::boot(data = DAT, statistic = WRAPPERFUNC, R = bootnum)
			### Guarding against non-converging bootstrap samples
		BOOT$t = ifelse(abs(BOOT$t)>5, matrix(NA,ncol = ncol(BOOT$t), nrow = nrow(BOOT$t)), BOOT$t)
		SDs = apply(BOOT$t,2,sd, na.rm=T)
		MEANS = apply(BOOT$t,2,mean, na.rm=T)
		SAVE_PARAM = rbind(SAVE_PARAM, MEANS)
		SAVE_VAR = rbind(SAVE_VAR,(SDs)^2)
		print(paste('Finished Imputation', i, sep = ' '))
	}
	
	###########################
	### Apply Rubin's Rules ###
	###########################
	OUT = RubinMe(means = t(SAVE_PARAM), vars = t(SAVE_VAR), impNum =IMPNUM)
	
	##############
	### Return ###
	##############
	
	if(ASSUME == 'ProportionalHazard'){
		TransCov$Trans14 = c(TransCov$Trans14, 'INT')
	}
	if(BASELINE == 'weib')
	{
		VARNAMES = c(paste(TransCov$Trans13, '_13', sep = ''),paste(TransCov$Trans24, '_24', sep = ''),
					paste(TransCov$Trans14, '_14', sep = ''), paste(TransCov$Trans34, '_34', sep = ''),
					paste(c('Intercept',TransCov$PNonCure), '_NonCure', sep = ''),
					paste('Scale', c('_13', '_24', '_14', '_34'), sep = ''), 
					paste('Shape', c('_13', '_24', '_14', '_34'), sep = ''))		
	}else if(BASELINE == 'cox'){
		VARNAMES = c(paste(TransCov$Trans13, '_13', sep = ''),paste(TransCov$Trans24, '_24', sep = ''),
					paste(TransCov$Trans14, '_14', sep = ''), paste(TransCov$Trans34, '_34', sep = ''),
					paste(c('Intercept',TransCov$PNonCure), '_NonCure', sep = ''))		
	}

	return(data.frame(VARNAMES,Estimate = OUT[,1], Variance = OUT[,2], v = OUT[,3]))
}
		


