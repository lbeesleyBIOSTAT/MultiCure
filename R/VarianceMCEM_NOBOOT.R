
#' VarianceMCEM_NOBOOT
#' @description The function VarianceMCEM_NOBOOT performs variance estimation for the Multistate Cure Model when the model is fit using a Monte Carlo EM Algorithm. For each imputed dataset, this function estimates parameter variances by fitting the multistate cure model the imputed dataset. Corresponding standard errors come are extracted from the corresponding logistic and proportional hazards model fits. Then, Rubin's rules are used to obtain a single set of parameter estimates and standard errors across imputed datasets. Note: This function has not yet been implemented for PENALTY values other than 'None'. Important!!!: The function ProperDraws_MC should be used to obtain the imputed datasets used in this function. This ensures that the imputed datasets are (roughly) proper imputations and that Rubin's rules can then be applied. 
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
#' \item v, the estimated degrees of freedom of the t-distribution of the parameter estimate from Rubin's Rules
#'}
#' @details This function provides parameter estimates and estimated variances. The parameter estimates are obtained using Rubin's rules, but an alternative estimate of the multistate cure model parameter can be obtained by averaging the parameter estimates from the last few iterations of the model fitting algorithm. In our experience, we found that the approach that averages across the last few iterations (rather than estimated using Rubin's rules) provides a better estimate of the parameter of interest.
#' @examples
#' attach(SimulateMultiCure(type = "UnequalCensoring"))
#' Cov = data.frame(X1,X2)
#' VARS = names(Cov)
#' TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#' datWIDE = data.frame( Y_R, Y_D, delta_R , delta_D, G)
#' fit = MultiCure(iternum = 100, datWIDE, Cov, ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "weib", IMPNUM = 10) ### Note: This will take a moment
#' Proper = ProperDraws_MC(datWIDE,Cov, CovImp = fit[[9]], GImp = fit[[10]], YRImp = fit[[11]], deltaRImp = fit[[12]], ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "weib") ### Note: This will take a moment
#' OUT = VarianceMCEM_NOBOOT(fit,datWIDE, CovImp = Proper[[1]], GImp = Proper[[2]], YRImp = Proper[[3]], deltaRImp = Proper[[4]],  ASSUME = "SameHazard", TransCov, BASELINE = "weib")
#' @export

VarianceMCEM_NOBOOT = function(fit,datWIDE, CovImp, GImp, YRImp, deltaRImp, ASSUME, TransCov, BASELINE, PENALTY = 'None'){			


	
	########################################
	### Initializing, Defining Functions ###
	########################################
	
	if(PENALTY != 'None'){stop('Not Implemented for Ridge or Lasso Penalties')}
	IMPNUM = length(CovImp)
	SAVE_VAR = c()
	SAVE_PARAM = c()

	ASSUME = match.arg(ASSUME, choices = c('SameHazard', 'AllSeparate', 'SameBaseHaz', 'ProportionalHazard'))
	BASELINE = match.arg(BASELINE, choices = c('weib','cox'))

	#####################################################
	### Fit Multistate Cure Model to Imputed Datasets ### (Estimate variances using coxph or survreg functions)
	#####################################################
	
	for(i in 1:IMPNUM){
		CovImpTEMP = CovImp[[i]]
		datWIDETEMP = datWIDE
		datWIDETEMP$G = GImp[,i]
		datWIDETEMP$delta_R = deltaRImp[,i]
		datWIDETEMP$Y_R = YRImp[,i]
		datWIDETEMP = subset(datWIDETEMP, select = c(Y_R, Y_D, delta_R, delta_D, G))

		ImputeDat = list(UnequalCens = NULL, CovMissing = NULL, CovImp = list(Cov= CovImpTEMP), GImp = matrix(datWIDETEMP$G), 
			YRImp =  matrix(datWIDETEMP$Y_R), deltaRImp =  matrix(datWIDETEMP$delta_R))
		if(BASELINE == 'weib')
		{
			param = MStep_WEIBVarEst(datWIDETEMP, CovImpTEMP, ImputeDat, ASSUME,  TransCov)
			beta = param[[1]]
			alpha = param[[2]]
			scale = param[[3]]
			shape = param[[4]]	
			VAR_beta = param[[5]]
			VAR_alpha = param[[6]]
			VAR_scale = param[[7]]
			VAR_shape = param[[8]]	
			SAVE_PARAM = rbind(SAVE_PARAM, c(beta, alpha, scale, shape))
			SAVE_VAR =rbind(SAVE_VAR, c(VAR_beta, VAR_alpha, VAR_scale, VAR_shape))
		}else{
			param = MStep_COXVarEst(datWIDETEMP, CovImpTEMP, ImputeDat, ASSUME,  TransCov)
			beta = param[[1]]
			alpha = param[[2]]
			VAR_beta = param[[3]]
			VAR_alpha = param[[4]]
			SAVE_PARAM = rbind(SAVE_PARAM, c(beta, alpha))
			SAVE_VAR =rbind(SAVE_VAR, c(VAR_beta, VAR_alpha))
		}
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
	if(BASELINE == 'weib' )
	{
		PARAM = OUT[,1]#c(fit[[1]], fit[[2]], fit[[3]], fit[[4]])	
		VARNAMES = c(paste(TransCov$Trans13, '_13', sep = ''),paste(TransCov$Trans24, '_24', sep = ''),
					paste(TransCov$Trans14, '_14', sep = ''), paste(TransCov$Trans34, '_34', sep = ''),
					paste(c('Intercept',TransCov$PNonCure), '_NonCure', sep = ''),
					paste('Scale', c('_13', '_24', '_14', '_34'), sep = ''), 
					paste('Shape', c('_13', '_24', '_14', '_34'), sep = ''))		
	}else if(BASELINE == 'cox' ){
		PARAM = OUT[,1]# c(fit[[1]], fit[[2]])	
		VARNAMES = c(paste(TransCov$Trans13, '_13', sep = ''),paste(TransCov$Trans24, '_24', sep = ''),
					paste(TransCov$Trans14, '_14', sep = ''), paste(TransCov$Trans34, '_34', sep = ''),
					paste(c('Intercept',TransCov$PNonCure), '_NonCure', sep = ''))		
	}
	

	return(data.frame(VARNAMES,Estimate = OUT[,1], Variance = OUT[,2], v = OUT[,3]))
}
		




### This function performs the variance estimation in a given imputed dataset under Weibull baseline hazards

#' @export

MStep_WEIBVarEst = function(datWIDE, Cov,ImputeDat, ASSUME, TransCov){
	###########################################
	### Transform the data into long format ###
	###########################################
	
	datLONG = CreateLong_MC(datWIDE, ImputeDat)
	datLONG_sub = datLONG[datLONG$w != 0,]		
	
	#######################################################
	### Transform the covariates into the proper format ###
	#######################################################
	
	Cov_long = subset(datLONG_sub, select = -c(id, from, to, trans, Tstart, Tstop, time, status,w))
	Cov_long_13 = Cov_long[,TransCov$Trans13]
	Cov_long_13[!(datLONG_sub$from == 1 & datLONG_sub$to == 3),] = 0
	Cov_long_34 = Cov_long[,TransCov$Trans34] 
	Cov_long_34[!(datLONG_sub$from == 3 & datLONG_sub$to == 4),] = 0
	if(ASSUME %in% c('SameHazard', 'SameCoeff', 'ProportionalHazard')){		
		Cov_long_1424 = Cov_long[,TransCov$Trans24] 
		Cov_long_1424[!((datLONG_sub$from == 2 | datLONG_sub$from == 1) & datLONG_sub$to == 4),] = 0
	}else if(ASSUME %in% c('AllSeparate', 'SameBaseHaz')){		
		Cov_long_14 = Cov_long[,TransCov$Trans14]
		Cov_long_14[!(datLONG_sub$from == 1 & datLONG_sub$to == 4),] = 0
		Cov_long_24 = Cov_long[,TransCov$Trans24] 
		Cov_long_24[!(datLONG_sub$from == 2 & datLONG_sub$to == 4),] = 0
	}		
	TRANS = I	
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)	
	
	###############################
	### Fit Weibull Regressions ###
	###############################
	
	fit13 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)), 
			weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 1))
	CONVERT13 = SurvRegCensCov::ConvertWeibull(fit13, conf.level = 0.95)
	mshape13 = CONVERT13[[1]][2,1]
	mscale13 = CONVERT13[[1]][1,1]
	coef13 = CONVERT13[[1]][3:length(CONVERT13[[1]][,1]),1]	
	
	fit34 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_34)),
		  	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 4))
	CONVERT34 = SurvRegCensCov::ConvertWeibull(fit34, conf.level = 0.95)
	mshape34 = CONVERT34[[1]][2,1]
	mscale34 = CONVERT34[[1]][1,1]
	coef34 = CONVERT34[[1]][3:length(CONVERT34[[1]][,1]),1]
			
	if(ASSUME %in% c('AllSeparate')){
		fit24 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_24)),
			 weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2))	
		CONVERT24 = SurvRegCensCov::ConvertWeibull(fit24, conf.level = 0.95)
		mshape24 = CONVERT24[[1]][2,1]
		mscale24 = CONVERT24[[1]][1,1]
		coef24 = CONVERT24[[1]][3:length(CONVERT24[[1]][,1]),1]	
		fit14 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_14)),
			 weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 3))
		CONVERT14 = SurvRegCensCov::ConvertWeibull(fit24, conf.level = 0.95)
		mshape14 = CONVERT14[[1]][2,1]
		mscale14 = CONVERT14[[1]][1,1]
		coef14 = CONVERT14[[1]][3:length(CONVERT14[[1]][,1]),1]
		beta = c(coef13, coef24, coef14, coef34)
		scale = as.numeric(c(mscale13, mscale24, mscale14, mscale34))
		shape = as.numeric(c(mshape13, mshape24, mshape14, mshape34))
		VAR_SHAPE =c( (CONVERT13[[1]][2,2])^2, (CONVERT24[[1]][2,2])^2, (CONVERT14[[1]][2,2])^2, (CONVERT34[[1]][2,2])^2 )
		VAR_SCALE =c( (CONVERT13[[1]][1,2])^2,  (CONVERT24[[1]][1,2])^2, (CONVERT14[[1]][1,2])^2, (CONVERT34[[1]][1,2])^2)
		VAR_BETA =c( (CONVERT13[[1]][3:length(CONVERT13[[1]][,1]),2])^2,  (CONVERT24[[1]][3:length(CONVERT24[[1]][,1]),2])^2,
					(CONVERT14[[1]][3:length(CONVERT14[[1]][,1]),2])^2, (CONVERT34[[1]][3:length(CONVERT34[[1]][,1]),2])^2)
			
	}else if(ASSUME %in% c('SameHazard')){
		fit2414 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_1424)),
		 	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2 | datLONG_sub$trans == 3))			
		CONVERT2414 = SurvRegCensCov::ConvertWeibull(fit2414, conf.level = 0.95)
		mshape2414 = CONVERT2414[[1]][2,1]
		mscale2414 = CONVERT2414[[1]][1,1]
		coef2414 = CONVERT2414[[1]][3:length(CONVERT2414[[1]][,1]),1]
		beta = c(coef13, coef2414, coef2414, coef34)
		scale = as.numeric(c(mscale13, mscale2414, mscale2414, mscale34)) 
		shape = as.numeric(c(mshape13, mshape2414, mshape2414, mshape34))
		
		VAR_SHAPE =c( (CONVERT13[[1]][2,2])^2, rep((CONVERT2414[[1]][2,2])^2,2),  (CONVERT34[[1]][2,2])^2 )
		VAR_SCALE =c( (CONVERT13[[1]][1,2])^2,  rep((CONVERT2414[[1]][1,2])^2,2), (CONVERT34[[1]][1,2])^2)
		VAR_BETA =c( (CONVERT13[[1]][3:length(CONVERT13[[1]][,1]),2])^2, rep( (CONVERT2414[[1]][3:length(CONVERT2414[[1]][,1]),2])^2,2),
					 (CONVERT34[[1]][3:length(CONVERT34[[1]][,1]),2])^2)
			
	}else if(ASSUME %in% c('SameBaseHaz')){
		fit2414 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_24)) + 
			TRANS(as.matrix(Cov_long_14)),
		 	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2 | datLONG_sub$trans == 3))				
		CONVERT2414 = SurvRegCensCov::ConvertWeibull(fit2414, conf.level = 0.95)
		mshape2414 = CONVERT2414[[1]][2,1]
		mscale2414 = CONVERT2414[[1]][1,1]
		coef2414 = CONVERT2414[[1]][3:length(CONVERT2414[[1]][,1]),1]
		A2 = length(TransCov$Trans24)
		A3 = length(TransCov$Trans14)
		TRANSITION = c(rep(2,A2), rep(3,A3))
		coef24 = coef2414[TRANSITION==2]
		coef14 = coef2414[TRANSITION==3]	
		beta = c(coef13, coef24, coef14, coef34)
		scale = as.numeric(c(mscale13, mscale2414, mscale2414, mscale34)) 
		shape = as.numeric(c(mshape13, mshape2414, mshape2414, mshape34))
		
		SEs = CONVERT2414[[1]][3:length(CONVERT2414[[1]][,1]),2]
		VAR_SHAPE =c( (CONVERT13[[1]][2,2])^2, rep((CONVERT2414[[1]][2,2])^2,2),  (CONVERT34[[1]][2,2])^2 )
		VAR_SCALE =c( (CONVERT13[[1]][1,2])^2,  rep((CONVERT2414[[1]][1,2])^2,2), (CONVERT34[[1]][1,2])^2)
		VAR_BETA =c( (CONVERT13[[1]][3:length(CONVERT13[[1]][,1]),2])^2,  
					(SEs[TRANSITION == 2])^2,(SEs[TRANSITION == 3])^2, 
					(CONVERT34[[1]][3:length(CONVERT34[[1]][,1]),2])^2)
			
		
	}else if(ASSUME %in% c('ProportionalHazard')){
		fit2414 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_1424)) + as.numeric(datLONG_sub$trans==3),
		 	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2 | datLONG_sub$trans == 3))	
		CONVERT2414 = SurvRegCensCov::ConvertWeibull(fit2414, conf.level = 0.95)
		mshape2414 = CONVERT2414[[1]][2,1]
		mscale2414 = CONVERT2414[[1]][1,1]
		coef2414 = CONVERT2414[[1]][3:length(CONVERT2414[[1]][,1]),1]		
		beta = c(coef13, coef2414[1:(length(coef2414)-1)], coef2414, coef34)
		scale = as.numeric(c(mscale13, mscale2414, mscale2414, mscale34)) 
		shape = as.numeric(c(mshape13, mshape2414, mshape2414, mshape34))
		
		SEs = CONVERT2414[[1]][3:length(CONVERT2414[[1]][,1]),2]
		VAR_SHAPE =c( (CONVERT13[[1]][2,2])^2, rep((CONVERT2414[[1]][2,2])^2,2),  (CONVERT34[[1]][2,2])^2 )
		VAR_SCALE =c( (CONVERT13[[1]][1,2])^2,  rep((CONVERT2414[[1]][1,2])^2,2), (CONVERT34[[1]][1,2])^2)
		VAR_BETA =c( (CONVERT13[[1]][3:length(CONVERT13[[1]][,1]),2])^2,  
					(SEs[1:(length(SEs)-1)])^2,(SEs)^2, 
					(CONVERT34[[1]][3:length(CONVERT34[[1]][,1]),2])^2)
	}
	
	###############################
	### Fit Logistic Regression ###
	###############################
	
	fitLogistic = stats::glm(as.numeric(ImputeDat[[4]])~as.matrix(Cov[,TransCov$PNonCure]), family = binomial(link = 'logit'))
	alpha = coef(fitLogistic)
	VARIANCEALPHA = diag(summary(fitLogistic)$cov.scaled)	
	
	##############
	### Return ###
	##############
	
	return(list(beta = as.numeric(beta),alpha = as.numeric(alpha), scale = as.numeric(scale), shape = as.numeric(shape),
		VAR_beta = as.numeric(VAR_BETA), VAR_alpha = as.numeric(VARIANCEALPHA) , 
		VAR_scale = as.numeric(VAR_SCALE) , VAR_shape = as.numeric(VAR_SHAPE) ))
}



### This function performs the variance estimation in a given imputed dataset under Cox baseline hazards
#' @export

MStep_COXVarEst = function(datWIDE, Cov, ImputeDat, ASSUME, TransCov){
	
	###########################################
	### Transform the data into long format ###
	###########################################
	
	datLONG = CreateLong_MC(datWIDE, ImputeDat)
	datLONG_sub = datLONG[datLONG$w != 0,]	
	
	#######################################################
	### Transform the covariates into the proper format ###
	#######################################################
	
	Cov_long = subset(datLONG_sub, select = -c(id, from, to, trans, Tstart, Tstop, time, status,w))
	Cov_long_13 = Cov_long[,TransCov$Trans13]
	Cov_long_13[!(datLONG_sub$from == 1 & datLONG_sub$to == 3),] = 0
	Cov_long_34 = Cov_long[,TransCov$Trans34] 
	Cov_long_34[!(datLONG_sub$from == 3 & datLONG_sub$to == 4),] = 0
	if(ASSUME %in% c( 'SameHazard', 'SameCoeff', 'ProportionalHazard')){		
		Cov_long_1424 = Cov_long[,TransCov$Trans24] 
		Cov_long_1424[!((datLONG_sub$from == 2 | datLONG_sub$from == 1) & datLONG_sub$to == 4),] = 0
	}else if(ASSUME %in% c('AllSeparate', 'SameBaseHaz')){		
		Cov_long_14 = Cov_long[,TransCov$Trans14]
		Cov_long_14[!(datLONG_sub$from == 1 & datLONG_sub$to == 4),] = 0
		Cov_long_24 = Cov_long[,TransCov$Trans24] 
		Cov_long_24[!(datLONG_sub$from == 2 & datLONG_sub$to == 4),] = 0
	}	
	TRANS = I
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	
	###########################
	### Fit Cox Regressions ###
	###########################
	
	strata = survival::strata
	if(ASSUME %in% c('AllSeparate')){
		fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
				TRANS(as.matrix(Cov_long_24)) + TRANS(as.matrix(Cov_long_14)) + TRANS(as.matrix(Cov_long_34)) +
				strata(datLONG_sub$trans), weights = datLONG_sub$w)
		beta = coef(fitCox)		
		VARIANCE = (summary(fitCox)$coefficients[,3])^2
	}else if(ASSUME %in% c('SameCoeff')){
		fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
				TRANS(as.matrix(Cov_long_1424)) + TRANS(as.matrix(Cov_long_34)) +
				strata(datLONG_sub$trans), weights = datLONG_sub$w)
		beta_short = coef(fitCox)
		VARIANCE_short = (summary(fitCox)$coefficients[,3])^2
		TRANSITION = c(rep(1,A1), rep(2,A2), rep(4,A4))
		beta = c(beta_short[TRANSITION==1], rep(beta_short[TRANSITION==2],2), beta_short[TRANSITION==4])	
		VARIANCE = 	c(VARIANCE_short[TRANSITION==1], rep(VARIANCE_short[TRANSITION==2],2), VARIANCE_short[TRANSITION==4])
	}else if(ASSUME %in% c('SameHazard')){
		Trans_MERGE = datLONG_sub$trans
		Trans_MERGE[Trans_MERGE==3] = 2
		fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
				TRANS(as.matrix(Cov_long_1424)) + TRANS(as.matrix(Cov_long_34)) +
				strata(Trans_MERGE), weights = datLONG_sub$w)
		beta_short = coef(fitCox)
		VARIANCE_short = (summary(fitCox)$coefficients[,3])^2
		TRANSITION = c(rep(1,A1), rep(2,A2), rep(4,A4))
		beta = c(beta_short[TRANSITION==1], rep(beta_short[TRANSITION==2],2), beta_short[TRANSITION==4])		
		VARIANCE = c(VARIANCE_short[TRANSITION==1], rep(VARIANCE_short[TRANSITION==2],2), VARIANCE_short[TRANSITION==4])		
	}else if(ASSUME %in% c('SameBaseHaz')){
		Trans_MERGE = datLONG_sub$trans
		Trans_MERGE[Trans_MERGE==3] = 2
		fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
				TRANS(as.matrix(Cov_long_24)) + TRANS(as.matrix(Cov_long_14)) + TRANS(as.matrix(Cov_long_34)) +
				strata(Trans_MERGE), weights = datLONG_sub$w)
		beta = coef(fitCox)
		VARIANCE = (summary(fitCox)$coefficients[,3])^2
	}else if(ASSUME %in% c('ProportionalHazard')){
		Trans_MERGE = datLONG_sub$trans
		Trans_MERGE[Trans_MERGE==3] = 2
		fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
				TRANS(as.matrix(Cov_long_1424)) + TRANS(as.matrix(Cov_long_34)) + as.numeric(datLONG_sub$trans==3)+
				strata(Trans_MERGE), weights = datLONG_sub$w)
		beta_short = coef(fitCox)
		VARIANCE_short = (summary(fitCox)$coefficients[,3])^2
		TRANSITION = c(rep(1,A1), rep(2,A2), rep(4,A4), 9)
		beta = c(beta_short[TRANSITION==1], rep(beta_short[TRANSITION==2],2), beta_short[TRANSITION==9], beta_short[TRANSITION==4])	
		VARIANCE = c(VARIANCE_short[TRANSITION==1], rep(VARIANCE_short[TRANSITION==2],2), VARIANCE_short[TRANSITION==9], VARIANCE_short[TRANSITION==4])
	}
	
	###############################
	### Fit Logistic Regression ###
	###############################
	
	fitLogistic = stats::glm(as.numeric(ImputeDat[[4]])~as.matrix(Cov[,TransCov$PNonCure]), family = binomial(link = 'logit'))
	alpha = coef(fitLogistic)
	VARIANCEALPHA = diag(summary(fitLogistic)$cov.scaled)	
		
	##############
	### Return ###
	##############
	
	return(list(beta = as.numeric(beta),alpha = as.numeric(alpha), VAR_beta = as.numeric(VARIANCE), VAR_alpha = as.numeric(VARIANCEALPHA)))
}


