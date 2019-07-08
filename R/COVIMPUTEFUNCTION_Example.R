


#' COVIMPUTEFUNCTION_Example
#' @description The function COVIMPUTEFUNCTION_Example creates a single imputed version of the covariate matrix. An example function is included in the package, but the user must specify their own function when applying MultiCure to datasets with covariate missingness. This function must have input/output as described below.

#' @param datWIDE defined as in MultiCure
#' @param param If baselines are WEIBULL, this is a vector containing beta, alpha, scale, and shape. If baselines are COX, this is a vector containing beta and alpha. 
#' @param ImputeDat This is a list with the following elements:
#' \itemize{
#' \item UnequalCens: A vector taking value 1 if the subject has unequal follow-up. Note: If subject is assumed cured in datWIDE, they are listed as UnequalCens = 0.
#' \item CovMissing: A matrix indicating which elements of Cov are missing. Not needed for this imputation.
#' \item CovImp: A list containing a single imputation of Cov
#' \item GImp: A vector with a recent single imputation of G
#' \item YRImp: A vector with a recent single imputation of Y_R
#' \item deltaRImp: A vector with a recent single imputation of delta_R
#' If baselines are COX, then this will also include
#' \item Basehaz13: A matrix containing the estimate of the baseline hazard function for the 1->3 transition specified intervals
#' \item Basehaz24: A matrix containing the estimate of the baseline hazard function for the 2->4 transition specified intervals
#' \item Basehaz14: A matrix containing the estimate of the baseline hazard function for the 1->4 transition specified intervals
#' \item Basehaz34: A matrix containing the estimate of the baseline hazard function for the 3->4 transition specified intervals
#' }
#' @param TransCov defined as in MultiCure
#'
#' @return CovImp a matrix with a SINGLE imputation of the covariate matrix
#' @details The example code included in the package imputes missing covariate X2 in the Multistate cure model example. In the example code, a normal covariate is imputed using an approach similar to SMC-FCS in Bartlett et al. (2014) and Metropolis-Hastings methods. In practice, this function can use any imputation method the user desires. For example, the user-written function can call 'mice' in R to perform the imputation. 
#' @export



COVIMPUTEFUNCTION_Example = function(datWIDE,param, ImputeDat, TransCov){
	BASELINE = ifelse(length(param)!=(length(reshape2::melt(TransCov)$value)+1), 'weib', 'cox')
	CovMissing = ImputeDat[[2]]
	CovImp = as.data.frame(ImputeDat[[3]])
	GImp = ImputeDat[[4]]
	YRImp = ImputeDat[[5]]
	deltaRImp = ImputeDat[[6]]	
	Nobs = length(datWIDE[,1])
	
	#################
	#################
	### Impute X2 ###
	#################
	#################

	#######################################
	### Estimate Parameters for X2 | X1 ###
	#######################################
	
	#Note: We do not draw the parameter values. We are doing improper imputation..
	fit = lm(X2~X1, data = CovImp)	
	theta = coef(fit)
 	sigma2 = (summary(fit)$sigma)^2


	mPropose = function(Y){
		return(rnorm(n=1, mean = Y, sd = 0.5))
	}
	current = CovImp[CovMissing[,'X2']==T,'X2']
	
	##################################
	### Propose New Imputed Values ###
	##################################
	proposal = sapply(current, mPropose)	
	
	XB = theta[1] + theta[2]*CovImp[CovMissing[,'X2']==T,'X1']
	dens_marg_CUR = dnorm(current, mean = XB, sd = sqrt(sigma2))
	dens_marg_PRO = dnorm(proposal, mean =XB, sd = sqrt(sigma2))

	TEMP = data.frame(datWIDE, CovImp, GImp, YRImp, deltaRImp)
	TEMP = TEMP[CovMissing[,'X2']==T,]	
	
	#######################################################################
	### Estimate Log-Likelihood for current and proposed imputed values ###
	#######################################################################
	
	if(BASELINE == 'weib'){
		TEMP[,c('X2')] = current
		LogLik_CUR = as.numeric(apply(TEMP,1,FUN = LogLikWEIB, TransCov = TransCov, param = param))
		TEMP[,c('X2')] = proposal
		LogLik_PRO = as.numeric(apply(TEMP,1,FUN = LogLikWEIB, TransCov = TransCov, param = param))				
	}else{
		TEMP[,c('X2')] = current
		LogLik_CUR = as.numeric(apply(TEMP,1,FUN = LogLikCOX, TransCov = TransCov, param = param, ImputeDat = ImputeDat))
		TEMP[,c('X2')] = proposal
		LogLik_PRO = as.numeric(apply(TEMP,1,FUN = LogLikCOX, TransCov = TransCov, param = param, ImputeDat = ImputeDat))				
	}

	alph<-runif(sum(CovMissing[,'X2']==T),0,1)

	####################################################################
	### Metropolis-Hastings, accept or reject proposed imputed value ###
	####################################################################
	
	ACCEPT = log(alph)<(LogLik_PRO +log(dens_marg_PRO)+log(dnorm(current, proposal,0.5))-LogLik_CUR-log(dens_marg_CUR)-log(dnorm(proposal, current,0.5)))
	CovImp[CovMissing[,'X2']==T,'X2'][ACCEPT] = proposal[ACCEPT]
	
	return(CovImp)
}





### These log-likelihood functions are used in the above implementation of covariate imputation, but they can be ignored for user-specified imputation functions


#' @export

LogLikWEIB = function(DAT, TransCov, param){	
	Y_D = as.numeric(DAT[c('Y_D')])
	delta_D = as.numeric(DAT[c('delta_D')])
	if('INT' %in% names(DAT)){
		CovImp = DAT[c('X1', 'X2', 'INT')]
	}else{
		CovImp = DAT[c('X1', 'X2')]		
	}
	GImp = as.numeric(DAT['GImp'])
	deltaRImp = as.numeric(DAT['deltaRImp'])
	YRImp = as.numeric(DAT['YRImp'])
	A = length(c(TransCov$Trans13, TransCov$Trans24, TransCov$Trans14, TransCov$Trans34))
	BreakParam = c(rep(1,A), rep(2,1+length(TransCov$PNonCure)), rep(3,4), rep(4,4))
	beta = param[BreakParam==1]
	alpha = param[BreakParam==2]
	scale = param[BreakParam==3]
	shape = param[BreakParam==4]
	XB_alpha = as.numeric(alpha %*% c(1, unlist(CovImp[TransCov$PNonCure])))
	prob_Noncure = exp(XB_alpha)/(1+exp(XB_alpha))	
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))
	XB_beta13 = as.numeric(beta[TRANS==1] %*% unlist(c(CovImp[TransCov$Trans13])))	
	XB_beta24 = as.numeric(beta[TRANS==2] %*% unlist(c(CovImp[TransCov$Trans24])))		
	XB_beta14 = as.numeric(beta[TRANS==3] %*% unlist(c(CovImp[TransCov$Trans14])))			
	XB_beta34 = as.numeric(beta[TRANS==4] %*% unlist(c(CovImp[TransCov$Trans34])))	
	S1_D = exp(- (scale[1]*((Y_D)^shape[1]) ) *exp(XB_beta13))*
		exp(-(scale[3]*((Y_D)^shape[3]) )*exp(XB_beta14))
	S1_R = exp(- (scale[1]*((YRImp)^shape[1]) ) *exp(XB_beta13))*exp(-(scale[3]*((YRImp)^shape[3]) )*exp(XB_beta14))
	S2_D = exp(-(scale[2]*((Y_D)^shape[2]) )*exp(XB_beta24))
	S3 = exp(-(scale[4]*((Y_D - YRImp)^shape[4]) )*exp(XB_beta34))
	h24_D = (scale[2]*shape[2]*((Y_D)^(shape[2]-1))      )*exp(XB_beta24)
	h14_D = (scale[3]*shape[3]*((Y_D)^(shape[3]-1))  )*exp(XB_beta14)
	h13_R = (scale[1]*shape[1]*((YRImp)^(shape[1]-1))  )*exp(XB_beta13)
	h34_D = (scale[4]*shape[4]*((Y_D-YRImp)^(shape[4]-1))  )*exp(XB_beta34)	
	L = (    (prob_Noncure*h13_R*S1_R*S3*(h34_D^delta_D))^as.numeric(GImp==1 & deltaRImp==1)  )*
	(    (prob_Noncure*(h14_D^delta_D)*S1_D)^as.numeric(GImp==1 & deltaRImp==0))*
	(	((1-prob_Noncure)*(h24_D^delta_D)*S2_D   )^as.numeric(GImp==0))	
	return(log(L))
}

#' @export

LogLikCOX = function(DAT, TransCov, param, ImputeDat){	
	Basehaz13 = ImputeDat[[7]]	
	Basehaz24 = ImputeDat[[8]]	
	Basehaz14 = ImputeDat[[9]]	
	Basehaz34 = ImputeDat[[10]]	
	BasehazFun_13 = stepfun(x= Basehaz13[,2], y = c(Basehaz13[,3],0), right = F)
	BasehazFun_24 = stepfun(x= Basehaz24[,2], y = c(Basehaz24[,3],0), right = F)
	BasehazFun_14 = stepfun(x= Basehaz14[,2], y = c(Basehaz14[,3],0), right = F)
	BasehazFun_34 = stepfun(x= Basehaz34[,2], y = c(Basehaz34[,3],0), right = F)
	Y_D = as.numeric(DAT[c('Y_D')])
	delta_D = as.numeric(DAT[c('delta_D')])
	if('INT' %in% names(DAT)){
		CovImp = DAT[c('X1', 'X2', 'INT')]
	}else{
		CovImp = DAT[c('X1', 'X2')]		
	}	
	GImp = as.numeric(DAT['GImp'])
	deltaRImp = as.numeric(DAT['deltaRImp'])
	YRImp = as.numeric(DAT['YRImp'])
	A = length(c(TransCov$Trans13, TransCov$Trans24, TransCov$Trans14, TransCov$Trans34))
	BreakParam = c(rep(1,A), rep(2,1+length(TransCov$PNonCure)))
	beta = param[BreakParam==1]
	alpha = param[BreakParam==2]
	XB_alpha = as.numeric(alpha %*% c(1, unlist(CovImp[TransCov$PNonCure])))
	prob_Noncure = exp(XB_alpha)/(1+exp(XB_alpha))	
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))
	XB_beta13 = as.numeric(beta[TRANS==1] %*% unlist(c(CovImp[TransCov$Trans13])))	
	XB_beta24 = as.numeric(beta[TRANS==2] %*% unlist(c(CovImp[TransCov$Trans24])))		
	XB_beta14 = as.numeric(beta[TRANS==3] %*% unlist(c(CovImp[TransCov$Trans14])))			
	XB_beta34 = as.numeric(beta[TRANS==4] %*% unlist(c(CovImp[TransCov$Trans34])))	
	S1_D = exp(-as.numeric(sapply(Y_D,Baseline_Hazard, Basehaz13))*exp(XB_beta13))*
		exp(-as.numeric(sapply(Y_D,Baseline_Hazard, Basehaz14))*exp(XB_beta14))
	S1_R = exp(-as.numeric(sapply(YRImp,Baseline_Hazard, Basehaz13))*exp(XB_beta13))*
		exp(-as.numeric(sapply(YRImp,Baseline_Hazard, Basehaz14))*exp(XB_beta14))	
	S2_D = exp(-as.numeric(sapply(Y_D,Baseline_Hazard, Basehaz24))*exp(XB_beta24))
	S3 = exp(-as.numeric(sapply(Y_D - YRImp,Baseline_Hazard, Basehaz34))*exp(XB_beta34))
	h24_D = BasehazFun_24(Y_D)*exp(XB_beta24)
	h14_D = BasehazFun_14(Y_D)*exp(XB_beta14)
	h13_R = BasehazFun_13(YRImp)*exp(XB_beta13)
	h34_D = BasehazFun_34(Y_D-YRImp)*exp(XB_beta34)	
	h34_D = ifelse(h34_D == 0, 0.001, h34_D)
	h13_R = ifelse(h13_R == 0, 0.001, h13_R)
	h14_D = ifelse(h14_D == 0, 0.001, h14_D)
	h24_D = ifelse(h24_D == 0, 0.001, h24_D)
	L = (    (prob_Noncure*h13_R*S1_R*S3*(h34_D^delta_D)                )^as.numeric(GImp==1 & deltaRImp==1)*
	(    (prob_Noncure*(h14_D^delta_D)*S1_D)                )^as.numeric(GImp==1 & deltaRImp==0))*
	(	((1-prob_Noncure)*(h24_D^delta_D)*S2_D   )^as.numeric(GImp==0))	
	return(log(L))
}


