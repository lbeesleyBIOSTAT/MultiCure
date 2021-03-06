

usethis::use_package("survival")
usethis::use_package("MASS")
usethis::use_package("reshape2")
usethis::use_package("cubature")
usethis::use_package("boot")
usethis::use_package("glmnet")
usethis::use_package("SurvRegCensCov")
usethis::use_package("shiny")
usethis::use_package("mice")




#' MultiCure
#' @description This function fits a Multistate Cure model using Expectation-Maximization (EM) and Monte Carlo Expectation-Maximization (MCEM) algorithms as in Beesley et al. (2018) in Biostatistics.
#'
#' @param iternum number of iterations for the EM or MCEM algorithm
#' @param datWIDE A data frame with the following columns (names must match): 
#' \itemize{
#' \item Y_R, the recurrence event/censoring time
#' \item delta_R, the recurrence event/censoring indicator
#'\item Y_D, the death event/censoring time 
#' \item delta_D, the death event/censoring indicator
#' \item G, the cure status variable. This takes value 1 for known non-cured, 0 for "known" cured and NA for unknown cur`e status
#'}
#' @param Cov A data frame containing the covariates used in the model fit. The columns must be named. Factors must be represented as dummy variables. If ridge or lasso penalties are being used, the covariates should be rescaled to have unit variances. 
#' @param ASSUME This variables indicates what equality assumptions we are making regarding the 24 and 14 transitions. The possible options are:
#' \itemize{
#' \item 'SameHazard': Lambda_14(t) = Lambda_24(t)
#' \item 'AllSeparate': No restrictions on Lambda_14(t) and Lambda_24(t)
#' \item  'ProportionalHazard': Lambda_14(t) = Lambda_24(t) exp(beta0)
#' \item  'SameBaseHaz': Lambda^0_14(t) = Lambda^0_24(t), No restrictions on beta_14 and beta_24
#' }
#' @param TransCov a list with elements: Trans13, Trans24, Trans14, Trans34, PNonCure. Each list element is a vector containing the names of the variables in Cov to be used in the model for the corresponding transition. 13 is NonCured -> Recurrence, 24 is Cured -> Death, 14 is NonCured -> Death, 34 is Recurrence -> Death. PNonCure contains the names of the covariates for the logistic regression for P(NonCure). 
#' @param IMPNUM number of imputed datasets. This is only used when covariates and/or outcome values are being imputed. 
#' @param BASELINE This variable indicates the assumptions about the baseline hazard form. This can take values 'weib' and 'cox'
#' @param PENALTY This variable indicates whether we are using any variable selection in the model fitting. The current code has been implemented and tested for option 'None' (no variable selection). Additional options include 'Ridge' (ridge regression for all covariates in all models) and 'Lasso' (lasso for all covariates in all models, only implemented for Cox baseline hazards), but these two additional options have not been rigorously tested.
#' @param PARAMINIT If desired, this can be a vector with initializations for the model parameters. The ordering of these parameters is c(beta, alpha, scale, shape) using the same ordering as in the output
#' @param UNEQUALCENSIMPUTE This is a function for imputing the outcome data in the unequal censoring (follow-up) setting. This only needs to be specified when we have unequal censoring. Several default options are included in this package, but this could also be a user-specified function. Inputs and outputs must match default versions.
#' @param COVIMPUTEFUNCTION This is a function for creating a single imputed version of the covariate set when covariate imputation is needed. This is user-specified. See COVIMPUTEFUNCTION_Example.R for an example of the input and output structure. 
#' @param COVIMPUTEINITIALIZE This is a function for initializing the missing values of the covariates. This is user-specified. See COVIMPUTEINITIALIZE_Example.R for an example of the input and output structure. 
#'
#' @return fit If BASELINE = 'weib', this is a list containing
#' \itemize{
#' \item beta estimate at final iteration. The ordering is: Beta for Transition 1->3, Beta for Transition 2->4, Beta for Transition 1->4, Beta_0 if BASELINE equals 'ProportionalHazard', Beta for Transition 3->4
#' \item alpha estimate at final iteration.
#' \item scale estimate at final iteration. The ordering is: Transition 1->3, 2->4, 1->4, 3->4 
#' \item shape estimate at final iteration. The ordering is: Transition 1->3, 2->4, 1->4, 3->4
#' \item fit will also contain estimates of beta, alpha, scale, and shape from each iteration and, if imputation is performed, the imputed values of Cov, G, Y_R, and delta_R from the last iteration. If imputation is not performed, includes the most recent weight p.
#'}
#' If BASELINE = 'cox', this is a list containing 
#' \itemize{
#' \item beta estimate at final iteration. The ordering is: Beta for Transition 1->3, Beta for Transition 2->4, Beta for Transition 1->4, Beta_0 if BASELINE equals 'ProportionalHazard', Beta for Transition 3->4
#' \item alpha estimate at final iterations. 
#' \item fit will also contain estimates of beta and alpha from each iteration and, if imputation is performed, the imputed values of Cov, G, Y_R, and delta_R from the last iteration. If imputation is not performed, includes the most recent weight p.
#'}
#' @details In order to fit a model with no covariates for one or more of the transitions or the logistic regression, include an all-zero covariate in Cov and list that covariate for the corresponding transition/s in TransCov. 
#'
#' In order to include recurrence time in the model for recurrence -> death, include covariate 'T_R' (initialized to equal the observed recurrence event/censoring time) in both Cov and Trans34 (in TransCov). If performing imputation for unequal follow-up, user must specify COVIMPUTEINITIALIZE and COVIMPUTEFUNCTION to update the values of 'T_R' based on the imputed outcome values.
#' @examples
#' attach(SimulateMultiCure(type = "NoMissingness"))
#' Cov = data.frame(X1,X2)
#' VARS = names(Cov)
#' TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#' datWIDE = data.frame( Y_R, Y_D, delta_R , delta_D, G)
#' fit = MultiCure(iternum = 100, datWIDE, Cov, ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "weib") 
#' @export



MultiCure = function(iternum, datWIDE, Cov, COVIMPUTEFUNCTION = NULL,  COVIMPUTEINITIALIZE = NULL, UNEQUALCENSIMPUTE = NULL, ASSUME = 'SameHazard', TransCov, IMPNUM = NULL, BASELINE = 'weib', PENALTY = 'None', PARAMINIT = NULL){
	
	###########################
	### Checking Arguments ####
	###########################
	Nobs = length(datWIDE[,1])
	if(is.null(TransCov) | is.null(Cov)){stop('TransCov and Cov must be specified. To fit a model without covariates, include all-zero covariate in the model for each transition')}
	ASSUME = match.arg(ASSUME, choices = c('SameHazard', 'AllSeparate', 'SameBaseHaz', 'ProportionalHazard'))
	BASELINE = match.arg(BASELINE, choices = c('weib','cox'))
	PENALTY = match.arg(PENALTY, choices = c('None', 'Ridge', 'Lasso'))
	UnequalCens = ifelse(datWIDE$Y_R < datWIDE$Y_D & datWIDE$delta_R == 0 & is.na(datWIDE$G), 1, 0)
	CovMissing = apply(Cov,1:2,is.na)
	NEEDTOIMPUTE = (sum(CovMissing)!=0 | sum(UnequalCens) != 0)	
	if(ASSUME %in% c('SameHazard', 'ProportionalHazard') & sum(TransCov$Trans24!= TransCov$Trans14) != 0 ){
		print('Transition 24 and 14 Covariates Do Not Match. Setting Trans14 Equal to Trans24 Covariates')
		TransCov$Trans24= TransCov$Trans14
	}
	if(ASSUME == 'ProportionalHazard'){
		Cov$INT = rep(1,length(Cov[,1]))
		TransCov$Trans14 = c(TransCov$Trans14, 'INT')
	}
	if(PENALTY == 'Lasso' & BASELINE == 'weib'){stop('Lasso penalization not yet implemented for Weibull baselines')}	
	if(NEEDTOIMPUTE & is.null(IMPNUM)){stop('Specify Number of Imputations')}
	if(sum(UnequalCens) != 0 & is.null(UNEQUALCENSIMPUTE)){
		if(BASELINE == 'weib' & 'T_R' %in% TransCov$Trans34){
			UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBINVERSION
		}else if(BASELINE == 'weib' & !('T_R' %in% TransCov$Trans34)){
			UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBREJECTION
		}else{UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTECOXMH}	
	}
	if((sum(CovMissing) != 0 | 'T_R' %in% TransCov$Trans34) & (is.null(COVIMPUTEFUNCTION) | is.null(COVIMPUTEINITIALIZE))){stop('Must Specify Covariate Initialization and Imputation Functions. Note: When "T_R" is in TransCov$Trans34 and we have unequal follow-up, these functions are used to update the values of T_R.')	}	
	
	
	####################
	### Initialize p ### (EM Algorithm)
	####################
	if(NEEDTOIMPUTE & sum(CovMissing)!=0){	
		CovTEMP =  COVIMPUTEINITIALIZE(Cov, CovMissing)
	}else{CovTEMP = Cov}	
	fitTemp = stats::glm(datWIDE$delta_R~.,data = CovTEMP, family = 'binomial') 
	predictions = as.numeric(stats::predict(fitTemp,type = 'response'))
	datWIDE$p = rep(NA,length(datWIDE[,1]))
	datWIDE$p = ifelse(is.na(datWIDE$G), predictions, datWIDE$G)	
		
	##########################
	### Initialize Imputes ###
	##########################
	
	if(NEEDTOIMPUTE ){
		GImp = replicate(IMPNUM,datWIDE$G)
		YRImp = replicate(IMPNUM,datWIDE$Y_R)
		deltaRImp = replicate(IMPNUM,datWIDE$delta_R)
		CovImp = replicate(IMPNUM,list(Cov))
		YRImpSAVE = replicate(IMPNUM,datWIDE$Y_R)	
		TAU_R = max(datWIDE$Y_R[datWIDE$delta_R==1]) #latest obseved recurrence											
		for(i in 1:IMPNUM){
			
			### Initialize G ###
			GImp[is.na(GImp[,i]) & datWIDE$Y_R > TAU_R ,i] = 0 #Set subjects at risk after TAU_R to be non-cured
			Draws = sapply(rep(mean(datWIDE$delta_R),Nobs),mSample)
			GImp[is.na(GImp[,i]),i] = Draws[is.na(GImp[,i])] #Draw remaining G status using logistic glm above
			
			### Initialize Unequal Censoring Imputations ###
			if(sum(datWIDE$UnequalCens)!=0){		
				### If GImp == 0, No Recurrence
					YRImp[UnequalCens == 1 & GImp[,i]==0,i] = datWIDE$Y_D[UnequalCens == 1 & GImp[,i]==0]
					deltaRImp[UnequalCens == 1 & GImp[,i]==0,i] = 0 #G=0, no recurrence				
				### If GImp == 1 and Yd is after last recurrence event, Recurrence
					deltaRImp[UnequalCens == 1 & GImp[,i]==1 & datWIDE$Y_D >= TAU_R,i] = 1 #G=1 and Yd >= last recurrence, deltaR = 1
				### If GImp == 1 and Yd is before last recurrence event, Draw recurrence supposing Recurrence is U(0, TAU_R). P(DeltaR = 1) is expressed as:
					prob = ((datWIDE$Y_D -datWIDE$Y_R)/(TAU_R-datWIDE$Y_R))[UnequalCens == 1 & GImp[,i]==1 & datWIDE$Y_D < TAU_R]
					Draws = sapply(prob,mSample)
					deltaRImp[UnequalCens == 1 & GImp[,i]==1 & datWIDE$Y_D < TAU_R,i] = Draws				
				### Initialize YrImpSAVE (only used in Metropolis-Hastings Cox imputation with unequal follow-up)
					IMPUTEYR = (UnequalCens ==1)
					MIN = datWIDE$Y_R[IMPUTEYR]
					MAX = pmin(datWIDE$Y_D[IMPUTEYR],TAU_R)
					MIN = ifelse(MIN >= MAX, MAX, MIN)
					U = apply(cbind(MIN, MAX),1, mHPropose) 
					YRImpSAVE[IMPUTEYR,i] = U				
				### Initialize Remaining YrImp
					YRImp[UnequalCens==1 & GImp[,i]==1 & deltaRImp[,i]==1,i] = 	YRImpSAVE[UnequalCens==1 & GImp[,i]==1 & deltaRImp[,i]==1,i]		
					YRImp[UnequalCens==1 & GImp[,i]==1 & deltaRImp[,i]==0,i] = 	datWIDE$Y_D[UnequalCens==1 & GImp[,i]==1 & deltaRImp[,i]==0]	
			}
			
			### Initialize Missing Covariates ###
			if(sum(CovMissing) != 0){
				CovImp[[i]] = COVIMPUTEINITIALIZE(Cov, CovMissing)
			}
		}		
		ImputeDat = list(UnequalCens, CovMissing, CovImp, GImp, YRImp, deltaRImp, YRImpSAVE )
	}else if(!NEEDTOIMPUTE){
		ImputeDat = c()
	}
	
	
	###########################################
	### Multistate Cure with No Missingness ### (This includes no missingness in G. This is used to fit the multistate cure model to complate data)
	###########################################

	if(!NEEDTOIMPUTE & sum(is.na(datWIDE$G))==0){ #NO ITERATION REQUIRED
		#print('No Iteration Needed')
		NEEDTOIMPUTE = TRUE #Uses data in ImputeDat for model fit
		ImputeDat = list(UnequalCens = NULL, CovMissing = NULL, CovImp = list(Cov= Cov), GImp = matrix(datWIDE$G), 
			YRImp = matrix(datWIDE$Y_R), deltaRImp = matrix(datWIDE$delta_R))
		if(BASELINE == 'weib'){
			param = MStep_WEIB(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)
			beta = param[[1]]
			alpha = param[[2]]
			scale = param[[3]]
			shape = param[[4]]	
			return(list(beta, alpha, scale, shape))	
		}else{
			param = MStep_COX(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)
			beta = param[[1]]
			alpha = param[[2]]
			return(list(beta, alpha))	
		}
	}

	
	#############################
	### Initialize Parameters ###
	#############################

	beta_save = c()
	alpha_save = c()
	scale_save = c()
	shape_save = c()
	p_save = c()
	l_save = c()
	iter = 1	

	if(BASELINE == 'weib'){
		if(!is.null(PARAMINIT)){
			beta = PARAMINIT$beta
			alpha = PARAMINIT$alpha
			scale = PARAMINIT$scale
			shape = PARAMINIT$shape				
		}else{
			param = MStep_WEIB(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)
			beta = param[[1]]
			alpha = param[[2]]
			scale = param[[3]]
			shape = param[[4]]	
		}
		beta_save = cbind(beta_save, beta)
		alpha_save = cbind(alpha_save, alpha)
		scale_save = cbind(scale_save, scale)
		shape_save = cbind(shape_save, shape)	
	}else{
		if(!is.null(PARAMINIT)){
			beta = PARAMINIT$beta
			alpha = PARAMINIT$alpha		
		}else{
			param = MStep_COX(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)			
			beta = param[[1]]
			alpha = param[[2]]	
		}
		beta_save = cbind(beta_save, beta)
		alpha_save = cbind(alpha_save, alpha)

	}
	
	######################		
	### EM for Weibull ###
	######################		
		
	if(BASELINE == 'weib'){
		while(iter <= iternum){ #ITERATE BETWEEN E-STEP AND M-STEP
			print(iter)
			### E Step
			if(NEEDTOIMPUTE){
				imputes = EStepWEIB_MC(datWIDE, beta, alpha, scale, shape, ImputeDat, COVIMPUTEFUNCTION, 
						UNEQUALCENSIMPUTE, TransCov)
				CovImp = imputes[[1]]
				GImp = imputes[[2]]
				YRImp = imputes[[3]]
				deltaRImp = imputes[[4]]
				ImputeDat[[3]] = CovImp
				ImputeDat[[4]] = GImp
				ImputeDat[[5]] = YRImp
				ImputeDat[[6]] = deltaRImp				
			}else{
				datWIDE = EStepWEIB(datWIDE, Cov, beta, alpha, scale, shape, TransCov)
			}
			
			### M Step
			param = MStep_WEIB(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)
			beta = param[[1]]
			alpha = param[[2]]
			scale = param[[3]]
			shape = param[[4]]
			beta_save = cbind(beta_save, beta)
			alpha_save = cbind(alpha_save, alpha)
			scale_save = cbind(scale_save, scale)
			shape_save = cbind(shape_save, shape)

			iter = iter + 1
			gc() #NEW
		}#end while loop	
	}
	
	##################		
	### EM for Cox ###
	##################	
	
	
	if(BASELINE == 'cox'){
		while(iter <= iternum){ #ITERATE BETWEEN E-STEP AND M-STEP
			print(iter)
			### E Step
			if(NEEDTOIMPUTE){
				imputes = EStepCOX_MC(datWIDE, beta, alpha, ImputeDat, COVIMPUTEFUNCTION, 
						UNEQUALCENSIMPUTE, TransCov, ASSUME)
				CovImp = imputes[[1]]
				GImp = imputes[[2]]
				YRImp = imputes[[3]]
				deltaRImp = imputes[[4]]
				YRImpSAVE = imputes[[9]]
				ImputeDat[[3]] = CovImp
				ImputeDat[[4]] = GImp
				ImputeDat[[5]] = YRImp
				ImputeDat[[6]] = deltaRImp		
				ImputeDat[[7]] = YRImpSAVE	
			}else{
				datWIDE = EStepCOX(datWIDE, Cov, beta, alpha, TransCov, ASSUME)
				p_save = cbind(p_save, datWIDE$p)
			}
			### M Step
			param = MStep_COX(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)		
			beta = param[[1]]
			alpha = param[[2]]
			beta_save = cbind(beta_save, beta)
			alpha_save = cbind(alpha_save, alpha)
			iter = iter + 1
			gc() #NEW
		}#end while loop
	}


	
	##############		
	### Return ###
	##############
	if(BASELINE == 'weib'){
		if(NEEDTOIMPUTE){
			return(list(beta = beta, alpha= alpha, scale = scale, shape = shape, 
						beta_save = beta_save, alpha_save = alpha_save, scale_save = scale_save, shape_save = shape_save, 
						CovImp = CovImp, GImp = GImp, YRImp = YRImp, deltaRImp = deltaRImp))
		}else if(!NEEDTOIMPUTE){	
			return(list(beta = beta, alpha= alpha, scale = scale, shape = shape, 
						beta_save = beta_save, alpha_save = alpha_save, scale_save = scale_save, shape_save = shape_save))
		}
	}else{
		if(NEEDTOIMPUTE){
			return(list(beta = beta, alpha= alpha, beta_save = beta_save, alpha_save = alpha_save, 
						CovImp = CovImp, GImp = GImp, YRImp = YRImp, deltaRImp = deltaRImp))
		}else if(!NEEDTOIMPUTE){	
			return(list(beta = beta, alpha= alpha, beta_save = beta_save, alpha_save = alpha_save, p_save = p_save))
		}			
	}#end ifelse
}#end function


