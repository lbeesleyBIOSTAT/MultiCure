
#' ProperDraws_MC
#' @description The function ProperDraws_MC is used to obtain "proper" imputations of the missing data after the MCEM algorithm is used to fit the multistate cure model. These proper imputations are then used in the functions MultiCure_VAREST_Imputation or MultiCure_VAREST_ImputationBOOT to estimate the parameter standard errors. 
#'
#' @param datWIDE A data frame with the following columns: 
#' \itemize{
#' \item Y_R, the recurrence event/censoring time
#' \item delta_R, the recurrence event/censoring indicator
#'\item Y_D, the death event/censoring time 
#' \item delta_D, the death event/censoring indicator
#' \item G, the cure status variable. This takes value 1 for known non-cured, 0 for "known" cured and NA for unknown cure status
#'}
#' @param Cov matrix of covariates used in MultiCure (may have missingness)
#' @param CovImp  A list with IMPNUM elements containing the imputations of Cov output from MultiCure
#' @param GImp  A matrix with IMPNUM elements containing the imputations of G output from MultiCure
#' @param YRImp  A matrix with IMPNUM elements containing the imputations of Y_R output from MultiCure
#' @param deltaRImp  A matrix with IMPNUM elements containing the imputations of delta_R output from MultiCure
#' @param COVIMPUTEFUNCTION This is a function for creating a single imputed version of the covariate set when covariate imputation is needed. This is user-specified. See XXXXXX for an example of the input and output structure. 
#' @param COVIMPUTEINITIALIZE This is a function for initializing the missing values of the covariates. This is user-specified. See XXXXXX for an example of the input and output structure. 
#' @param UNEQUALCENSIMPUTE This is a function for imputing the outcome data in the unequal censoring (follow-up) setting. This only needs to be specified when we have unequal censoring. Several default options exist, but this could also be a user-specified function. Inputs and outputs must match default versions.
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
#' \item CovImp  A list with IMPNUM elements containing "proper" imputations of Cov
#' \item GImp  A list with IMPNUM elements containing "proper" imputations of G
#' \item YRImp  A list with IMPNUM elements containing "proper" imputations of Y_R
#' \item deltaRImp  A list with IMPNUM elements containing "proper" imputations of delta_R
#'}
#' @details In order to output the imputed data from MultiCure, one must use the trace = TRUE option in MultiCure.
#'
#' @author Lauren J Beesley, \email{lbeesley@umich.edu}
#' @export

ProperDraws_MC = function( datWIDE,Cov,CovImp, GImp, YRImp, deltaRImp, COVIMPUTEFUNCTION = NULL,  COVIMPUTEINITIALIZE = NULL,
			UNEQUALCENSIMPUTE = NULL, ASSUME = 'SameHazard', TransCov, BASELINE, PENALTY = 'None',POSTITER = 5){
	Nobs = length(datWIDE[,1])
	ASSUME = match.arg(ASSUME, choices = c('SameHazard', 'AllSeparate', 'SameBaseHaz',  'ProportionalHazard'))
	BASELINE = match.arg(BASELINE, choices = c('weib','cox'))
	PENALTY = match.arg(PENALTY, choices = c('None', 'Ridge', 'Lasso'))
	UnequalCens = ifelse(datWIDE$Y_R < datWIDE$Y_D & datWIDE$delta_R == 0 & is.na(datWIDE$G), 1, 0)
	CovMissing = apply(Cov,1:2,is.na)
	NEEDTOIMPUTE = TRUE

	IMPNUM = length(CovImp)
	if(sum(UnequalCens) != 0 & is.null(UNEQUALCENSIMPUTE)){
		if(BASELINE == 'weib'){
			UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIB
		}else{UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTECOXMH}
	}
	if(sum(CovMissing) != 0 & (is.null(COVIMPUTEFUNCTION) | is.null(COVIMPUTEINITIALIZE))){stop('Must Specify Covariate Initialization and Imputation Functions')	}		
	if(ASSUME == 'ProportionalHazard'){
		Cov$INT = rep(1,length(Cov[,1]))
		TransCov$Trans14 = c(TransCov$Trans14, 'INT')
	}	
	
	FunTEMP = function(x){
			return(x[whichboot,])
	}
		
	#################################
	### Obtain Proper Imputations ###	
	#################################
	YRImpSAVE = YRImp
	for(i in 1:IMPNUM){
		TAU_R = max(datWIDE$Y_R[datWIDE$delta_R==1])
 		TAU_D = max(datWIDE$Y_D[datWIDE$delta_D==1])		
 		DIFFMAX =  max((datWIDE$Y_D -datWIDE$Y_R)[datWIDE$delta_D==1 & datWIDE$delta_R==1]) #max distance of Yr and Yd such that both events occur
		IMPUTEYR = (UnequalCens ==1)
		MIN = pmax(datWIDE$Y_R[IMPUTEYR],datWIDE$Y_D[IMPUTEYR] -DIFFMAX)
		MAX = pmin(datWIDE$Y_D[IMPUTEYR], TAU_R)
		MIN = ifelse(MIN >= MAX, MAX, MIN)
		U = runif(sum(IMPUTEYR),min = MIN, max = MAX)
		YRImpSAVE[IMPUTEYR,i] = ifelse(deltaRImp[IMPUTEYR,i]==1, YRImp[IMPUTEYR,i], U)	
	}
	iter = 1
	while(iter <= POSTITER)
	{
		for(i in 1:IMPNUM){
			whichboot = sample(x=c(1:Nobs), size = Nobs, replace = TRUE, prob = rep(1/Nobs, Nobs))
			ImputeDatBOOT = list(UnequalCens=UnequalCens[whichboot], CovMissing=CovMissing[whichboot,], 
							CovImp= list( CovImp[[i]][whichboot,]), GImp= matrix(GImp[whichboot,i]), 
							YRImp= matrix(YRImp[whichboot,i]), deltaRImp= matrix(deltaRImp[whichboot,i]), YRImpSAVE = matrix(YRImpSAVE[whichboot,i])  )	
			ImputeDat = list(UnequalCens= UnequalCens, CovMissing= CovMissing, CovImp  = list(CovImp[[i]]), GImp  = matrix(GImp[,i]), YRImp  = matrix(YRImp[,i]), 
							deltaRImp  = matrix(deltaRImp[,i]), YRImpSAVE = matrix(YRImpSAVE[,i]))		
			if(BASELINE == 'weib'){		
				param = MStep_WEIB(datWIDE = datWIDE[whichboot,], Cov = Cov[whichboot,], ImputeDat = ImputeDatBOOT, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)
				beta = param[[1]]
				alpha = param[[2]]
				scale = param[[3]]
				shape = param[[4]]
				imputes = EStepWEIB_MC(datWIDE, beta, alpha, scale, shape, ImputeDat, COVIMPUTEFUNCTION, 
						UNEQUALCENSIMPUTE, TransCov)
				CovImp[[i]] = data.frame(imputes[[1]])
				GImp[,i] = imputes[[2]]
				YRImp[,i] = imputes[[3]]
				deltaRImp[,i] = imputes[[4]]		
			}else{		
				param = MStep_COX(datWIDE = datWIDE[whichboot,], Cov=Cov[whichboot,], ImputeDat = ImputeDatBOOT, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY)
				beta = param[[1]]
				alpha = param[[2]]
				imputes = EStepCOX_MC(datWIDE, beta, alpha, ImputeDat, COVIMPUTEFUNCTION, 
						UNEQUALCENSIMPUTE, TransCov, ASSUME)
				CovImp[[i]] = data.frame(imputes[[1]])
				GImp[,i] = imputes[[2]]
				YRImp[,i] = imputes[[3]]
				deltaRImp[,i] = imputes[[4]]
				YRImpSAVE[,i] = imputes[[9]]
			}			
		}
		iter = iter+1
	}
	return(list(CovImp, GImp, YRImp, deltaRImp))
}#end function





