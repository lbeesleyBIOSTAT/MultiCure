
#' UNEQUALCENSIMPUTECOXNESTEDWEIB
#' @description The function UNEQUALCENSIMPUTECOXNESTEDWEIB will perform an imputation algorithm to handle unequal follow-up for recurrence and death. This function can be applied when we assume COX baseline hazards. Even though the main multistate cure model assumes Cox baseline hazards, this function (temporarily) uses Weibull baseline hazard assumptions to perform imputation of Y_R and delta_R. This function makes use of a rejection sampling algorithm.

#' @param datWIDE defined as in MultiCure
#' @param beta A vector containing the most recent estimates of beta
#' @param alpha A vector containing the most recent estimates of alpha
#' @param ImputeDat This is a list with the following elements:
#' \itemize{
#' \item UnequalCens: A vector taking value 1 if the subject has unequal follow-up. Note: If subject is assumed cured in datWIDE, they are listed as UnequalCens = 0.
#' \item CovMissing: A matrix indicating which elements of Cov are missing. Not needed for this imputation.
#' \item CovImp: A list containing a single imputation of Cov
#' \item GImp: A vector with a recent single imputation of G
#' \item YRImp: A vector with a recent single imputation of Y_R
#' \item deltaRImp: A vector with a recent single imputation of delta_R
#' \item y: The integral of the target kernel over Yr0 to Yd
#' \item Basehaz13: A matrix containing the estimate of the baseline hazard function for the 1->3 transition specified intervals (not used)
#' \item Basehaz24: A matrix containing the estimate of the baseline hazard function for the 2->4 transition specified intervals (not used)
#' \item Basehaz14: A matrix containing the estimate of the baseline hazard function for the 1->4 transition specified intervals (not used)
#' \item Basehaz34: A matrix containing the estimate of the baseline hazard function for the 3->4 transition specified intervals (not used)
#' }
#' @param TransCov defined as in MultiCure
#'
#' @return a list containing 
#' \itemize{
#' \item [[1]]: deltaRImp, A single imputation of delta_R
#' \item [[2]]: YRImp, A single imputation of Y_R
#'}
#' @export





UNEQUALCENSIMPUTECOXNESTEDWEIB = function(datWIDE, beta, alpha, ImputeDat, TransCov){
	
	##################
	### Initialize ###
	##################
	
	UnequalCens = ImputeDat[[1]]
	CovImp = as.data.frame(ImputeDat[[3]])
	GImp = ImputeDat[[4]]
	YRImp = ImputeDat[[5]]
	deltaRImp = ImputeDat[[6]]
	y = ImputeDat[[7]]
	Basehaz13 = ImputeDat[[8]]	
	Basehaz24 = ImputeDat[[9]]	
	Basehaz14 = ImputeDat[[10]]	
	Basehaz34 = ImputeDat[[11]]	
	GImpSAVE = ImputeDat[[13]]
	Nobs = length(datWIDE[,1])
	TAU_R = max(Basehaz13[,1])
	
	###################################################
	### Fit model assuming Weibull baseline hazards ###
	###################################################
	
	ImputeDatSHORT = list(UnequalCens = matrix(UnequalCens), CovMissing = matrix(CovMissing), CovImp = list(CovImp), GImp = matrix(GImpSAVE), 
						YRImp = matrix(YRImp), deltaRImp = matrix(deltaRImp))					
	param = MStep_WEIB(datWIDE= datWIDE, Cov = CovImp, ImputeDat  = ImputeDatSHORT, ASSUME, TransCov, NEEDTOIMPUTE = T, PENALTY = 'None')
	beta = param[[1]]
	alpha = param[[2]]
	scale = param[[3]]
	shape = param[[4]]

	###########################
	### Impute Delta_R, Y_R ###
	###########################
	
	### Use original datWIDE to update YRImp and deltaRImp
	if('T_R' %in% TransCov$Trans34){
		ImputedOutcomes = UNEQUALCENSIMPUTEWEIBINVERSION(datWIDE, beta, alpha, scale, shape, ImputeDat = ImputeDat, TransCov)
	}else{
		ImputedOutcomes = UNEQUALCENSIMPUTEWEIBREJECTION(datWIDE, beta, alpha, scale, shape, ImputeDat = ImputeDat, TransCov)
	}
	deltaRImp = ImputedOutcomes[[1]]
	YRImp = ImputedOutcomes[[2]]	
		
	return(list(deltaRImp, YRImp))
}

