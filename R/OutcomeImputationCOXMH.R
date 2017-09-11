
#' UNEQUALCENSIMPUTECOXMH
#' @description The function UNEQUALCENSIMPUTECOXMH will perform an imputation algorithm to handle unequal follow-up for recurrence and death. This function can be applied when we assume COX baseline hazards. This function performs imputation using a Metropolis-Hastings algorithm. The proposal distribution is Uniform with bounds such that the target kernel is nonzero. 

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
#' \item Basehaz13: A matrix containing the estimate of the baseline hazard function for the 1->3 transition specified intervals
#' \item Basehaz24: A matrix containing the estimate of the baseline hazard function for the 2->4 transition specified intervals
#' \item Basehaz14: A matrix containing the estimate of the baseline hazard function for the 1->4 transition specified intervals
#' \item Basehaz34: A matrix containing the estimate of the baseline hazard function for the 3->4 transition specified intervals
#' \item YRImpSAVE: A vecotr with the most recent ACCEPTED values of Y_R from the Metropolis-Hastings algorithm
#' }
#' @param TransCov defined as in MultiCure
#'
#' @return a list containing 
#' \itemize{
#' \item [[1]]: deltaRImp, A single imputation of delta_R
#' \item [[2]]: YRImp, A single imputation of Y_R
#'}
#' @author Lauren J Beesley, \email{lbeesley@umich.edu}
#' @export


UNEQUALCENSIMPUTECOXMH = function(datWIDE, beta, alpha, ImputeDat, TransCov){
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
	YRImpSAVE = ImputeDat[[12]]
	Nobs = length(datWIDE[,1])
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))
	XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[,TransCov$Trans13])))	
	XB_beta24 = as.numeric(beta[TRANS==2] %*% t(cbind(CovImp[,TransCov$Trans24])))	
	XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[,TransCov$Trans14])))	
	XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[,TransCov$Trans34])))

	BasehazFun_13 = stepfun(x= Basehaz13[,2], y = c(Basehaz13[,3],0), right = F)
	BasehazFun_24 = stepfun(x= Basehaz24[,2], y = c(Basehaz24[,3],0), right = F)
	BasehazFun_14 = stepfun(x= Basehaz14[,2], y = c(Basehaz14[,3],0), right = F)
	BasehazFun_34 = stepfun(x= Basehaz34[,2], y = c(Basehaz34[,3],0), right = F)

	S1_D = exp(-as.numeric(sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz13))*exp(XB_beta13))*
			exp(-as.numeric(sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz14))*exp(XB_beta14))
	h14_D = BasehazFun_14(datWIDE$Y_D)*exp(XB_beta14)
		
	YRImp = ifelse(GImp==0,datWIDE$Y_D, ifelse(GImp==1 & UnequalCens == 0,datWIDE$Y_R,rep(NA,Nobs) ))
	deltaRImp = ifelse(GImp==0,rep(0,Nobs), ifelse(GImp==1 & UnequalCens == 0,datWIDE$delta_R,rep(NA,Nobs) ))

	######################
	### Impute Delta R ###
	######################			
	num = y
	denom = (h14_D^datWIDE$delta_D)*S1_D
	ratio = ifelse(num==0,num,num/(num + denom))	[GImp==1 & UnequalCens == 1]			
	deltaRImp[GImp==1 & UnequalCens == 1] = apply(matrix(ratio), 1,mSample)   
	YRImp[GImp==1 & UnequalCens == 1 & deltaRImp==0] = datWIDE$Y_D[GImp==1 & UnequalCens == 1 & deltaRImp==0]

	INDICES = which(is.na(YRImp))
	
	if('T_R' %in% TransCov$Trans34){
		fdCOX<-function(x){	
			v = x[1]
			m = x[2]
			XB_beta34MOD = as.numeric(beta[TRANS==4][TransCov$Trans34!= 'T_R'] %*% t(cbind(CovImp[[i]][m,TransCov$Trans34[TransCov$Trans34!='T_R']])))
			XB_beta34MOD = XB_beta34MOD + as.numeric(beta[TRANS==4][TransCov$Trans34== 'T_R'] %*% t(cbind(v)))
			Cumhazard13_temp = exp(XB_beta13[m])*as.numeric(Baseline_Hazard(v, Basehaz13 ))
			Cumhazard14_temp = exp(XB_beta14[m])*as.numeric(Baseline_Hazard(v, Basehaz14 ))
			Cumhazard34_temp = exp(XB_beta34MOD)*as.numeric(Baseline_Hazard(datWIDE$Y_D[m]-v, Basehaz34) )		
			Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
			Surv3_temp = exp(-Cumhazard34_temp)		
			hazard13_temp = exp(XB_beta13[m])*BasehazFun_13(v) 
			hazard34_temp = ifelse(v == datWIDE$Y_D[m],0,exp(XB_beta34MOD)*BasehazFun_34(datWIDE$Y_D[m]-v))		
			return(hazard13_temp*Surv1_temp* Surv3_temp*((hazard34_temp)^datWIDE$delta_D[m]))   
		}  				
	}else{
		fdCOX<-function(x){	
			v = x[1]
			m = x[2]
			Cumhazard13_temp = exp(XB_beta13[m])*as.numeric(Baseline_Hazard(v, Basehaz13 ))
			Cumhazard14_temp = exp(XB_beta14[m])*as.numeric(Baseline_Hazard(v, Basehaz14 ))
			Cumhazard34_temp = exp(XB_beta34[m])*as.numeric(Baseline_Hazard(datWIDE$Y_D[m]-v, Basehaz34) )		
			Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
			Surv3_temp = exp(-Cumhazard34_temp)		
			hazard13_temp = exp(XB_beta13[m])*BasehazFun_13(v) 
			hazard34_temp = ifelse(v == datWIDE$Y_D[m],0,exp(XB_beta34[m])*BasehazFun_34(datWIDE$Y_D[m]-v))		
			return(hazard13_temp*Surv1_temp* Surv3_temp*((hazard34_temp)^datWIDE$delta_D[m]))   
		} 
	}

	TAU_R = max(Basehaz13[,1])

	current = YRImpSAVE[INDICES]	
	### Limits of proposal distribution determined so the baseline hazard and survival functions in fdCOX are nonzero. 
	MIN = datWIDE$Y_R[INDICES]
	MAX = pmin(datWIDE$Y_D[INDICES], TAU_R) 
	#For subjects in INDICES, MIN <= MAX. This is a result of setting lambda13(TAU_R) = 0, which assigns subjects at risk after TAU_R to G=0
	proposal = apply(cbind(MIN, MAX),1, mHPropose)	
	
	logdens_CUR = log(as.numeric(apply(cbind(current, INDICES),1, fdCOX)))
	logdens_PRO = log(as.numeric(apply(cbind(proposal, INDICES),1, fdCOX)))

	alph<-runif(length(INDICES),0,1)

	ACCEPT = log(alph)<(logdens_PRO +log(dunif(current, min = MIN, max = MAX))-logdens_CUR-log(dunif(proposal, min = MIN, max = MAX)))
	ACCEPT[is.na(ACCEPT)] = TRUE #should not ever be used. This is to catch errors in which fdCOX is infinite	
	
	YRImpSAVE[INDICES][!ACCEPT] = current[!ACCEPT]
	YRImpSAVE[INDICES][ACCEPT] = proposal[ACCEPT]
	
	YRImp[INDICES][ACCEPT] = proposal[ACCEPT]	
	YRImp[INDICES][!ACCEPT] = current[!ACCEPT]
	
	
	return(list(deltaRImp, YRImp, YRImpSAVE))
}


