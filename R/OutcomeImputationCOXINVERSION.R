
#' UNEQUALCENSIMPUTECOXINVERSION
#' @description The function UNEQUALCENSIMPUTECOXINVERSION will perform an imputation algorithm to handle unequal follow-up for recurrence and death. This function can be applied when we assume COX baseline hazards. This function performs imputation through inverting the survival function of the target distribution. 

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
#' }
#' @param TransCov defined as in MultiCure
#'
#' @return a list containing 
#' \itemize{
#' \item [[1]]: deltaRImp, A single imputation of delta_R
#' \item [[2]]: YRImp, A single imputation of Y_R
#'}
#' @export




UNEQUALCENSIMPUTECOXINVERSION = function(datWIDE, beta, alpha, ImputeDat, TransCov){
	
	
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
	
	TAU = max(datWIDE$Y_R[datWIDE$delta_R==1])
	
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
	
	
	########################
	### Define Functions ###
	########################
	
	if('T_R' %in% TransCov$Trans34){
		fdCOX<-function(v, m){	
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
		fdCOX<-function(v, m){	
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
	DrawVAL = function(TIME, U, m){		
			g=cubature::adaptIntegrate(Vectorize(fdCOX), lowerLimit = TIME, upperLimit = datWIDE$Y_D[m],m, maxEval=10)  
			ZERO=(g$integral/y[m])-U
		return(ZERO)
	}	

	##################
	### Impute T_R ### (By inverting the survival function of T_R)
	##################
	
	DrawVALWRAPPER = function(s){
		m = INDICES[s]	
		U1 = runif(n=1, min = 0, max = 1)
		draw = stats::uniroot(DrawVAL, interval = c(datWIDE$Y_R[m], datWIDE$Y_D[m]),U1, m, tol = 0.01, maxiter = 20)$root
		if(draw >= datWIDE$Y_D[m] ){draw = datWIDE$Y_D[m] - (datWIDE$Y_D[m]/1000)}
		if(draw <= datWIDE$Y_R[m] ){draw = datWIDE$Y_R[m] +  (datWIDE$Y_R[m]/1000)}
		#print(m)	
		return(draw)
	}		
	
	DRAWS = sapply(as.numeric(c(1:length(INDICES))), DrawVALWRAPPER)
	YRImp[INDICES] = DRAWS	
	return(list(deltaRImp, YRImp))
}


