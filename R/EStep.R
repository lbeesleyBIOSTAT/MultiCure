 
### These functions perform either the E step or Imputation Step of the EM or MonteCarlo EM (respectively) algorithms


########################################################
### E-Step for Monte Carlo EM with Weibull Baselines ### (Imputation Step)
########################################################
#' @export

EStepWEIB_MC = function(datWIDE, beta, alpha, scale, shape, ImputeDat, COVIMPUTEFUNCTION, 
					UNEQUALCENSIMPUTE, TransCov){
	
	UnequalCens = ImputeDat[[1]]
	CovMissing = ImputeDat[[2]]
	CovImp = ImputeDat[[3]]
	GImp = ImputeDat[[4]]
	YRImp = ImputeDat[[5]]
	deltaRImp = ImputeDat[[6]]
	Nobs = length(datWIDE[,1])
	IMPNUM = length(CovImp)	
	if('T_R' %in% TransCov$Trans34){
		fd<-function(v){
			XB_beta34MOD = as.numeric(beta[TRANS==4][TransCov$Trans34!= 'T_R'] %*% t(cbind(CovImp[[i]][m,TransCov$Trans34[TransCov$Trans34!='T_R']])))
			XB_beta34MOD = XB_beta34MOD + as.numeric(beta[TRANS==4][TransCov$Trans34== 'T_R'] %*% t(cbind(v)))	
			S1 = exp(-(scale[1]*((v)^shape[1]) ) *exp(XB_beta13[m]))*exp(-(scale[3]*((v)^shape[3]) )*exp(XB_beta14[m]))
			S3 = exp(-(scale[4]*((datWIDE$Y_D[m]-v)^shape[4]) )*exp(XB_beta34MOD))		
			h13 = (scale[1]*shape[1]*((v)^(shape[1]-1))  )*exp(XB_beta13[m])
			h34 = ifelse(v==datWIDE$Y_D[m],0,(scale[4]*shape[4]*((datWIDE$Y_D[m]-v)^(shape[4]-1))  )*exp(XB_beta34MOD)	)
			return(h13*S1*S3*((h34)^datWIDE$delta_D[m]))   
		}  				
	}else{
		fd<-function(v){
			S1 = exp(-(scale[1]*((v)^shape[1]) ) *exp(XB_beta13[m]))*exp(-(scale[3]*((v)^shape[3]) )*exp(XB_beta14[m]))
			S3 = exp(-(scale[4]*((datWIDE$Y_D[m]-v)^shape[4]) )*exp(XB_beta34[m]))		
			h13 = (scale[1]*shape[1]*((v)^(shape[1]-1))  )*exp(XB_beta13[m])
			h34 = ifelse(v==datWIDE$Y_D[m],0,(scale[4]*shape[4]*((datWIDE$Y_D[m]-v)^(shape[4]-1))  )*exp(XB_beta34[m])	)
			return(h13*S1*S3*((h34)^datWIDE$delta_D[m]))   
		}  		
	}

	###################
	### Calculate Y ###
	###################
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))	
	y = matrix(0,ncol = IMPNUM, nrow = Nobs)		
	if(sum(UnequalCens)!=0){ #have unequal censoring		
		### Initialize XBs
		i = 1
		XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[[i]][,TransCov$Trans13])))	
		XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[[i]][,TransCov$Trans14])))	
		XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[[i]][,TransCov$Trans34])))				
		options(warn = -1)		
		if(sum(CovMissing) == 0){	 #No covariate missingness
			for(m in 1:Nobs){
				if(UnequalCens[m]==1){
					g=cubature::adaptIntegrate(Vectorize(fd), lowerLimit = datWIDE$Y_R[m], upperLimit = datWIDE$Y_D[m],maxEval=15)  
					y[m,1:IMPNUM]=g$integral
				} 
			}
		}else{	 ###Covariate missingness for some subjects			
			for(m in 1:Nobs){	
				g=cubature::adaptIntegrate(Vectorize(fd), lowerLimit = datWIDE$Y_R[m], upperLimit = datWIDE$Y_D[m],maxEval=15)  
				y[m,1:IMPNUM]=g$integral				
				if(sum(CovMissing[m,]) == 1 & UnequalCens[m]==1 & IMPNUM > 1){
					for(i in 2:IMPNUM){
						XB_beta13[m] = as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[[i]][m,TransCov$Trans13])))	
						XB_beta14[m] = as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[[i]][m,TransCov$Trans14])))	
						XB_beta34[m] = as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[[i]][m,TransCov$Trans34])))		
						g= cubature::adaptIntegrate(Vectorize(fd), lowerLimit = datWIDE$Y_R[m], upperLimit = datWIDE$Y_D[m],maxEval=15)  
						y[m,i]=g$integral
					}#end for loop
				}#end if	
			}#end for loop					
		}#end ifelse		
		options(warn = 1)	
	}#end if
	
	
	################
	### Impute G ###
	################
			
	for(i in 1:IMPNUM){
		XB_alpha = as.numeric(alpha[1:(1+length(TransCov$PNonCure))] %*% t(cbind(rep(1,Nobs), CovImp[[i]][,TransCov$PNonCure])))
		prob_Noncure = exp(XB_alpha)/(1+exp(XB_alpha))	
		XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[[i]][,TransCov$Trans13])))	
		XB_beta24 = as.numeric(beta[TRANS==2] %*% t(cbind(CovImp[[i]][,TransCov$Trans24])))	
		XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[[i]][,TransCov$Trans14])))	
		XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[[i]][,TransCov$Trans34])))		
		S1_D = exp(- (scale[1]*((datWIDE$Y_D)^shape[1]) ) *exp(XB_beta13))*
			exp(-(scale[3]*((datWIDE$Y_D)^shape[3]) )*exp(XB_beta14))
		S2_D = exp(-(scale[2]*((datWIDE$Y_D)^shape[2]) )*exp(XB_beta24))
		h24_D = (scale[2]*shape[2]*((datWIDE$Y_D)^(shape[2]-1))      )*exp(XB_beta24)
		h14_D = (scale[3]*shape[3]*((datWIDE$Y_D)^(shape[3]-1))  )*exp(XB_beta14)
		A = (prob_Noncure*S1_D*ifelse(datWIDE$delta_D == 1, h14_D, rep(1,Nobs) ) ) + (prob_Noncure*UnequalCens*y[,i])
		B = (1-prob_Noncure)*S2_D*ifelse(datWIDE$delta_D == 1, h24_D, rep(1,Nobs) )
		C = ifelse(A==0, rep(0,Nobs), A/(A+B))
		Draws = sapply(C,mSample)
		GImp[is.na(datWIDE$G),i] = Draws[is.na(datWIDE$G)]		
	}
	
	###############################
	### Impute Unequal Outcomes ###
	###############################	
	
	if(sum(UnequalCens) != 0 ){
		for(i in 1:IMPNUM){
			ImputeDatSHORT = list(UnequalCens, CovMissing, CovImp[[i]], GImp[,i], YRImp[,i], deltaRImp[,i], y=y[,i])
			ImputedOutcomes = UNEQUALCENSIMPUTE(datWIDE, beta, alpha, scale, shape, ImputeDatSHORT, TransCov)
			deltaRImp[,i] = ImputedOutcomes[[1]]
			YRImp[,i] = ImputedOutcomes[[2]]		
		}		
	}
		
	#########################
	### Impute Covariates ###
	#########################
	
	if(sum(CovMissing) != 0 | 'T_R' %in% TransCov$Trans34){
		for(i in 1:IMPNUM){
			ImputeDatSHORT = list(UnequalCens, CovMissing, CovImp[[i]], GImp[,i], YRImp[,i], deltaRImp[,i])
			CovImp[[i]] = COVIMPUTEFUNCTION(datWIDE, c(beta, alpha, scale, shape), ImputeDatSHORT, TransCov)	
		}		
	}		
	
	return(list(CovImp, GImp, YRImp, deltaRImp))
}





####################################################
### E-Step for Monte Carlo EM with Cox Baselines ### (Imputation Step)
####################################################


 #' @export


EStepCOX_MC = function(datWIDE, beta, alpha, ImputeDat, COVIMPUTEFUNCTION, 
					UNEQUALCENSIMPUTE, TransCov, ASSUME){
	UnequalCens = ImputeDat[[1]]
	CovMissing = ImputeDat[[2]]
	CovImp = ImputeDat[[3]]
	GImp = ImputeDat[[4]]
	YRImp = ImputeDat[[5]]
	deltaRImp = ImputeDat[[6]]
	YRImpSAVE = ImputeDat[[7]]
	Nobs = length(datWIDE[,1])	
 	
	if('T_R' %in% TransCov$Trans34){
		fdCOX<-function(v){	
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
		fdCOX<-function(v){	
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
	IMPNUM = length(CovImp)
	
	#################################
	### Estimate Baseline Hazards ###
	#################################
	
	Ests = Baselinehazard_IMP(datWIDE, CovImp,GImp, YRImp,deltaRImp, beta, alpha, TransCov, ASSUME)
	Basehaz13 = Ests[[1]]
	Basehaz24 = Ests[[2]]
	Basehaz14 = Ests[[3]]
	Basehaz34 = Ests[[4]]
	BasehazFun_13 = stepfun(x= Basehaz13[,2], y = c(Basehaz13[,3],0), right = F)
	BasehazFun_24 = stepfun(x= Basehaz24[,2], y = c(Basehaz24[,3],0), right = F)
	BasehazFun_14 = stepfun(x= Basehaz14[,2], y = c(Basehaz14[,3],0), right = F)
	BasehazFun_34 = stepfun(x= Basehaz34[,2], y = c(Basehaz34[,3],0), right = F)



	###################
	### Calculate Y ###
	###################
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))	
	y = matrix(0,ncol = IMPNUM, nrow = Nobs)		
	if(sum(UnequalCens)!=0){ #have unequal censoring		
		### Initialize XBs
		i = 1
		XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[[i]][,TransCov$Trans13])))	
		XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[[i]][,TransCov$Trans14])))	
		XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[[i]][,TransCov$Trans34])))				
		options(warn = -1)		
		if(sum(CovMissing) == 0){	 #No covariate missingness
			for(m in 1:Nobs){
				if(UnequalCens[m]==1){
					g=cubature::adaptIntegrate(Vectorize(fdCOX), lowerLimit = datWIDE$Y_R[m], upperLimit = datWIDE$Y_D[m],maxEval=15)  
					y[m,1:IMPNUM]=g$integral
				} 
			}
		}else{	 ###Covariate missingness for some subjects			
			for(m in 1:Nobs){	
				g=cubature::adaptIntegrate(Vectorize(fdCOX), lowerLimit = datWIDE$Y_R[m], upperLimit = datWIDE$Y_D[m],maxEval=15)  
				y[m,1:IMPNUM]=g$integral				
				if(sum(CovMissing[m,]) == 1 & UnequalCens[m]==1){
					for(i in 2:IMPNUM){
						XB_beta13[m] = as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[[i]][m,TransCov$Trans13])))	
						XB_beta14[m] = as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[[i]][m,TransCov$Trans14])))	
						XB_beta34[m] = as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[[i]][m,TransCov$Trans34])))		
						g=cubature::adaptIntegrate(Vectorize(fdCOX), lowerLimit = datWIDE$Y_R[m], upperLimit = datWIDE$Y_D[m],maxEval=15)  
						y[m,i]=g$integral
					}#end for loop
				}#end if	
			}#end for loop					
		}#end ifelse		
		options(warn = 1)	
	}#end if
	
	
	################
	### Impute G ###
	################
	GImpSAVE = GImp
	A13 = sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz13)
	A14 = sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz14)
	A24 = sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz24)
	B24 = BasehazFun_24(datWIDE$Y_D)
	B14 = BasehazFun_14(datWIDE$Y_D)

	for(i in 1:IMPNUM){
		XB_alpha = as.numeric(alpha[1:(1+length(TransCov$PNonCure))] %*% t(cbind(rep(1,Nobs), CovImp[[i]][,TransCov$PNonCure])))
		prob_Noncure = exp(XB_alpha)/(1+exp(XB_alpha))	
		XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[[i]][,TransCov$Trans13])))	
		XB_beta24 = as.numeric(beta[TRANS==2] %*% t(cbind(CovImp[[i]][,TransCov$Trans24])))	
		XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[[i]][,TransCov$Trans14])))	
		XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[[i]][,TransCov$Trans34])))		
		S1_D = exp(-as.numeric(A13)*exp(XB_beta13))*exp(-as.numeric(A14)*exp(XB_beta14))
		S2_D = exp(-as.numeric(A24)*exp(XB_beta24))
		h24_D = B24*exp(XB_beta24)
		h14_D = B14*exp(XB_beta14)
		A = (prob_Noncure*S1_D*ifelse(datWIDE$delta_D == 1, h14_D, rep(1,Nobs) ) ) + (prob_Noncure*UnequalCens*y[,i])
		B = (1-prob_Noncure)*S2_D*ifelse(datWIDE$delta_D == 1, h24_D, rep(1,Nobs) )
		C = ifelse(A==0, rep(0,Nobs), A/(A+B))
		Draws = sapply(C,mSample)
		GImp[is.na(datWIDE$G),i] = Draws[is.na(datWIDE$G)]	
	}
	
	###############################
	### Impute Unequal Outcomes ###
	###############################	
	
	if(sum(UnequalCens) != 0 ){
		for(i in 1:IMPNUM){
			ImputeDatSHORT = list(UnequalCens, CovMissing, CovImp[[i]], GImp[,i], YRImp[,i], deltaRImp[,i], y=y[,i], 
				Basehaz13, Basehaz24, Basehaz14, Basehaz34, YRImpSAVE[,i], GImpSAVE[,i])
			ImputedOutcomes = UNEQUALCENSIMPUTE(datWIDE, beta, alpha, ImputeDat = ImputeDatSHORT, TransCov)
			deltaRImp[,i] = ImputedOutcomes[[1]]
			YRImp[,i] = ImputedOutcomes[[2]]
			if(length(ImputedOutcomes)==3){
				YRImpSAVE[,i] = ImputedOutcomes[[3]]
			}		
		}	
	}
	
	#########################
	### Impute Covariates ###
	#########################
	
	if(sum(CovMissing) != 0 | 'T_R' %in% TransCov$Trans34){
		for(i in 1:IMPNUM){
			ImputeDatSHORT = list(UnequalCens, CovMissing, CovImp[[i]], GImp[,i], YRImp[,i], deltaRImp[,i], Basehaz13, Basehaz24, Basehaz14, Basehaz34)
			CovImp[[i]] = COVIMPUTEFUNCTION(datWIDE, c(beta, alpha), ImputeDat= ImputeDatSHORT, TransCov)	
		}		
	}		
	
	return(list(CovImp, GImp, YRImp, deltaRImp, Basehaz13, Basehaz24, Basehaz14, Basehaz34, YRImpSAVE))
}



#' @export



############################################
### E-Step for EM with Weibull Baselines ###
############################################
 

EStepWEIB = function(datWIDE, Cov, beta, alpha, scale, shape, TransCov){	
	Nobs = length(datWIDE[,1])	
	XB_alpha = as.numeric(alpha[1:(1+length(TransCov$PNonCure))] %*% t(cbind(rep(1,Nobs), Cov[,TransCov$PNonCure])))
	prob_Noncure = exp(XB_alpha)/(1+exp(XB_alpha))	
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))	
	XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(Cov[,TransCov$Trans13])))	
	XB_beta24 = as.numeric(beta[TRANS==2] %*% t(cbind(Cov[,TransCov$Trans24])))	
	XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(Cov[,TransCov$Trans14])))	
	XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(Cov[,TransCov$Trans34])))							
	p = rep(NA,Nobs)	
	S1_D = exp(- (scale[1]*((datWIDE$Y_D)^shape[1]) ) *exp(XB_beta13))*exp(-(scale[3]*((datWIDE$Y_D)^shape[3]) )*exp(XB_beta14))
	S2_D = exp(-(scale[2]*((datWIDE$Y_D)^shape[2]) )*exp(XB_beta24))
	h24_D = (scale[2]*shape[2]*((datWIDE$Y_D)^(shape[2]-1))      )*exp(XB_beta24)
	h14_D = (scale[3]*shape[3]*((datWIDE$Y_D)^(shape[3]-1))  )*exp(XB_beta14)
	A = prob_Noncure*S1_D*ifelse(datWIDE$delta_D == 1, h14_D, rep(1,Nobs) )
	B = (1-prob_Noncure)*S2_D*ifelse(datWIDE$delta_D == 1, h24_D, rep(1,Nobs) )
	C = ifelse(A==0, rep(0,Nobs),A/(A+B))
	datWIDE$p = ifelse(is.na(datWIDE$G), C, datWIDE$G)
	return(datWIDE)
}



########################################
### E-Step for EM with Cox Baselines ###
########################################


#' @export

 

EStepCOX = function(datWIDE, Cov, beta, alpha, TransCov, ASSUME){	
	Nobs = length(datWIDE[,1])	
	XB_alpha = as.numeric(alpha[1:(1+length(TransCov$PNonCure))] %*% t(cbind(rep(1,Nobs), Cov[,TransCov$PNonCure])))
	prob_Noncure = exp(XB_alpha)/(1+exp(XB_alpha))	
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))
	XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(Cov[,TransCov$Trans13])))	
	XB_beta24 = as.numeric(beta[TRANS==2] %*% t(cbind(Cov[,TransCov$Trans24])))	
	XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(Cov[,TransCov$Trans14])))	
	XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(Cov[,TransCov$Trans34])))							
	
	#################################
	### Estimate Baseline Hazards ###
	#################################
	Ests = BaselineHazard_NOIMP(datWIDE, Cov, beta, alpha, TransCov, ASSUME, p = datWIDE$p)
	Haz_13 = Ests[[1]]
	Haz_24 = Ests[[2]]
	Haz_14 = Ests[[3]]
	Haz_34 = Ests[[4]]
	
	p = rep(NA,Nobs)	
	dOtherCauses = ifelse(datWIDE$delta_D==1 & datWIDE$delta_R == 0, 1, 0)
	S1_D = exp(-Haz_13(datWIDE$Y_D)*exp(XB_beta13))*exp(-Haz_14(datWIDE$Y_D)*exp(XB_beta14))
	S2_D = exp(-Haz_24(datWIDE$Y_D)*exp(XB_beta24))
	h24_D = ifelse(dOtherCauses==1,(Haz_24(datWIDE$Y_D+0.0001)-  Haz_24(datWIDE$Y_D-0.0001)    )*exp(XB_beta24), rep(0,Nobs))
	h14_D = ifelse(dOtherCauses==1,(Haz_14(datWIDE$Y_D+0.0001)-  Haz_14(datWIDE$Y_D-0.0001)    )*exp(XB_beta14), rep(0,Nobs))
	A = prob_Noncure*S1_D*ifelse(datWIDE$delta_D == 1, h14_D, rep(1,Nobs) )
	B = (1-prob_Noncure)*S2_D*ifelse(datWIDE$delta_D == 1, h24_D, rep(1,Nobs) )
	C = ifelse(A==0, rep(0,Nobs),A/(A+B))
	datWIDE$p = ifelse(is.na(datWIDE$G), C, datWIDE$G)
	return(datWIDE)
}




