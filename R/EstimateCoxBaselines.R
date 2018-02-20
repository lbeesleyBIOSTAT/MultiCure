
#' Baselinehazard_IMP
#' @description The function Baselinehazard_IMP is used to estimate the baseline hazard functions within the MCEM algorithm under COX baseline hazards
#' @param datWIDE A data frame with the following columns: 
#' \itemize{
#' \item Y_R, the recurrence event/censoring time
#' \item delta_R, the recurrence event/censoring indicator
#'\item Y_D, the death event/censoring time 
#' \item delta_D, the death event/censoring indicator
#' \item G, the cure status variable. This takes value 1 for known non-cured, 0 for "known" cured and NA for unknown cure status
#'}
#' @param CovImp  A list with IMPNUM elements containing the imputations of Cov output from MultiCure
#' @param GImp  A matrix with IMPNUM elements containing the imputations of G output from MultiCure
#' @param YRImp  A matrix with IMPNUM elements containing the imputations of Y_R output from MultiCure
#' @param deltaRImp  A matrix with IMPNUM elements containing the imputations of delta_R output from MultiCure
#' @param beta Current estimate of beta
#' @param alpha Current estimate of alpha
#' @param TransCov a list with elements: Trans13, Trans24, Trans14, Trans34, PNonCure. Each list element is a vector containing the names of the variables in Cov to be used in the model for the corresponding transition. 13 is NonCured -> Recurrence, 24 is Cured -> Death, 14 is NonCured -> Death, 34 is Recurrence -> Death. PNonCure contains the names of the covariates for the logistic regression for P(NonCure). 
#' @param ASSUME This variables indicates what equality assumptions we are making regarding the 24 and 14 transitions. The possible options are:
#' \itemize{
#' \item 'SameHazard': Lambda_14(t) = Lambda_24(t)
#' \item 'AllSeparate': No restrictions on Lambda_14(t) and Lambda_24(t)
#' \item  'ProportionalHazard': Lambda_14(t) = Lambda_24(t) exp(Beta0)
#' \item  'SameBaseHaz': Lambda^0_14(t) = Lambda^0_24(t), No restrictions on beta_14 and beta_24
#' }
#' @return EST a list containing the estimates for the baseline hazard function (and cumulative baseline hazard function) for the 1->3, 2->4, 1->4, and 3->4 transitions.
#'
#'
#' @export




Baselinehazard_IMP = function(datWIDE, CovImp,GImp, YRImp,deltaRImp, beta, alpha, TransCov, ASSUME){	
	
	UnequalCens = ifelse(datWIDE$Y_R < datWIDE$Y_D & datWIDE$delta_R == 0 & is.na(datWIDE$G), 1, 0)
	IMPNUM = length(CovImp)
	Nobs = length(datWIDE[,1])
	#################################
	### Estimate Baseline Hazards ###
	#################################
	
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))
	
	XB_beta13TEMP = XB_beta24TEMP  = XB_beta14TEMP  = XB_beta34TEMP  =rep(NA,Nobs*IMPNUM)
	for(i in 1:IMPNUM){
		XB_beta13TEMP[(i-1)*Nobs+(1:Nobs)] =as.numeric(beta[TRANS==1] %*% t(cbind(CovImp[[i]][,TransCov$Trans13])))
		XB_beta24TEMP[(i-1)*Nobs+(1:Nobs)] =as.numeric(beta[TRANS==2] %*% t(cbind(CovImp[[i]][,TransCov$Trans24])))	
		XB_beta14TEMP[(i-1)*Nobs+(1:Nobs)] =as.numeric(beta[TRANS==3] %*% t(cbind(CovImp[[i]][,TransCov$Trans14])))	
		XB_beta34TEMP[(i-1)*Nobs+(1:Nobs)] =as.numeric(beta[TRANS==4] %*% t(cbind(CovImp[[i]][,TransCov$Trans34])))	
	}
	
	Basehaz13 = Baseline_hazard(Y = reshape2::melt(YRImp)$value, d = reshape2::melt(deltaRImp)$value, w = reshape2::melt(GImp)$value/IMPNUM, 
		XB = XB_beta13TEMP, 
		cuts = 'grouped')
	### Set last element of 	Basehaz13 hazard to zero. With usually be zero already, but just as a check (zero because GImp = 0 for Yr0 > tauR)
	Basehaz13[length(Basehaz13[,3]),3] = 0
	Basehaz34 = Baseline_hazard(Y = rep(datWIDE$Y_D,IMPNUM) - reshape2::melt(YRImp)$value, d = rep(datWIDE$delta_D,IMPNUM), 
			w = reshape2::melt(deltaRImp)$value/IMPNUM, XB = XB_beta34TEMP, cuts = 'grouped')		
	### Set last element of 	Basehaz34 hazard to be nonzero. This results in improved performance when we have unequal follow-up
	Basehaz34$hazard[length(Basehaz34$hazard)] = Basehaz34$hazard[length(Basehaz34$hazard)-1]
			
	#When estimating the 1424 baseline hazards, we will artificially censor subjects with G=1 and Unequal follow-up. This may build in more stability in the estimators
	YRMOD = YRImp
	deltaRMOD = deltaRImp
	deltaDMOD = datWIDE$delta_D
	dOtherCausesMOD = ifelse(reshape2::melt(deltaDMOD)$value==1 & reshape2::melt(deltaRMOD)$value == 0, 1, 0)

	if(ASSUME %in% c('SameHazard')){
		Basehaz1424 = Baseline_hazard(Y = reshape2::melt(YRMOD)$value, d = dOtherCausesMOD, 
			w = rep(1,Nobs*IMPNUM)/IMPNUM, XB = XB_beta14TEMP, cuts = 'grouped')
		Basehaz14 = Basehaz1424
		Basehaz24 = Basehaz1424	
	}else if(ASSUME %in% c('SameBaseHaz', 'ProportionalHazard')){
		Basehaz1424 = Baseline_hazardSEPARATE2414(Y = reshape2::melt(YRMOD)$value, d = dOtherCausesMOD, 
			w = reshape2::melt(GImp)$value/IMPNUM, wcomp = (1- reshape2::melt(GImp)$value)/IMPNUM,
			XB1 = XB_beta14TEMP, XB2 = XB_beta24TEMP, cuts = 'grouped')		
		Basehaz14 = Basehaz1424
		Basehaz24 = Basehaz1424
	}else{
		Basehaz24 = Baseline_hazard(Y = reshape2::melt(YRMOD)$value, d = dOtherCausesMOD, w = (1-reshape2::melt(GImp)$value)/IMPNUM, 
			XB = XB_beta24TEMP, cuts = 'grouped')
		Basehaz14 = Baseline_hazard(Y = reshape2::melt(YRMOD)$value, d = dOtherCausesMOD, w = reshape2::melt(GImp)$value/IMPNUM, 
			XB = XB_beta14TEMP, cuts = 'grouped')
	}	
	
	Basehaz13$Hazard_lower = as.numeric(sapply(Basehaz13[,1], Baseline_Hazard_Slow, Basehaz13))	
	Basehaz24$Hazard_lower = as.numeric(sapply(Basehaz24[,1], Baseline_Hazard_Slow, Basehaz24))	
	Basehaz14$Hazard_lower = as.numeric(sapply(Basehaz14[,1], Baseline_Hazard_Slow, Basehaz14))	
	Basehaz34$Hazard_lower = as.numeric(sapply(Basehaz34[,1], Baseline_Hazard_Slow, Basehaz34))	

	return(list(Basehaz13, Basehaz24, Basehaz14, Basehaz34))
}





#' BaselineHazard_NOIMP
#' @description The function BaselineHazard_NOIMP is used to estimate the baseline hazard functions within the EM algorithm under COX baseline hazards
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
#' @param beta Current estimate of beta
#' @param alpha Current estimate of alpha
#' @param TransCov a list with elements: Trans13, Trans24, Trans14, Trans34, PNonCure. Each list element is a vector containing the names of the variables in Cov to be used in the model for the corresponding transition. 13 is NonCured -> Recurrence, 24 is Cured -> Death, 14 is NonCured -> Death, 34 is Recurrence -> Death. PNonCure contains the names of the covariates for the logistic regression for P(NonCure). 
#' @param ASSUME This variables indicates what equality assumptions we are making regarding the 24 and 14 transitions. The possible options are:
#' \itemize{
#' \item 'SameHazard': Lambda_14(t) = Lambda_24(t)
#' \item 'AllSeparate': No restrictions on Lambda_14(t) and Lambda_24(t)
#' \item  'ProportionalHazard': Lambda_14(t) = Lambda_24(t) exp(Beta0)
#' \item  'SameBaseHaz': Lambda^0_14(t) = Lambda^0_24(t), No restrictions on beta_14 and beta_24
#' }
#' @param p The current estimate of the E-step weights
#' @return EST a list containing a step function estimate for the CUMULATIVE baseline hazard function for each transition in the following order: 1->3, 2->4, 1->4, 3->4
#'
#' @export




BaselineHazard_NOIMP = function(datWIDE, Cov, beta, alpha, TransCov, ASSUME, p){	
	Nobs = length(datWIDE[,1])	
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
	# Assume that subjects still at risk after the last OBSERVED recurrence event are cured
	Baseline_13 = Baseline_HazardEM(Y = datWIDE$Y_R, d = datWIDE$delta_R, w = p, XB = XB_beta13)
	Baseline_34 = Baseline_HazardEM(Y = datWIDE$Y_D - datWIDE$Y_R, d = datWIDE$delta_D, w = datWIDE$delta_R, XB = XB_beta34)						
	dOtherCauses = ifelse(datWIDE$delta_D==1 & datWIDE$delta_R == 0, 1, 0)
	if(ASSUME %in% c('SameHazard')){
		Baseline_1424 = Baseline_HazardEM(Y = datWIDE$Y_R, d = dOtherCauses, w = rep(1,Nobs), XB = XB_beta14)
		Baseline_14 = Baseline_1424
		Baseline_24 = Baseline_1424	
	}else if(ASSUME %in% c('SameBaseHaz', 'ProportionalHazard')){
		Baseline_1424 = Baseline_HazardSEPARATE2414EM(Y = datWIDE$Y_R, d = dOtherCauses, w = p, 
			wcomp = (1-p), XB1 = XB_beta14, XB2 = XB_beta24)		
		Baseline_14 = Baseline_1424
		Baseline_24 = Baseline_1424
	}else{
		Baseline_24 = Baseline_HazardEM(Y = datWIDE$Y_D, d = dOtherCauses, w = (1-p), XB = XB_beta24)
		Baseline_14 = Baseline_HazardEM(Y = datWIDE$Y_R, d = dOtherCauses, w = p, XB = XB_beta14)
	}		
	Haz_13 = stepfun(x=Baseline_13$event_times, y = c(0,Baseline_13$Hazard), right = F)
	Haz_24 = stepfun(x=Baseline_24$event_times, y = c(0,Baseline_24$Hazard), right = F)
	Haz_14 = stepfun(x=Baseline_14$event_times, y = c(0,Baseline_14$Hazard), right = F)
	Haz_34 = stepfun(x=Baseline_34$event_times, y = c(0,Baseline_34$Hazard), right = F)

	return(list(Haz_13, Haz_24, Haz_14, Haz_34))
}






####################################################
### Functions called within Baselinehazard_NOIMP ###
####################################################



 #' @export


Baseline_hazard_iEM = function(TIME,Y,d,w,XB){
	at_risk = which(Y >= TIME)
	had_event = which(Y == TIME & d == 1)
	hazard = ifelse(sum(w[had_event]) ==0, 0, sum(w[had_event]) /sum(w[at_risk]*exp(XB[at_risk])))
	return(hazard)
}


 #' @export

Baseline_HazardEM = function(Y,d,w,XB){
	event_times = sort(unique(Y[d==1]))	
	hazard_save = sapply(event_times, Baseline_hazard_iEM, Y,d,w,XB)
	Hazard = cumsum(hazard_save)
	return(data.frame(event_times, Hazard, hazard = hazard_save))
}
#' @export

 
Baseline_hazardSEPARATE2414_iEM = function(TIME,Y,d,w,wcomp, XB1, XB2){
	at_risk = which(Y >= TIME)
	had_event = which(Y == TIME & d == 1)
	hazard = ifelse(sum(w[had_event]+wcomp[had_event]) == 0, 0,(sum(w[had_event]+wcomp[had_event])) /(sum(w[at_risk]*exp(XB1[at_risk])) + sum(wcomp[at_risk]*exp(XB2[at_risk])))  )
	return(hazard)
}

#' @export

 
Baseline_HazardSEPARATE2414EM = function(Y,d,w, wcomp ,XB1,XB2){
	event_times = sort(unique(Y[d==1]))	
	hazard_save = sapply(event_times, Baseline_hazardSEPARATE2414_iEM, Y,d,w, wcomp ,XB1, XB2)
	Hazard = cumsum(hazard_save)
	return(data.frame(event_times, Hazard, hazard = hazard_save))
}





##################################################
### Functions called within Baselinehazard_IMP ###
##################################################



### These next two functions estimate baseline hazards for a single transition ###

 #' @export
Baseline_hazard_i = function(TIME,Y,d,w,XB, cutoffs){
	index = which(cutoffs == TIME)
	#MAXTIME = max(Y)
	TIME_LOWER = ifelse(index == 1, 0, cutoffs[index-1])
	TIME_UPPER = cutoffs[index]
	
	#Interval open on right, closed on left [t_index-1, t_index)
	at_risk = which(Y >= TIME_LOWER)
	had_event = which(Y >= TIME_LOWER & Y < TIME_UPPER & d==1)

	NUM = sum(w[had_event])
	DENOM = sum(  (  w*exp(XB)*(pmin(Y,rep(TIME_UPPER, length(Y)))-TIME_LOWER)  )[at_risk[!(at_risk %in% had_event)] ])# +  sum((  w*exp(XB)*(TIME_UPPER-TIME_LOWER)  )[had_event]  )
	hazard = ifelse(NUM==0, 0, NUM/DENOM )
	hazard = ifelse(is.infinite(hazard),0,hazard)
	return(hazard)
}
 
#' @export

Baseline_hazard = function(Y,d,w,XB, cuts = 'events'){	
	MAXTIME = max(Y)
	Y = Y[w!=0]
	d = d[w!=0]
	XB = XB[w!=0]
	w = w[w!=0]	
	event_times = c(sort(unique(Y[d==1])),MAXTIME)
	if(cuts == 'events'){
		cutoffs = event_times
	}else if(cuts == 'grouped'){
		toosmall = which(diff(event_times)<0.005) #was 0.5
		toosmall = toosmall + 1
		if(length(toosmall)==0){
			cutoffs = event_times
		}else{
			cutoffs = event_times[-toosmall]
		}
	}
	options(warn = -1)	
	hazard_save = sapply(cutoffs, Baseline_hazard_i, Y,d,w, XB, cutoffs)
	options(warn = 1)
	return(data.frame(lower = c(0,cutoffs[1:(length(cutoffs)-1)]), upper = cutoffs, hazard = hazard_save))
}





### These next two functions estimate baseline hazard for the 2->4 and 1->4 transitions when they have the same baselines and perhaps different X*Beta ##

#' @export
Baseline_hazardSEPARATE2414 = function(Y,d,w,wcomp,XB1, XB2, cuts = 'events'){
	MAXTIME = max(Y)
	event_times = c(sort(unique(Y[d==1 & w!=0])),MAXTIME)
	if(cuts == 'events'){
		cutoffs = event_times
	}else if(cuts == 'grouped'){
		toosmall = which(diff(event_times)<0.005)#was 0.5
		toosmall = toosmall + 1
		if(length(toosmall)==0){
			cutoffs = event_times
		}else{
			cutoffs = event_times[-toosmall]
		}	
	}
	options(warn = -1)	
	hazard_save = sapply(cutoffs, Baseline_hazardSEPARATE2414_i, Y,d,w,wcomp, XB1, XB2, cutoffs)
	options(warn = 1)
	return(data.frame(lower = c(0,cutoffs[1:(length(cutoffs)-1)]), upper = cutoffs, hazard = hazard_save))
}

 #' @export

Baseline_hazardSEPARATE2414_i = function(TIME,Y,d,w, wcomp,XB1, XB2, cutoffs){
	index = which(cutoffs == TIME)
	
	TIME_LOWER = ifelse(index == 1, 0, cutoffs[index-1])
	TIME_UPPER = cutoffs[index]
	
	#Interval open on right, closed on left [t_index-1, t_index)
	at_risk = which(Y >= TIME_LOWER)
	had_event = which(Y >= TIME_LOWER & Y < TIME_UPPER & d==1)

	NUM = sum(w[had_event]+wcomp[had_event])
	DENOM = sum(  (  w*exp(XB1)*(pmin(Y,rep(TIME_UPPER, length(Y)))-TIME_LOWER)  )[at_risk[!(at_risk %in% had_event)]] +
			 wcomp*exp(XB2)*(pmin(Y,rep(TIME_UPPER, length(Y)))-TIME_LOWER)[at_risk[!(at_risk %in% had_event)]]  )  
			#+ sum(  (  w*exp(XB1)*(TIME_UPPER-TIME_LOWER)  )[had_event] +
			# wcomp*exp(XB2)*(TIME_UPPER-TIME_LOWER)[had_event]  ) 
	hazard = ifelse(NUM==0, 0, NUM/DENOM )
	hazard = ifelse(is.infinite(hazard),0,hazard)
	return(hazard)
}



### This function estimates the cumulative baseline hazard given the baseline hazard function
#' @export

Baseline_Hazard = function(TIME, Basehaz){
	j = max(which(Basehaz[,1]<=TIME))
	return(Basehaz$Hazard_lower[j] + (TIME-Basehaz$lower[j])*Basehaz$hazard[j])
}



#' @export

Baseline_Hazard_Slow = function(TIME,Basehaz){
	TIME_LOWER = Basehaz[,1]
	TIME_UPPER = Basehaz[,2]
	TIME_LOWER = c(TIME_LOWER,Basehaz[length(Basehaz[,2]),2])
	TIME_UPPER = c(TIME_UPPER,10^6)
	HAZARD = c(Basehaz[,3],0)
	#data.frame(TIME_LOWER, TIME_UPPER, HAZARD)
	j = c(1:length(TIME_LOWER)) #which interval are we in
	Li = j[TIME_LOWER <= TIME & TIME_UPPER > TIME] #interval is closed on the left and open on the right
	Hazard  = sum(  (HAZARD*(TIME_UPPER-TIME_LOWER))[which(j <= (Li-1))]) + (HAZARD*(TIME-TIME_LOWER))[which(j==Li)]
	return(Hazard)
}
