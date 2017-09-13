
#' UNEQUALCENSIMPUTECOXREJECTION
#' @description The function UNEQUALCENSIMPUTECOXREJECTION will perform an imputation algorithm to handle unequal follow-up for recurrence and death. This function can be applied when we assume COX baseline hazards. This function performs imputation through rejection sampling. 

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
#' @details This function cannot be applied when 'T_R' is included in the model for Recurrence -> Death. In this case, use one of the following functions instead UNEQUALCENSIMPUTECOXNESTEDWEIBULL, UNEQUALCENSIMPUTECOXMH or UNEQUALCENSIMPUTECOXINVERSION.
#' @export





### This function should only be used when T_R is not in Trans34
UNEQUALCENSIMPUTECOXREJECTION = function(datWIDE, beta, alpha, ImputeDat, TransCov){
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


	sum(diff(Basehaz34[,1])>0.01)
	BasehazFun_13 = stepfun(x= Basehaz13[,2], y = c(Basehaz13[,3],0), right = F)
	BasehazFun_24 = stepfun(x= Basehaz24[,2], y = c(Basehaz24[,3],0), right = F)
	BasehazFun_14 = stepfun(x= Basehaz14[,2], y = c(Basehaz14[,3],0), right = F)
	BasehazFun_34 = stepfun(x= Basehaz34[,2], y = c(Basehaz34[,3],0), right = F)

	S1_D = exp(-as.numeric(sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz13))*exp(XB_beta13))*
			exp(-as.numeric(sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz14))*exp(XB_beta14))
	S2_D = exp(-as.numeric(sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz24))*exp(XB_beta24))
	h24_D = BasehazFun_24(datWIDE$Y_D)*exp(XB_beta24)
	h14_D = BasehazFun_14(datWIDE$Y_D)*exp(XB_beta14)
	
	H13_R = exp(XB_beta13)*as.numeric(sapply(datWIDE$Y_R,Baseline_Hazard, Basehaz13))
	H13_D = exp(XB_beta13)*as.numeric(sapply(datWIDE$Y_D,Baseline_Hazard, Basehaz13))
	H34_R = exp(XB_beta34)*as.numeric(sapply(datWIDE$Y_D-datWIDE$Y_R,Baseline_Hazard, Basehaz34))
	H34_D = rep(0,Nobs)


	### For use with uniroot (SLOW!)
	Draw_Trunc13 = function(TIME, U,YR, YD){
		
		S_R = exp(-exp(XB_beta13[m])* Baseline_Hazard(YR, Basehaz13))
		S_D = exp(-exp(XB_beta13[m])* Baseline_Hazard(min(YD,TAU), Basehaz13))
		ZERO =Baseline_Hazard(TIME, Basehaz13) + (log( U*(S_R - S_D) + S_D)*(1/exp(XB_beta13[m])))
		return(ZERO)
	}
	Draw_Trunc34 = function(TIME, U,YR, YD){
		
		S_R = exp(-exp(XB_beta34[m])* Baseline_Hazard(YD-YR, Basehaz34))
		S_D = exp(-exp(XB_beta34[m])* Baseline_Hazard(YD - min(YD,TAU), Basehaz34))
		ZERO =Baseline_Hazard(YD-TIME, Basehaz34) + (log( U*(S_D - S_R) + S_R)*(1/exp(XB_beta34[m])))
		return(ZERO)
	}


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


	DrawVALWRAPPER = function(s){
		m = INDICES[s]	
		if(datWIDE$delta_D[m]==0){
				ACCEPT = FALSE
				counter = 1
				
				#This will give an approximate K (should even be an over-estimate)		
				timesLOWER = sort(unique(c(Basehaz14[,1], datWIDE$Y_D[m]-Basehaz34[,1], datWIDE$Y_R[m], datWIDE$Y_D[m]) ))
			  	timesLOWER = timesLOWER[timesLOWER>=datWIDE$Y_R[m] & timesLOWER <=datWIDE$Y_D[m] ]
				K = max(((exp(XB_beta34[m])*BasehazFun_34(datWIDE$Y_D[m]-timesLOWER))^datWIDE$delta_D[m])*
					exp(-exp(XB_beta34[m])*as.numeric(sapply(datWIDE$Y_D[m]-timesLOWER,Baseline_Hazard, Basehaz34 )))*
					exp(-exp(XB_beta14[m])*as.numeric(sapply(timesLOWER,Baseline_Hazard, Basehaz14 ))))
						
				while(ACCEPT == FALSE){				
					U1 = runif(n=1, min = 0, max = 1)
					draw = stats::uniroot(Draw_Trunc13, interval = c(datWIDE$Y_R[m], min(datWIDE$Y_D[m],TAU)), 
						U1,datWIDE$Y_R[m], min(datWIDE$Y_D[m],TAU))$root
					if(abs(draw-datWIDE$Y_D[m]) < 0.001){draw = draw - 0.001}
					H14_draw = exp(XB_beta14[m])*as.numeric(sapply(draw,Baseline_Hazard, Basehaz14 )) #14
					H34_draw = exp(XB_beta34[m])*as.numeric(sapply(datWIDE$Y_D[m]-draw,Baseline_Hazard, Basehaz34 ))	 
					h34_draw = exp(XB_beta34[m])*BasehazFun_34(datWIDE$Y_D[m]-draw)  #34
					Surv3_draw = exp(-H34_draw) #34																				
					
					U2 = runif(n=1, min = 0, max = 1)
					if(U2 <= ((1/K)*exp(-H14_draw))*Surv3_draw*(h34_draw^datWIDE$delta_D[m])){
						ACCEPT = TRUE
					}	
					counter = counter + 1
					print(counter)					
				}
		}else if(datWIDE$delta_D[m]==1){
			ACCEPT = FALSE
			counter = 1		
						
			#This will give an approximate K (should even be an over-estimate)		
			timesLOWER = sort(unique(c(Basehaz14[,1], Basehaz13[,1], datWIDE$Y_R[m], datWIDE$Y_D[m]) ))
		  	timesLOWER = timesLOWER[timesLOWER>=datWIDE$Y_R[m] & timesLOWER <=datWIDE$Y_D[m] ]
			K = max((exp(XB_beta13[m])*BasehazFun_13(timesLOWER))*exp(-exp(XB_beta13[m])*as.numeric(sapply(timesLOWER,Baseline_Hazard, Basehaz13 )))*
			exp(-exp(XB_beta14[m])*as.numeric(sapply(timesLOWER,Baseline_Hazard, Basehaz14 ))))

			while(ACCEPT == FALSE){
				
				#Draw from truncated weibull distribution
				U1 = runif(n=1, min = 0, max = 1)
				draw = stats::uniroot(Draw_Trunc34, interval = c(datWIDE$Y_R[m],min(datWIDE$Y_D[m],TAU)), 
					U1,datWIDE$Y_R[m], min(datWIDE$Y_D[m],TAU))$root
				if(abs(draw-datWIDE$Y_D[m]) < 0.001){draw = draw - 0.001}
				H14_draw = exp(XB_beta14[m])*as.numeric(sapply(draw,Baseline_Hazard, Basehaz14 )) #14	
				H13_draw = exp(XB_beta13[m])*as.numeric(sapply(draw,Baseline_Hazard, Basehaz13 )) #13
				h13_draw = exp(XB_beta13[m])*BasehazFun_13(draw) #13
				Surv13_draw = exp(-H13_draw) #13
																	
				U2 = runif(n=1, min = 0, max = 1)
				if(U2 <= ((1/K)*exp(-H14_draw))*Surv13_draw*h13_draw){
					ACCEPT = TRUE
				}	
				counter = counter + 1
				print(counter)					
			}
		}#End of ifelse
		if(draw >= datWIDE$Y_D[m] ){draw = datWIDE$Y_D[m] - (datWIDE$Y_D[m]/1000)}
		if(draw <= datWIDE$Y_R[m] ){draw = datWIDE$Y_R[m] +  (datWIDE$Y_R[m]/1000)}
		#print(m)	
		return(draw)
	}
	DRAWS = sapply(as.numeric(c(1:length(INDICES))), DrawVALWRAPPER)
	YRImp[INDICES] = DRAWS	

	return(list(deltaRImp, YRImp))
}



