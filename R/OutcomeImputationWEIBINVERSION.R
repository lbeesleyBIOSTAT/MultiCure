#' UNEQUALCENSIMPUTEWEIBINVERSION
#' @description The function UNEQUALCENSIMPUTEWEIB_TR will perform an imputation algorithm to handle unequal follow-up for recurrence and death. This function can be applied when we assume WEIBULL baseline hazards and 'T_R' is included in the model for Recurrence -> Death. Rather than imputing new values of YR using rejection sampling, this function performs imputation through inverting the survival function of the target distribution. 

#' @param datWIDE defined as in MultiCure
#' @param beta A vector containing the most recent estimates of beta
#' @param alpha A vector containing the most recent estimates of alpha
#' @param scale A vector containing the most recent estimates of scale
#' @param shape A vector containing the most recent estimates of shape
#' @param ImputeDat This is a list with the following elements:
#' \itemize{
#' \item UnequalCens: A vector taking value 1 if the subject has unequal follow-up. Note: If subject is assumed cured in datWIDE, they are listed as UnequalCens = 0.
#' \item CovMissing: A matrix indicating which elements of Cov are missing. Not needed for this imputation.
#' \item CovImp: A list containing a single imputation of Cov
#' \item GImp: A vector with a recent single imputation of G
#' \item YRImp: A vector with a recent single imputation of Y_R
#' \item deltaRImp: A vector with a recent single imputation of delta_R
#' \item y: The integral of the target kernel over Yr0 to Yd
#' }
#' @param TransCov defined as in MultiCure
#'
#' @return a list containing 
#' \itemize{
#' \item [[1]]: deltaRImp, A single imputation of delta_R
#' \item [[2]]: YRImp, A single imputation of Y_R
#'}
#' @export



UNEQUALCENSIMPUTEWEIBINVERSION = function(datWIDE, beta, alpha, scale, shape, ImputeDat, TransCov){
	UnequalCens = ImputeDat[[1]]
	CovImp = as.data.frame(ImputeDat[[3]])
	GImp = ImputeDat[[4]]
	YRImp = ImputeDat[[5]]
	deltaRImp = ImputeDat[[6]]
	y = ImputeDat[[7]]	
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
	
	S1_D = exp(- (scale[1]*((datWIDE$Y_D)^shape[1]) ) *exp(XB_beta13))*
			exp(-(scale[3]*((datWIDE$Y_D)^shape[3]) )*exp(XB_beta14))
	h14_D = (scale[3]*shape[3]*((datWIDE$Y_D)^(shape[3]-1))  )*exp(XB_beta14)

	YRImp = ifelse(GImp==0,datWIDE$Y_D, ifelse(GImp==1 & UnequalCens == 0,datWIDE$Y_R,rep(NA,Nobs) ))
	deltaRImp = ifelse(GImp==0,rep(0,Nobs), ifelse(GImp==1 & UnequalCens == 0,datWIDE$delta_R,rep(NA,Nobs) ))
	#missing iff GImp = 1 and UnequalCens == 1

	######################
	### Impute Delta R ###
	######################			
	num = y
	denom = (h14_D^datWIDE$delta_D)*S1_D
	ratio = ifelse(num==0,num,num/(num + denom))[GImp==1 & UnequalCens == 1]			
	deltaRImp[GImp==1 & UnequalCens == 1] = apply(matrix(ratio), 1,mSample)   
	YRImp[GImp==1 & UnequalCens == 1 & deltaRImp==0] = datWIDE$Y_D[GImp==1 & UnequalCens == 1 & deltaRImp==0]

	INDICES = which(is.na(YRImp))

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
	DrawVAL = function(TIME, U){		
			g=cubature::adaptIntegrate(Vectorize(fd), lowerLimit = TIME, upperLimit = datWIDE$Y_D[m],maxEval=15)  
			ZERO=(g$integral/y[m])-U
		return(ZERO)
	}	
	
	for(s in 1:length(INDICES)){
		m = INDICES[s]
		U1 = runif(n=1, min = 0, max = 1)
		draw = stats::uniroot(DrawVAL, interval = c(datWIDE$Y_R[m], datWIDE$Y_D[m]),U1)$root
		YRImp[m] = draw	
		#print(s)	
	}#End of for loop
		
	return(list(deltaRImp, YRImp))
}


