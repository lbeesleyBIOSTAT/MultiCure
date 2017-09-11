#' UNEQUALCENSIMPUTEWEIB
#' @description The function UNEQUALCENSIMPUTEWEIB will perform a rejection sampling algorithm to handle unequal follow-up for recurrence and death. This function is the default imputation method when we assume WEIBULL baseline hazards. We allow the user to specify the function UNEQUALCENSIMPUTE in MultiCure if desired, but we provide default method. When baseline hazards are WEIBULL, the input/output of user-written functions must match the input/output of this default function. 

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
#' @details This function cannot be applied when 'T_R' is included in the model for Recurrence -> Death. In this case, use the function UNEQUALCENSIMPUTEWEIB_TR.
#' @author Lauren J Beesley, \email{lbeesley@umich.edu}
#' @export







UNEQUALCENSIMPUTEWEIB = function(datWIDE, beta, alpha, scale, shape, ImputeDat, TransCov){
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
	
	H13_R = scale[1]*exp(XB_beta13)*datWIDE$Y_R^(shape[1]) 
	H13_D = scale[1]*exp(XB_beta13)*datWIDE$Y_D^(shape[1])		
	H34_R = scale[4]*exp(XB_beta34)*(datWIDE$Y_D-datWIDE$Y_R)^(shape[4]) 
	H34_D = rep(0,Nobs)
	S1_D = exp(- (scale[1]*((datWIDE$Y_D)^shape[1]) ) *exp(XB_beta13))*
			exp(-(scale[3]*((datWIDE$Y_D)^shape[3]) )*exp(XB_beta14))
	S2_D = exp(-(scale[2]*((datWIDE$Y_D)^shape[2]) )*exp(XB_beta24))
	h24_D = (scale[2]*shape[2]*((datWIDE$Y_D)^(shape[2]-1))      )*exp(XB_beta24)
	h14_D = (scale[3]*shape[3]*((datWIDE$Y_D)^(shape[3]-1))  )*exp(XB_beta14)

		#which(GImp == 1 & UnequalCens==1 & datWIDE$delta_D==0)

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

		fd<-function(v){
			S1 = exp(-(scale[1]*((v)^shape[1]) ) *exp(XB_beta13[m]))*exp(-(scale[3]*((v)^shape[3]) )*exp(XB_beta14[m]))
			S3 = exp(-(scale[4]*((datWIDE$Y_D[m]-v)^shape[4]) )*exp(XB_beta34[m]))		
			h13 = (scale[1]*shape[1]*((v)^(shape[1]-1))  )*exp(XB_beta13[m])
			h34 = ifelse(v==datWIDE$Y_D[m],0,(scale[4]*shape[4]*((datWIDE$Y_D[m]-v)^(shape[4]-1))  )*exp(XB_beta34[m])	)
			return(h13*S1*S3*((h34)^datWIDE$delta_D[m]))   
		}
	INDICES = which(is.na(YRImp))
	#which(datWIDE$delta_D==0 & c(1:Nobs)%in% INDICES)
	for(s in 1:length(INDICES)){
		m = INDICES[s]
		if(datWIDE$delta_D[m]==0){
				ACCEPT = FALSE
				counter = 1
				while(ACCEPT == FALSE){				
					#Draw from truncated weibull distribution
					U1 = runif(n=1, min = 0, max = 1)
					draw = (-log( U1*(exp(-H13_R[m]) - exp(-H13_D[m])) + exp(-H13_D[m])    )*
						(1/(exp(XB_beta13[m])*scale[1]))  )^(1/shape[1])
					H14_draw = scale[3]*exp(XB_beta14[m])*draw^(shape[3]) #14
					H34_draw = scale[4]*exp(XB_beta34[m])*(datWIDE$Y_D[m]-draw)^(shape[4]) #34
					h34_draw = shape[4]*scale[4]*exp(XB_beta34[m])*(datWIDE$Y_D[m]-draw)^(shape[4]-1) #34
					Surv3_draw = exp(-H34_draw) #34									
					
					### K is the maximum value  lambda_34(Yd-t)^delta_d S_3(Yd-t) S_1(t) takes in Yr^0 to Yd
					U = c(runif(n = 1000, min = datWIDE$Y_R[m], max = datWIDE$Y_D[m]))
					U = sort(U)
					K = ((shape[4]*scale[4]*exp(XB_beta34[m])*((datWIDE$Y_D[m] - U)^(shape[4]-1)))^datWIDE$delta_D[m])*
						exp(-scale[4]*exp(XB_beta34[m])*((datWIDE$Y_D[m] - U)^shape[4]))*
						exp(-scale[3]*exp(XB_beta14[m])*(U^shape[3]))
					K = max(K)
					
					# y = fd(U)
					# z = K* (dweibull(U,shape =shape[1], scale = ((scale[1]*exp(XB_beta13[m]))^(-1/shape[1])) ) ) 	
					# plot(U,y, ylim = c(0,max(z)), xlab = 'Possible Tr Values', ylab = 'Kernel Value', 
						 # main = 'Comparing Target and Dominating Kernels', cex.main = 0.8, cex.lab = 0.8)
					# lines(U,z, col = 'red', lwd = 2)
					# legend(x='topright', fill = c('red', 'black'), legend=c('Dominating Kernel', 'Target Kernel'), cex = 0.8)
									
					U2 = runif(n=1, min = 0, max = 1)
					if(U2 <= ((1/K)*exp(-H14_draw))*Surv3_draw*(h34_draw^datWIDE$delta_D[m])){
						ACCEPT = TRUE
					}	
					counter = counter + 1
					if(counter == 1000){
						#print('High counter')
						ACCEPT = TRUE
					}
					#print(counter)					
				}
				YRImp[m] = draw
		}else if(datWIDE$delta_D[m]==1){
			ACCEPT = FALSE
			counter = 1
			while(ACCEPT == FALSE){
				#Draw from truncated weibull distribution
					U1 = runif(n=1, min = 0, max = 1)
					draw = (-log( U1*(exp(-H34_D[m]) - exp(-H34_R[m])) + exp(-H34_R[m])    )*
						(1/(exp(XB_beta34[m])*scale[4]))  )^(1/shape[4])
						
					# draw = (-log( -U1*(1 - exp(-H34_R[m])) + 1    )*
						# (1/(exp(XB_beta34[m])*scale[4]))  )^(1/shape[4]) #CHECKED THAT THESE ARE EQUIVLENT!!!!!!!!

					draw = datWIDE$Y_D[m] - draw
					H14_draw = scale[3]*exp(XB_beta14[m])*(draw^(shape[3])) #14	
					H13_draw = scale[1]*exp(XB_beta13[m])*(draw)^(shape[1]) #13
					h13_draw = shape[1]*scale[1]*exp(XB_beta13[m])*(draw^(shape[1]-1)) #13
					Surv13_draw = exp(-H13_draw) #13
					
					### K is the maximum value the function lambda_34(Yd-t)^delta_d S_3(Yd-t) S_1(t) takes in Yr^0 to Yd
					U = c(runif(n = 1000, min = datWIDE$Y_R[m], max = datWIDE$Y_D[m]))
					U = sort(U)
					K = ((shape[1]*scale[1]*exp(XB_beta13[m])*(U)^(shape[1]-1)))*
						exp(-scale[1]*exp(XB_beta13[m])*(U^shape[1]))*
						exp(-scale[3]*exp(XB_beta14[m])*(U^shape[3]))
					K = max(K)										
					
					# x  = fd(U)
					# z = K*(dweibull(datWIDE$Y_D[m]-U,shape =shape[4], 
						# scale = (scale[4]*exp(XB_beta34[m]))^(-1/shape[4]) ) ) 
					# plot(U,x, ylim = c(0,max(z)), xlab = 'Possible Tr Values', ylab = 'Kernel Value', 
						# main = 'Comparing Target and Dominating Kernels', cex.main = 0.8, cex.lab = 0.8)
					# lines(U,z, col = 'red', lwd = 2)
					# legend(x='bottomleft', fill = c('red', 'black'), legend=c('Dominating Kernel', 'Target Kernel'), cex=0.8)
					
					U2 = runif(n=1, min = 0, max = 1)
					if(U2 <= ((1/K)*exp(-H14_draw))*Surv13_draw*h13_draw){
						ACCEPT = TRUE
					}	
					counter = counter + 1
					#print(counter)	
					if(counter == 1000){
						print('High counter')
						ACCEPT = TRUE
					}
			}							
			YRImp[m] = draw
		}#End of ifelse
		
	}#End of for loop

	return(list(deltaRImp, YRImp))
}

