
#' SimulateMultiCure
#' @description This function will create a simulated dataset under a multistate cure model.
#'
#' @param type Indicates whether we should simulate data with missingness only in cure status, missingness in cure status and covariates, or missingness in cure status and unequal censoring. Possible values include "No Missingness", "CovariateMissingness" and "UnequalCensoring".
#'
#' @return OUT This is a data frame containing
#' \itemize{
#' \item Y_R: Event/Censoring Time for Recurrence
#' \item delta_R: Event/Censoring Indicator for Recurrence
#' \item Y_D: Event/Censoring Time for Death
#' \item delta_D: Event/Censoring Indicator for Death
#' \item G: Non-Cure Status
#' \item X1: First Covariate
#' \item X2: Second Covariate
#'}
#'
#' @examples
#' SIMS = SimulateMultiCure(type = "NoMissingness")
#' @export



SimulateMultiCure = function(type = 'NoMissingness'){
	type = match.arg(type, choices = c('NoMissingness','CovariateMissingness', 'UnequalCensoring'))
	mySample = function(Y){
		return(sample(x=c(0,1), size = 1, prob = c(1-Y,Y)))
	}
	expit = function(x){
		return(exp(x)/(1+exp(x)))
	}
	
	Nobs = 2000
	censoringInterval = c(10,80) 
	cureCutoff = 50 #Assume subjects still at risk after 50 are cured
	
	########################
	### Parameter Values ###
	########################
	
	beta_24 = c(0.5,0.5)
	beta_14 = c(0.5,0.5)
	beta_34 = c(0.5,0.5)
	beta_13 = c(0.5,0.5)
	alpha_NC = c(0.5, 0.5, 0.5)
	#H0 = scale*t^shape
	scale_13 = 0.005
	shape_13 = 2
	scale_14 = 0.002
	shape_14 = 1.4
	scale_34 = 0.08
	shape_34 = 1.9
	scale_24 = 0.002
	shape_24 = 1.4
	scale_C = 0.01
	shape_C = 1.5
	C = Y_R = Y_D = delta_R = delta_D = X1 = X2 = T_14 = T_13 = T_24 = T_34 = G = GTrue = X2Missing = matrix(NA,Nobs,1)
	
	###########################
	### Generate Covariates ###
	###########################
	
	Z = MASS::mvrnorm(n=Nobs, mu = c(0,0), Sigma = cbind(c(1, 0.5), c(0.5, 1)))
	X1 = Z[,1]
	X2 = Z[,2]
	
	############################
	### Generate Cure Status ###
	############################
	
	Prob = expit(alpha_NC[1] + alpha_NC[2]*X1 + alpha_NC[3]*X2)
	GTrue = sapply(Prob, mySample) 
	
	###############################
	### Generate Censoring Time ###
	###############################
	
	C = runif(n=Nobs, min=censoringInterval[1], max=censoringInterval[2])
	
	#########################
	### Generate Outcomes ###
	#########################
	
	XB_14 = beta_14[1]*X1 + beta_14[2]*X2
	XB_13 = beta_13[1]*X1 + beta_13[2]*X2
	XB_34 = beta_34[1]*X1 + beta_34[2]*X2
	XB_24 = beta_24[1]*X1 + beta_24[2]*X2
	
	for(i in 1:Nobs)
	{
		if(GTrue[i]==1) #Can have recurrence
		{
			U = runif(n=Nobs, min = 0, max = 1)
			T_14[i] = ( (-log(U[i]))  /  (scale_14*exp(XB_14[i]))  )^(1/shape_14)
			U = runif(n=Nobs, min = 0, max = 1)
			T_13[i] = ( (-log(U[i]))  /  (scale_13*exp(XB_13[i]))  )^(1/shape_13)
			U = runif(n=Nobs, min = 0, max = 1)
			T_34[i] = ( (-log(U[i]))  /  (scale_34*exp(XB_34[i]))  )^(1/shape_34)
			T_24[i] = 10^6				
			if(T_14[i]<T_13[i]) #die or censored before recurrence
			{
				Y_D[i] = min(T_14[i],C[i])
				delta_D[i] = ifelse(which.min(c(T_14[i],C[i]))==1, 1, 0)	
				Y_R[i] = Y_D[i]
				delta_R[i] = 0
			}else{
				if(T_13[i]<C[i]) #observe a recurrence
				{
					Y_R[i] = T_13[i]
					delta_R[i] = 1
					Y_D[i] = min(T_13[i]+T_34[i],C[i])
					delta_D[i] = ifelse(which.min(c(T_13[i]+T_34[i],C[i]))==1, 1, 0)							
				}else{
					Y_R[i] = C[i]
					Y_D[i] = C[i]
					delta_R[i] = 0
					delta_D[i] = 0				
				}		
			}	
		}else if(GTrue[i]==0)
		{
			U = runif(n=Nobs, min = 0, max = 1)
			T_24[i] = ( (-log(U[i]))  /  (scale_24*exp(XB_24[i]))  )^(1/shape_24)
			T_14[i] = 10^6
			T_13[i] = 10^6
			T_34[i] = 10^6		
			Y_D[i] = min(T_24[i],C[i])
			delta_D[i] = ifelse(which.min(c(T_24[i],C[i]))==1, 1, 0)	
			Y_R[i] = Y_D[i]
			delta_R[i] = 0
		}
	}
	
	####################################
	### Impose Covariate Missingness ### 
	####################################
     
	missP2 = rep(0.3, Nobs)
	NmissInd = sapply(missP2, mySample) 
	X2Miss = ifelse(NmissInd == 1, rep(NA,Nobs), X2)
	
	################################
	### Impose Unequal Follow-up ###
	################################
	
	#(Assume first 500 subjects have equal censoring. Early censoring is the minimum of the censoring time for death--perhaps administrative censoring-- and a random dropout rate )
	U = runif(n=Nobs, min=0, max=1)
	temp = runif(n=Nobs, min=10, max=40)
	ID = c(1:Nobs)
	C_early = ifelse(C<temp  | ID < 750, C, temp)
	delta_R_early = ifelse(Y_R >= C_early, rep(0,Nobs), delta_R)
	Y_R_early = ifelse(Y_R >= C_early, C_early, Y_R)
	
	############################
	### Observed Cure Status ###
	############################
	
	G = ifelse(delta_R==1, 1, ifelse(Y_D>cureCutoff & Y_R>cureCutoff & delta_R==0, 0, NA))
	G_early = ifelse(delta_R_early==1, 1, ifelse(Y_D>cureCutoff & Y_R_early>cureCutoff & delta_R_early==0, 0, NA))

	##############		
	### Return ###
	##############
	if(type == 'NoMissingness'){
		return(data.frame(Y_R, delta_R, Y_D, delta_D, G, X1, X2))	
	}else if(type == 'CovariateMissingness'){
		return(data.frame(Y_R, delta_R, Y_D, delta_D, G, X1, X2 = X2Miss))					
	}else{
		return(data.frame(Y_R = Y_R_early, delta_R = delta_R_early, Y_D, delta_D, G = G_early, X1, X2))		
	}
}#end function


