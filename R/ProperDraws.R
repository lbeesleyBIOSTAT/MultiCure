

#' @export

ProperDraws = function( datWIDE,Cov,CovImp, GImp, YRImp, deltaRImp, COVIMPUTEFUNCTION = NULL,  COVIMPUTEINITIALIZE = NULL,
			UNEQUALCENSIMPUTE = NULL, ASSUME = 'SameHazard', TransCov, BASELINE, PENALTY = 'None',POSTITER = 5){
	
	##################
	### Initialize ###
	##################
	
	Nobs = length(datWIDE[,1])
	ASSUME = match.arg(ASSUME, choices = c('SameHazard', 'AllSeparate', 'SameBaseHaz',  'ProportionalHazard'))
	BASELINE = match.arg(BASELINE, choices = c('weib','cox'))
	PENALTY = match.arg(PENALTY, choices = c('None', 'Ridge', 'Lasso'))
	UnequalCens = ifelse(datWIDE$Y_R < datWIDE$Y_D & datWIDE$delta_R == 0 & is.na(datWIDE$G), 1, 0)
	CovMissing = apply(Cov,1:2,is.na)
	NEEDTOIMPUTE = TRUE

	IMPNUM = length(CovImp)
	if(sum(UnequalCens) != 0 & is.null(UNEQUALCENSIMPUTE)){
		if(BASELINE == 'weib' & 'T_R' %in% TransCov$Trans34){
			UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBINVERSION
		}else if(BASELINE == 'weib' & !('T_R' %in% TransCov$Trans34)){
			UNEQUALCENSIMPUTE = UNEQUALCENSIMPUTEWEIBREJECTION
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

	
	###############################################
	### Initialize Missing Values for YRImpSAVE ###
	###############################################
	
	# This is only needed for unequal censoring imputation via Metropolis-Hastings (OutcomeImputationCOXMH)
	YRImpSAVE = YRImp
	for(i in 1:IMPNUM){
		TAU_R = max(datWIDE$Y_R[datWIDE$delta_R==1])
		IMPUTEYR = (UnequalCens ==1)
		MIN = datWIDE$Y_R[IMPUTEYR]
		MAX = pmin(datWIDE$Y_D[IMPUTEYR], TAU_R)
		MIN = ifelse(MIN >= MAX, MAX, MIN)
		U = apply(cbind(MIN, MAX),1, mHPropose) 
		YRImpSAVE[IMPUTEYR,i] = ifelse(deltaRImp[IMPUTEYR,i]==1, YRImp[IMPUTEYR,i], U)	
	}
	
			
	#################################
	### Obtain Proper Imputations ###	
	#################################
	
	iter = 1
	while(iter <= POSTITER)
	{
		for(i in 1:IMPNUM){
			
			################################################
			### Draw Bootstrap Sample of Imputed Dataset ###
			################################################
			
			whichboot = sample(x=c(1:Nobs), size = Nobs, replace = TRUE, prob = rep(1/Nobs, Nobs))
			ImputeDatBOOT = list(UnequalCens=UnequalCens[whichboot], CovMissing=CovMissing[whichboot,], 
							CovImp= list( CovImp[[i]][whichboot,]), GImp= matrix(GImp[whichboot,i]), 
							YRImp= matrix(YRImp[whichboot,i]), deltaRImp= matrix(deltaRImp[whichboot,i]), YRImpSAVE = matrix(YRImpSAVE[whichboot,i])  )	
			ImputeDat = list(UnequalCens= UnequalCens, CovMissing= CovMissing, CovImp  = list(CovImp[[i]]), GImp  = matrix(GImp[,i]), YRImp  = matrix(YRImp[,i]), 
							deltaRImp  = matrix(deltaRImp[,i]), YRImpSAVE = matrix(YRImpSAVE[,i]))		
			
			########################################################
			### Fit model to bootstrap sample and then re-impute ###
			########################################################
			
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





