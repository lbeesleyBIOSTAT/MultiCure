### These functions perform the M Step in the Multistate Cure EM Algorithm
		
 

	#' @export

MStep_WEIB = function(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY){
	### Transform the data into long format
	if(NEEDTOIMPUTE){
		datLONG = CreateLong_MC(datWIDE, ImputeDat)
	}else{
		datLONG = CreateLong(datWIDE, Cov)
	}
	datLONG_sub = datLONG[datLONG$w != 0 & datLONG$time != 0,]		
	### Transform the covariates into the proper format
	Cov_long = subset(datLONG_sub, select = -c(id, from, to, trans, Tstart, Tstop, time, status,w))
	Cov_long_13 = data.frame(Cov_long[,TransCov$Trans13])
	Cov_long_13[!(datLONG_sub$from == 1 & datLONG_sub$to == 3),] = 0
	Cov_long_34 = data.frame(Cov_long[,TransCov$Trans34] )
	Cov_long_34[!(datLONG_sub$from == 3 & datLONG_sub$to == 4),] = 0
	if(ASSUME %in% c('SameHazard',  'ProportionalHazard')){		
		Cov_long_1424 = data.frame(Cov_long[,TransCov$Trans24] )
		Cov_long_1424[!((datLONG_sub$from == 2 | datLONG_sub$from == 1) & datLONG_sub$to == 4),] = 0
	}else if(ASSUME %in% c('AllSeparate', 'SameBaseHaz')){		
		Cov_long_14 = data.frame(Cov_long[,TransCov$Trans14])
		Cov_long_14[!(datLONG_sub$from == 1 & datLONG_sub$to == 4),] = 0
		Cov_long_24 = data.frame(Cov_long[,TransCov$Trans24] )
		Cov_long_24[!(datLONG_sub$from == 2 & datLONG_sub$to == 4),] = 0
	}		
	if(PENALTY == 'None'){
		TRANS = I
	}else if(PENALTY == 'Ridge'){
		TRANS = survival::ridge
	}else{
		stop('Invalid PENALTY Argument: Only Implemented for Ridge Penalty and No Penalty')
	}	
	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)	
	
	datLONG_sub$w[datLONG_sub$trans==3 & datLONG_sub$status==1]
	fit13 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)), 
			weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 1))
	mshape13 = 1/fit13$scale
	mscale13 = exp(-coef(fit13)[1]*mshape13)
	coef13 = -coef(fit13)[2:length(coef(fit13))]*mshape13	
	fit34 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_34)),
		  	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 4))				  	  	
	mshape34 = 1/fit34$scale
	mscale34 = exp(-coef(fit34)[1]*mshape34)
	coef34 = -coef(fit34)[2:length(coef(fit34))]*mshape34			
	if(ASSUME %in% c('AllSeparate')){
		fit24 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_24)),
			 weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2))	
		mshape24 = 1/fit24$scale
		mscale24 = exp(-coef(fit24)[1]*mshape24)
		coef24 = -coef(fit24)[2:length(coef(fit24))]*mshape24
		fit14 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_14)),
			 weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 3))
		mshape14 = 1/fit14$scale
		mscale14 = exp(-coef(fit14)[1]*mshape14)
		coef14 = -coef(fit14)[2:length(coef(fit14))]*mshape14
		beta = c(coef13, coef24, coef14, coef34)
		scale = as.numeric(c(mscale13, mscale24, mscale14, mscale34))
		shape = as.numeric(c(mshape13, mshape24, mshape14, mshape34))
	}else if(ASSUME %in% c('SameHazard')){
		fit2414 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_1424)),
		 	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2 | datLONG_sub$trans == 3))	
		mshape2414 = 1/fit2414$scale
		mscale2414 = exp(-coef(fit2414)[1]*mshape2414)
		coef2414 = -coef(fit2414)[2:length(coef(fit2414))]*mshape2414
		beta = c(coef13, coef2414, coef2414, coef34)
		scale = as.numeric(c(mscale13, mscale2414, mscale2414, mscale34)) 
		shape = as.numeric(c(mshape13, mshape2414, mshape2414, mshape34))
	}else if(ASSUME %in% c('SameBaseHaz')){
		fit2414 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_24)) + 
			TRANS(as.matrix(Cov_long_14)),
		 	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2 | datLONG_sub$trans == 3))		
		mshape2414 = 1/fit2414$scale
		mscale2414 = exp(-coef(fit2414)[1]*mshape2414)
		coef2414 = -coef(fit2414)[2:length(coef(fit2414))]*mshape2414
		A2 = length(TransCov$Trans24)
		A3 = length(TransCov$Trans14)
		TRANSITION = c(rep(2,A2), rep(3,A3))
		coef24 = coef2414[TRANSITION==2]
		coef14 = coef2414[TRANSITION==3]	
		beta = c(coef13, coef24, coef14, coef34)
		scale = as.numeric(c(mscale13, mscale2414, mscale2414, mscale34)) 
		shape = as.numeric(c(mshape13, mshape2414, mshape2414, mshape34))
	}else if(ASSUME %in% c('ProportionalHazard')){
		fit2414 = survival::survreg(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_1424)) + as.numeric(datLONG_sub$trans==3),
		 	weights = datLONG_sub$w, dist = 'weib', subset = (datLONG_sub$trans == 2 | datLONG_sub$trans == 3))	
		mshape2414 = 1/fit2414$scale
		mscale2414 = exp(-coef(fit2414)[1]*mshape2414)
		coef2414 = -coef(fit2414)[2:length(coef(fit2414))]*mshape2414
		beta = c(coef13, coef2414[1:(length(coef2414)-1)], coef2414, coef34)
		scale = as.numeric(c(mscale13, mscale2414, mscale2414, mscale34)) 
		shape = as.numeric(c(mshape13, mshape2414, mshape2414, mshape34))
	}

	if(NEEDTOIMPUTE){
		TEMPCOV = c()
		CovImp = ImputeDat[[3]]
		for(k in 1:length(CovImp)){
			TEMPCOV = rbind(TEMPCOV,CovImp[[k]][,TransCov$PNonCure])			
		}
		Out = reshape2::melt(ImputeDat[[4]])$value		
		if(PENALTY == 'None'){
			options(warn=-1)
			fitLogistic = stats::glm(Out~as.matrix(TEMPCOV), family = binomial(link = 'logit'), weights = rep(1/length(CovImp), length(Out)))
			options(warn=1)
			alpha = coef(fitLogistic)			
		}else if(PENALTY %in% c('Ridge', 'Lasso')){
			ALPHA = ifelse(PENALTY == 'Ridge', 0,1)
			cvfit = glmnet::cv.glmnet(x=as.matrix(TEMPCOV), y=Out, family= 'binomial', weights = rep(1/length(CovImp), length(Out)), alpha = ALPHA) 
			RESULTS = coef(cvfit, s = "lambda.1se") 
			alpha = rep(0,1+length(TransCov$PNonCure))
			alpha[as.numeric(slot(RESULTS,'i'))+1] = slot(RESULTS,'x')
			alpha = as.numeric(alpha)
		}
	}else{
		if(PENALTY == 'None'){
			options(warn=-1)		
			fitLogistic = stats::glm(datWIDE$p~as.matrix(Cov[,TransCov$PNonCure]), family = binomial(link = 'logit'))
			options(warn=1)	
			alpha = coef(fitLogistic)
		}else{
			Out = cbind(1-datWIDE$p, datWIDE$p)	
			ALPHA = ifelse(PENALTY == 'Ridge', 0,1)
			cvfit = glmnet::cv.glmnet(x=as.matrix(Cov), y=Out, family= 'binomial', alpha = ALPHA) 
			RESULTS = coef(cvfit, s = "lambda.1se") 
			alpha = rep(0,1+length(TransCov$PNonCure))
			alpha[as.numeric(slot(RESULTS,'i'))+1] = slot(RESULTS,'x')
			alpha = as.numeric(alpha)		
		}
	}
	alpha[is.na(alpha)] = 0
	return(list(beta = as.numeric(beta),alpha = as.numeric(alpha), scale = as.numeric(scale), shape = as.numeric(shape)))
}

	


#' @export

MStep_COX = function(datWIDE, Cov, ImputeDat, ASSUME, TransCov, NEEDTOIMPUTE, PENALTY){
	### Transform the data into long format
	if(NEEDTOIMPUTE){
		datLONG = CreateLong_MC(datWIDE, ImputeDat)
	}else{
		datLONG = CreateLong(datWIDE, Cov)
	}
	datLONG_sub = datLONG[datLONG$w != 0 & datLONG$time != 0,]
	### Transform the covariates into the proper format
	Cov_long = subset(datLONG_sub, select = -c(id, from, to, trans, Tstart, Tstop, time, status,w))
	Cov_long_13 = data.frame(Cov_long[,TransCov$Trans13])
	Cov_long_13[!(datLONG_sub$from == 1 & datLONG_sub$to == 3),] = 0
	Cov_long_34 = data.frame(Cov_long[,TransCov$Trans34] )
	Cov_long_34[!(datLONG_sub$from == 3 & datLONG_sub$to == 4),] = 0
	if(ASSUME %in% c('SameHazard',  'ProportionalHazard')){		
		Cov_long_1424 = data.frame(Cov_long[,TransCov$Trans24] )
		Cov_long_1424[!((datLONG_sub$from == 2 | datLONG_sub$from == 1) & datLONG_sub$to == 4),] = 0
	}else if(ASSUME %in% c('AllSeparate', 'SameBaseHaz')){		
		Cov_long_14 = data.frame(Cov_long[,TransCov$Trans14])
		Cov_long_14[!(datLONG_sub$from == 1 & datLONG_sub$to == 4),] = 0
		Cov_long_24 = data.frame(Cov_long[,TransCov$Trans24] )
		Cov_long_24[!(datLONG_sub$from == 2 & datLONG_sub$to == 4),] = 0
	}		

	A1 = length(TransCov$Trans13)
	A2 = length(TransCov$Trans24)
	A3 = length(TransCov$Trans14)
	A4 = length(TransCov$Trans34)
	if(PENALTY == 'None'){
		TRANS = I
		strata = survival::strata
		if(ASSUME %in% c('AllSeparate')){
			fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
					TRANS(as.matrix(Cov_long_24)) + TRANS(as.matrix(Cov_long_14)) + TRANS(as.matrix(Cov_long_34)) +
					strata(datLONG_sub$trans), weights = datLONG_sub$w)
			beta = coef(fitCox)				
		}else if(ASSUME %in% c('SameHazard')){
			Trans_MERGE = datLONG_sub$trans
			Trans_MERGE[Trans_MERGE==3] = 2
			fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
					TRANS(as.matrix(Cov_long_1424)) + TRANS(as.matrix(Cov_long_34)) +
					strata(Trans_MERGE), weights = datLONG_sub$w)
			beta_short = coef(fitCox)
			TRANSITION = c(rep(1,A1), rep(2,A2), rep(4,A4))
			beta = c(beta_short[TRANSITION==1], rep(beta_short[TRANSITION==2],2), beta_short[TRANSITION==4])		
		}else if(ASSUME %in% c('SameBaseHaz')){
			Trans_MERGE = datLONG_sub$trans
			Trans_MERGE[Trans_MERGE==3] = 2
			fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
					TRANS(as.matrix(Cov_long_24)) + TRANS(as.matrix(Cov_long_14)) + TRANS(as.matrix(Cov_long_34)) +
					strata(Trans_MERGE), weights = datLONG_sub$w)
			beta = coef(fitCox)
		}else if(ASSUME %in% c('ProportionalHazard')){
			Trans_MERGE = datLONG_sub$trans
			Trans_MERGE[Trans_MERGE==3] = 2
			fitCox = survival::coxph(survival::Surv(datLONG_sub$time,datLONG_sub$status)~TRANS(as.matrix(Cov_long_13)) + 
					TRANS(as.matrix(Cov_long_1424)) + TRANS(as.matrix(Cov_long_34)) + as.numeric(datLONG_sub$trans==3)+
					strata(Trans_MERGE), weights = datLONG_sub$w)
			beta_short = coef(fitCox)
			TRANSITION = c(rep(1,A1), rep(2,A2), rep(4,A4), 9)
			beta = c(beta_short[TRANSITION==1], rep(beta_short[TRANSITION==2],2), beta_short[TRANSITION==9], beta_short[TRANSITION==4])		
		}
	}else{
		ALPHA = ifelse(PENALTY == 'Ridge', 0,1)		
		if(ASSUME %in% c('AllSeparate')){
			cvfit13 = glmnet::cv.glmnet(x=as.matrix(Cov_long_13[datLONG_sub$trans==1,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==1],datLONG_sub$status[datLONG_sub$trans==1]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==1], alpha = ALPHA) 
			cvfit24 = glmnet::cv.glmnet(x=as.matrix(Cov_long_24[datLONG_sub$trans==2,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==2],datLONG_sub$status[datLONG_sub$trans==2]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==2], alpha = ALPHA) 
			cvfit14 = glmnet::cv.glmnet(x=as.matrix(Cov_long_14[datLONG_sub$trans==3,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==3],datLONG_sub$status[datLONG_sub$trans==3]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==3], alpha = ALPHA) 
			cvfit34 = glmnet::cv.glmnet(x=as.matrix(Cov_long_34[datLONG_sub$trans==4,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==4],datLONG_sub$status[datLONG_sub$trans==4]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==4], alpha = ALPHA) 
			RESULTS13 = coef(cvfit13, s = "lambda.1se") 
			RESULTS24 = coef(cvfit24, s = "lambda.1se") 
			RESULTS14 = coef(cvfit14, s = "lambda.1se") 
			RESULTS34 = coef(cvfit34, s = "lambda.1se") 			
			beta13 = rep(0,A1)
			beta13[as.numeric(slot(RESULTS13,'i'))+1] = slot(RESULTS13,'x')
			beta13 = as.numeric(beta13)
			beta24 = rep(0,A2)
			beta24[as.numeric(slot(RESULTS24,'i'))+1] = slot(RESULTS24,'x')
			beta24 = as.numeric(beta24)			
			beta14 = rep(0,A3)
			beta14[as.numeric(slot(RESULTS14,'i'))+1] = slot(RESULTS14,'x')
			beta14 = as.numeric(beta14)		
			beta34 = rep(0,A4)
			beta34[as.numeric(slot(RESULTS34,'i'))+1] = slot(RESULTS34,'x')
			beta34 = as.numeric(beta34)
			beta = c(beta13, beta24, beta14, beta34)					
		}else if(ASSUME %in% c('SameHazard')){
			cvfit13 = glmnet::cv.glmnet(x=as.matrix(Cov_long_13[datLONG_sub$trans==1,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==1],datLONG_sub$status[datLONG_sub$trans==1]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==1], alpha = ALPHA) 
			cvfit24 = glmnet::cv.glmnet(x=as.matrix(Cov_long_1424[datLONG_sub$trans==2 | datLONG_sub$trans==3,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==2 | datLONG_sub$trans==3],datLONG_sub$status[datLONG_sub$trans==2 | datLONG_sub$trans==3]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==2 | datLONG_sub$trans==3], alpha = ALPHA) 
			cvfit34 = glmnet::cv.glmnet(x=as.matrix(Cov_long_34[datLONG_sub$trans==4,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==4],datLONG_sub$status[datLONG_sub$trans==4]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==4], alpha = ALPHA) 			
			RESULTS13 = coef(cvfit13, s = "lambda.1se") 
			RESULTS24 = coef(cvfit24, s = "lambda.1se") 
			RESULTS34 = coef(cvfit34, s = "lambda.1se") 			
			beta13 = rep(0,A1)
			beta13[as.numeric(slot(RESULTS13,'i'))+1] = slot(RESULTS13,'x')
			beta13 = as.numeric(beta13)
			beta24 = rep(0,A2)
			beta24[as.numeric(slot(RESULTS24,'i'))+1] = slot(RESULTS24,'x')
			beta24 = as.numeric(beta24)			
			beta14 = beta24	
			beta34 = rep(0,A4)
			beta34[as.numeric(slot(RESULTS34,'i'))+1] = slot(RESULTS34,'x')
			beta34 = as.numeric(beta34)
			beta = c(beta13, beta24, beta14, beta34)				
		}else if(ASSUME %in% c('SameBaseHaz')){			
			cvfit13 = glmnet::cv.glmnet(x=as.matrix(Cov_long_13[datLONG_sub$trans==1,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==1],datLONG_sub$status[datLONG_sub$trans==1]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==1], alpha = ALPHA) 
			cvfit24 = glmnet::cv.glmnet(x=as.matrix(cbind(Cov_long_24[datLONG_sub$trans==2 | datLONG_sub$trans==3,], Cov_long_14[datLONG_sub$trans==2 | datLONG_sub$trans==3,])), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==2 | datLONG_sub$trans==3],datLONG_sub$status[datLONG_sub$trans==2 | datLONG_sub$trans==3]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==2 | datLONG_sub$trans==3], alpha = ALPHA) 
			cvfit34 = glmnet::cv.glmnet(x=as.matrix(Cov_long_34[datLONG_sub$trans==4,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==4],datLONG_sub$status[datLONG_sub$trans==4]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==4], alpha = ALPHA) 			
			RESULTS13 = coef(cvfit13, s = "lambda.1se") 
			RESULTS2414 = coef(cvfit24, s = "lambda.1se") 
			RESULTS34 = coef(cvfit34, s = "lambda.1se") 			
			beta13 = rep(0,A1)
			beta13[as.numeric(slot(RESULTS13,'i'))+1] = slot(RESULTS13,'x')
			beta13 = as.numeric(beta13)
			beta1424 = rep(0,A2 + A3)
			beta1424[as.numeric(slot(RESULTS2414,'i'))+1] = slot(RESULTS2414,'x')
			beta1424 = as.numeric(beta1424)
			beta24 = beta1424[1:A2]		
			beta14 = beta1424[(A2+1:length(beta1424))+1]
			beta34 = rep(0,A4)
			beta34[as.numeric(slot(RESULTS34,'i'))+1] = slot(RESULTS34,'x')
			beta34 = as.numeric(beta34)
			beta = c(beta13, beta24, beta14, beta34)	
		}else if(ASSUME %in% c('ProportionalHazard')){
			cvfit13 = glmnet::cv.glmnet(x=as.matrix(Cov_long_13[datLONG_sub$trans==1,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==1],datLONG_sub$status[datLONG_sub$trans==1]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==1], alpha = ALPHA) 
			cvfit24 = glmnet::cv.glmnet(x=as.matrix(cbind(Cov_long_1424[datLONG_sub$trans==2 | datLONG_sub$trans==3,], 
								as.numeric(datLONG_sub$trans==3)[datLONG_sub$trans==2 | datLONG_sub$trans==3,])), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==2 | datLONG_sub$trans==3],datLONG_sub$status[datLONG_sub$trans==2 | datLONG_sub$trans==3]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==2 | datLONG_sub$trans==3], alpha = ALPHA) 
			cvfit34 = glmnet::cv.glmnet(x=as.matrix(Cov_long_34[datLONG_sub$trans==4,]), 
					y=survival::Surv(datLONG_sub$time[datLONG_sub$trans==4],datLONG_sub$status[datLONG_sub$trans==4]), family= 'cox', 
					weights = datLONG_sub$w[datLONG_sub$trans==4], alpha = ALPHA) 			
			RESULTS13 = coef(cvfit13, s = "lambda.1se") 
			RESULTS2414 = coef(cvfit24, s = "lambda.1se") 
			RESULTS34 = coef(cvfit34, s = "lambda.1se") 			
			beta13 = rep(0,A1)
			beta13[as.numeric(slot(RESULTS13,'i'))+1] = slot(RESULTS13,'x')
			beta13 = as.numeric(beta13)
			beta1424 = rep(0,A2 + A3)
			beta1424[as.numeric(slot(RESULTS2414,'i'))+1] = slot(RESULTS2414,'x')
			beta1424 = as.numeric(beta1424)
			beta24 = beta1424[1:A2]		
			beta14 = beta1424[(A2+1:length(beta1424))]
			beta34 = rep(0,A4)
			beta34[as.numeric(slot(RESULTS34,'i'))+1] = slot(RESULTS34,'x')
			beta34 = as.numeric(beta34)
			beta = c(beta13, beta24, beta14, beta34)		
		}
		
	}#end of PENALTY ifelse


	if(NEEDTOIMPUTE){
		TEMPCOV = c()
		CovImp = ImputeDat[[3]]
		for(k in 1:length(CovImp)){
			TEMPCOV = rbind(TEMPCOV,CovImp[[k]][,TransCov$PNonCure])			
		}
		Out = reshape2::melt(ImputeDat[[4]])$value		
		if(PENALTY == 'None'){
			options(warn=-1)
			fitLogistic = stats::glm(Out~as.matrix(TEMPCOV), family = binomial(link = 'logit'), weights = rep(1/length(CovImp), length(Out)))
			options(warn=1)
			alpha = coef(fitLogistic)			
		}else if(PENALTY %in% c('Ridge', 'Lasso')){
			ALPHA = ifelse(PENALTY == 'Ridge', 0,1)
			cvfit = glmnet::cv.glmnet(x=as.matrix(TEMPCOV), y=Out, family= 'binomial', weights = rep(1/length(CovImp), length(Out)), alpha = ALPHA) 
			RESULTS = coef(cvfit, s = "lambda.1se") 
			alpha = rep(0,1+length(TransCov$PNonCure))
			alpha[as.numeric(slot(RESULTS,'i'))+1] = slot(RESULTS,'x')
			alpha = as.numeric(alpha)
		}
	}else{
		if(PENALTY == 'None'){
			options(warn=-1)		
			fitLogistic = stats::glm(datWIDE$p~as.matrix(Cov[,TransCov$PNonCure]), family = binomial(link = 'logit'))
			options(warn=1)	
			alpha = coef(fitLogistic)
		}else{
			Out = cbind(1-datWIDE$p, datWIDE$p)	
			ALPHA = ifelse(PENALTY == 'Ridge', 0,1)
			cvfit = glmnet::cv.glmnet(x=as.matrix(Cov), y=Out, family= 'binomial', alpha = ALPHA) 
			RESULTS = coef(cvfit, s = "lambda.1se") 
			alpha = rep(0,1+length(TransCov$PNonCure))
			alpha[as.numeric(slot(RESULTS,'i'))+1] = slot(RESULTS,'x')
			alpha = as.numeric(alpha)		
		}
	}
	alpha[is.na(alpha)] = 0
	return(list(beta = as.numeric(beta),alpha = as.numeric(alpha)))
}


