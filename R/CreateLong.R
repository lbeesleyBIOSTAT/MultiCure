
### These functions convert the wide outcome dataset to a stacked long outcome dataset (augmented data). We reformat the wide dataset into a long outcome dataset for each of the IMPNUM imputed datasets and stack these on top of eachother. Each row corresponds to time at risk for a particular transition for a particular individual in a particular imputed dataset. 

#' @export

CreateLong_MC = function(datWIDE, ImputeDat){
	CovImp = ImputeDat[[3]]
	GImp = ImputeDat[[4]]
	YRImp = ImputeDat[[5]]
	deltaRImp = ImputeDat[[6]]	
	NOBS = length(datWIDE[,1])
	N = 4*NOBS	
	datLONGImp = c()
	IMPNUM = length(ImputeDat[[3]]) 
	for(i in 1:IMPNUM)
	{
		datWIDE_temp = data.frame(id = rep(1:NOBS), 
			Tstart1 = rep(NA,NOBS), Tstop1 = rep(NA,NOBS), time1 = rep(NA,NOBS), status1 = rep(NA,NOBS), w1 = rep(NA,NOBS),
			Tstart2 = rep(NA,NOBS), Tstop2 = rep(NA,NOBS), time2 = rep(NA,NOBS), status2 = rep(NA,NOBS), w2 = rep(NA,NOBS),
			Tstart3 = rep(NA,NOBS), Tstop3 = rep(NA,NOBS), time3 = rep(NA,NOBS), status3 = rep(NA,NOBS), w3 = rep(NA,NOBS),
			Tstart4 = rep(NA,NOBS), Tstop4 = rep(NA,NOBS), time4 = rep(NA,NOBS), status4 = rep(NA,NOBS), w4 = rep(NA,NOBS),
			data.frame(CovImp[[i]]))
		#13
		datWIDE_temp$Tstart1 = rep(0,NOBS)
		datWIDE_temp$Tstop1 = YRImp[,i]
		datWIDE_temp$time1 = datWIDE_temp$Tstop1-datWIDE_temp$Tstart1
		datWIDE_temp$status1 = deltaRImp[,i]
		datWIDE_temp$w1 = (GImp[,i])/IMPNUM
		#24
		datWIDE_temp$Tstart2 = ifelse(deltaRImp[,i] == 1,rep(NA,NOBS), rep(0,NOBS))
		datWIDE_temp$Tstop2 = ifelse(deltaRImp[,i] == 1,rep(NA,NOBS), YRImp[,i])
		datWIDE_temp$status2 = ifelse(deltaRImp[,i] == 1,rep(NA,NOBS), datWIDE$delta_D)
		datWIDE_temp$time2 = datWIDE_temp$Tstop2-datWIDE_temp$Tstart2
		datWIDE_temp$w2 = (1-GImp[,i])/IMPNUM		
		#14
		datWIDE_temp$Tstart3 = rep(0,NOBS)
		datWIDE_temp$Tstop3 = YRImp[,i]
		datWIDE_temp$status3 = ifelse(deltaRImp[,i] == 1, rep(0,NOBS), datWIDE$delta_D)
		datWIDE_temp$time3 = datWIDE_temp$Tstop3-datWIDE_temp$Tstart3
		datWIDE_temp$w3 = (GImp[,i])/IMPNUM
		#34
		datWIDE_temp$Tstart4 = ifelse(deltaRImp[,i] == 0,rep(NA,NOBS), YRImp[,i])
		datWIDE_temp$Tstop4 = ifelse(deltaRImp[,i] == 0,rep(NA,NOBS), datWIDE$Y_D)
		datWIDE_temp$status4 = ifelse(deltaRImp[,i] == 0,rep(NA,NOBS), datWIDE$delta_D)
		datWIDE_temp$time4 = datWIDE_temp$Tstop4-datWIDE_temp$Tstart4
		datWIDE_temp$w4 = (deltaRImp[,i])/IMPNUM		
		datLONG = stats::reshape(datWIDE_temp, varying = list( c("Tstart1", "Tstart2", "Tstart3", "Tstart4"), 
						 c("Tstop1", "Tstop2", "Tstop3", "Tstop4"),
						 c("time1", "time2", "time3", "time4"),
						 c("status1", "status2", "status3", "status4"),
						 c("w1", "w2", "w3", "w4")), 
						 v.names = c('Tstart', 'Tstop', 'time', 'status', 'w'),
						idvar = "id", timevar = "trans", times = c(1,2,3,4), direction = "long", sep = "")
		datLONG$from = rep(NA,N)
		datLONG$to = rep(NA,N)		
		datLONG$from[datLONG$trans == 1] = 1
		datLONG$to[datLONG$trans == 1] = 3
		datLONG$from[datLONG$trans == 2] = 2
		datLONG$to[datLONG$trans == 2] = 4
		datLONG$from[datLONG$trans == 3] = 1
		datLONG$to[datLONG$trans == 3] = 4
		datLONG$from[datLONG$trans == 4] = 3
		datLONG$to[datLONG$trans == 4] = 4		
		datLONGImp= rbind(datLONGImp, datLONG)
	}
	datLONGImp = datLONGImp[complete.cases(datLONGImp),]
	return(datLONGImp)
}






### These functions convert the wide outcome dataset to a long outcome dataset. Each row corresponds to time at risk for a particular transition for a particular individual. 

#' @export

CreateLong = function(datWIDE, Cov){
	NOBS = length(datWIDE[,1])
	N = 4*NOBS	

	datWIDE_temp = data.frame(id = rep(1:NOBS), 
		Tstart1 = rep(NA,NOBS), Tstop1 = rep(NA,NOBS), time1 = rep(NA,NOBS), status1 = rep(NA,NOBS), w1 = rep(NA,NOBS),
		Tstart2 = rep(NA,NOBS), Tstop2 = rep(NA,NOBS), time2 = rep(NA,NOBS), status2 = rep(NA,NOBS), w2 = rep(NA,NOBS),
		Tstart3 = rep(NA,NOBS), Tstop3 = rep(NA,NOBS), time3 = rep(NA,NOBS), status3 = rep(NA,NOBS), w3 = rep(NA,NOBS),
		Tstart4 = rep(NA,NOBS), Tstop4 = rep(NA,NOBS), time4 = rep(NA,NOBS), status4 = rep(NA,NOBS), w4 = rep(NA,NOBS),
		data.frame(Cov))
	#13
	datWIDE_temp$Tstart1 = rep(0,NOBS)
	datWIDE_temp$Tstop1 = datWIDE$Y_R
	datWIDE_temp$time1 = datWIDE_temp$Tstop1-datWIDE_temp$Tstart1
	datWIDE_temp$status1 = datWIDE$delta_R
	datWIDE_temp$w1 = datWIDE$p
	#24
	datWIDE_temp$Tstart2 = ifelse(datWIDE$delta_R == 1,rep(NA,NOBS), rep(0,NOBS))
	datWIDE_temp$Tstop2 = ifelse(datWIDE$delta_R == 1,rep(NA,NOBS), datWIDE$Y_R)
	datWIDE_temp$status2 = ifelse(datWIDE$delta_R == 1,rep(NA,NOBS), datWIDE$delta_D)
	datWIDE_temp$time2 = datWIDE_temp$Tstop2-datWIDE_temp$Tstart2
	datWIDE_temp$w2 = 1-datWIDE$p		
	#14
	datWIDE_temp$Tstart3 = rep(0,NOBS)
	datWIDE_temp$Tstop3 = datWIDE$Y_R
	datWIDE_temp$status3 = ifelse(datWIDE$delta_R == 1, rep(0,NOBS), datWIDE$delta_D)
	datWIDE_temp$time3 = datWIDE_temp$Tstop3-datWIDE_temp$Tstart3
	datWIDE_temp$w3 = datWIDE$p
	#34
	datWIDE_temp$Tstart4 = ifelse(datWIDE$delta_R == 0,rep(NA,NOBS),datWIDE$Y_R)
	datWIDE_temp$Tstop4 = ifelse(datWIDE$delta_R == 0,rep(NA,NOBS), datWIDE$Y_D)
	datWIDE_temp$status4 = ifelse(datWIDE$delta_R == 0,rep(NA,NOBS), datWIDE$delta_D)
	datWIDE_temp$time4 = datWIDE_temp$Tstop4-datWIDE_temp$Tstart4
	datWIDE_temp$w4 = datWIDE$delta_R	
	datLONG = stats::reshape(datWIDE_temp, varying = list( c("Tstart1", "Tstart2", "Tstart3", "Tstart4"), 
					 c("Tstop1", "Tstop2", "Tstop3", "Tstop4"),
					 c("time1", "time2", "time3", "time4"),
					 c("status1", "status2", "status3", "status4"),
					 c("w1", "w2", "w3", "w4")), 
					 v.names = c('Tstart', 'Tstop', 'time', 'status', 'w'),
					idvar = "id", timevar = "trans", times = c(1,2,3,4), direction = "long", sep = "")
	datLONG$from = rep(NA,N)
	datLONG$to = rep(NA,N)		
	datLONG$from[datLONG$trans == 1] = 1
	datLONG$to[datLONG$trans == 1] = 3
	datLONG$from[datLONG$trans == 2] = 2
	datLONG$to[datLONG$trans == 2] = 4
	datLONG$from[datLONG$trans == 3] = 1
	datLONG$to[datLONG$trans == 3] = 4
	datLONG$from[datLONG$trans == 4] = 3
	datLONG$to[datLONG$trans == 4] = 4		

	datLONG = datLONG[complete.cases(datLONG),]
	return(datLONG)
}




