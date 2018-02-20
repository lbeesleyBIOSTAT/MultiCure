

#' STATEOCCUPANCYCOX_IMP
#' @description This function uses RShiny to plot predicted state occupancy probabilities, overall survival probabilities, event-free probabilities, and the transition cumulative hazards for a particular set of covariate values. These probabilities are estimated based on a Multistate cure model fit with COX baseline hazards fit using the MCEM algorithm. The estimated probabilities incorporate only baseline information. This function cannot be applied for prediction when recurrence time is included in the model from recurrence to death. 

#'
#' @param times number of iterations for the EM or MCEM algorithm
#' @param newCov A dataframe with columns corresponding to newCovariates used in the model fit. The rows should correspond to the newCovariate values for which we will make our predictions.
#' @param TransCov a list with elements: Trans13, Trans24, Trans14, Trans34, PNonCure. Each list element is a vector containing the names of the variables in newCov to be used in the model for the corresponding transition. 13 is NonCured -> Recurrence, 24 is Cured -> Death, 14 is NonCured -> Death, 34 is Recurrence -> Death. PNonCure contains the names of the newCovariates for the logistic regression for P(NonCure). 
#' @param beta Estimate from multistate cure model fit
#' @param alpha Estimate from multistate cure model fit
#' @param Basehaz13 Estimate of the baseline hazard for the 1->3 transition. This is estimated using function BaselineHazard_IMP.
#' @param Basehaz24 Estimate of the baseline hazard for the 2->4 transition. This is estimated using function BaselineHazard_IMP.
#' @param Basehaz14 Estimate of the baseline hazard for the 1->4 transition. This is estimated using function BaselineHazard_IMP.
#' @param Basehaz34 Estimate of the baseline hazard for the 3->4 transition. This is estimated using function BaselineHazard_IMP.
#' @param shape Estimate from multistate cure model fit
#'
#' @return NULL
#' @examples
#' attach(SimulateMultiCure(type = "CovariateMissingness"))
#' Cov = data.frame(X1,X2)
#' VARS = names(Cov)
#' TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#' datWIDE = data.frame( Y_R, Y_D, delta_R , delta_D, G)
#' fit = MultiCure(iternum = 100, datWIDE, Cov, ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "cox", IMPNUM = 10, COVIMPUTEFUNCTION = COVIMPUTEFUNCTION, COVIMPUTEINITIALIZE = COVIMPUTEINITIALIZE) ### Note: This will take awhile
#' beta = apply(fit[[3]][,90:100], 1, mean)
#' alpha =  apply(fit[[4]][,90:100], 1, mean)
#' haz = Baselinehazard_IMP(datWIDE, CovImp = fit[[5]],GImp = fit[[6]], YRImp = fit[[7]],deltaRImp = fit[[8]], beta, alpha, TransCov, ASSUME = "SameHazard")	
#' STATEOCCUPANCYCOX_IMP(times = seq(0,100,5), TransCov, newCov = data.frame(X1 = 0, X2 = 0), beta, alpha, Basehaz13 = haz[[1]], Basehaz24 = haz[[2]], Basehaz14 = haz[[3]], Basehaz34 = haz[[4]])  ### Note: This will take awhile
#' @export


 
 
 
STATEOCCUPANCYCOX_IMP <- function(times, TransCov, newCov, beta, alpha, Basehaz13, Basehaz24, Basehaz14, Basehaz34) {

	app = list(	
	ui = shiny::fluidPage(
		shiny::titlePanel("Plotting Options"),
		shiny::sidebarLayout(
			shiny::sidebarPanel(
				shiny::selectInput("plottype", "Plot Type:",
		                c("Transition Hazards" = "HAZ",
		                	"State Occupancy Probabilities" = "S",
		                  "Overall Survival Probability" = "OS",
		                  "Event Free Probability" = "EF")),
		       	shiny::checkboxInput("SavePlot", label = "Save Plots in Working Directory", value = FALSE)  
		    ),
		    shiny::mainPanel(
		       shiny::uiOutput("plots")
		    )
	  	)
	),
	server = function(input, output ) {
		plotInput <- shiny::reactive({
			plot_output_list <- lapply(1:length(newCov[,1]), function(i) {			
			    plotname <- paste("Subject", i, sep="")			
			    shiny::plotOutput(plotname, inline=TRUE)})
			do.call(tagList, plot_output_list)
		})			    			    
		if(length(beta) > length(c(TransCov$Trans14, TransCov$Trans24, TransCov$Trans13, TransCov$Trans34)) & !('INT' %in% TransCov$Trans14)){
	 	 	newCov$INT = rep(1,length(newCov[,1]))
			TransCov$Trans14 = c(TransCov$Trans14, 'INT')
		}
		for(i in 1:length(newCov[,1])){	
			 local({
			 	my_i <- i
    				plotname <- paste("Subject", my_i, sep="")
				output[[plotname]]   <- renderPlot({
					
					#######################################
					### Initialize and Define Functions ###
					#######################################
					
					newCovTEMP = newCov[my_i,]				
					Prob1Save = rep(0,length(times))
					Prob2Save = rep(0,length(times))
					Prob3Save = rep(0,length(times))
					Prob4Save = rep(0,length(times))											
					integrate1_short = function(vr)
					{
						Cumhazard13_temp = exp(XB_beta13)*as.numeric(sapply(vr,Baseline_Hazard, Basehaz13) )
						Cumhazard14_temp = exp(XB_beta14)*as.numeric(sapply(vr,Baseline_Hazard, Basehaz14))
						Cumhazard34_temp = exp(XB_beta34)*as.numeric(sapply(t-vr,Baseline_Hazard, Basehaz34))			
						Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
						Surv3_temp = exp(-Cumhazard34_temp)
						hazard13_temp = exp(XB_beta13)*as.numeric(BasehazFun_13(vr) )
						return(hazard13_temp*Surv1_temp*(1-Surv3_temp))   
					}
					integrate1_shortMOD = function(vr)
					{
						Cumhazard13_temp = exp(XB_beta13)*as.numeric(sapply(vr,Baseline_Hazard, Basehaz13) )
						Cumhazard14_temp = exp(XB_beta14)*as.numeric(sapply(vr,Baseline_Hazard, Basehaz14))
						Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
						hazard13_temp = exp(XB_beta13)*as.numeric(BasehazFun_13(vr) )
						return(hazard13_temp*Surv1_temp)   
					}
					integrate2_short = function(vr)
					{
						Cumhazard13_temp = exp(XB_beta13)*as.numeric(sapply(vr,Baseline_Hazard, Basehaz13) )
						Cumhazard14_temp = exp(XB_beta14)*as.numeric(sapply(vr,Baseline_Hazard, Basehaz14))
						Cumhazard34_temp = exp(XB_beta34)*as.numeric(sapply(t-vr,Baseline_Hazard, Basehaz34))			
						Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
						Surv3_temp = exp(-Cumhazard34_temp)
						hazard13_temp = exp(XB_beta13)*as.numeric(BasehazFun_13(vr))
						return(hazard13_temp*Surv1_temp*Surv3_temp)   
					}
					integrate4_short = function(vd)
					{	
						Cumhazard13_temp = exp(XB_beta13)*as.numeric(sapply(vd, Baseline_Hazard,Basehaz13 ))
						Cumhazard14_temp = exp(XB_beta14)*as.numeric(sapply(vd, Baseline_Hazard,Basehaz14 ))
						Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
						hazard14_temp = exp(XB_beta14)*as.numeric(BasehazFun_14(vd))
						return(hazard14_temp*Surv1_temp)   
					}
					
					XB_alpha = as.numeric(alpha %*% t(cbind(1, newCovTEMP[,TransCov$PNonCure])))
					prob_Noncure = exp(XB_alpha)/(1+exp(XB_alpha))
					
					A1 = length(TransCov$Trans13)
					A2 = length(TransCov$Trans24)
					A3 = length(TransCov$Trans14)
					A4 = length(TransCov$Trans34)
					TRANS = c(rep(1,A1), rep(2,A2), rep(3,A3), rep(4,A4))
					XB_beta13 = as.numeric(beta[TRANS==1] %*% t(cbind(newCovTEMP[,TransCov$Trans13])))	
					XB_beta24 = as.numeric(beta[TRANS==2] %*% t(cbind(newCovTEMP[,TransCov$Trans24])))	
					XB_beta14 = as.numeric(beta[TRANS==3] %*% t(cbind(newCovTEMP[,TransCov$Trans14])))	
					XB_beta34 = as.numeric(beta[TRANS==4] %*% t(cbind(newCovTEMP[,TransCov$Trans34])))		
					
					BasehazFun_13 = stepfun(x= Basehaz13[,2], y = c(Basehaz13[,3],0), right = F)
					BasehazFun_24 = stepfun(x= Basehaz24[,2], y = c(Basehaz24[,3],0), right = F)
					BasehazFun_14 = stepfun(x= Basehaz14[,2], y = c(Basehaz14[,3],0), right = F)
					BasehazFun_34 = stepfun(x= Basehaz34[,2], y = c(Basehaz34[,3],0), right = F)												


					#######################################################
					### Estimate and Plot State Occupancy Probabilities ###
					#######################################################
					
  					cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
					if(input$plottype=='S'){
						changePoints1 = sort(unique(c(0,Basehaz13[,1], Basehaz24[,1], Basehaz14[,1], max(times))))
						changePoints1 = changePoints1[changePoints1 <= max(times)]							
						Group1_tempMOD = Group4_temp = rep(0, length(changePoints1)-1)
						INT1MOD = integrate1_shortMOD(changePoints1)
						INT1MOD_Minus = integrate1_shortMOD(pmax(changePoints1-1e-10,0))
						INT4 = integrate4_short(changePoints1)
						INT4_Minus = integrate4_short(pmax(changePoints1-1e-10,0))
						if(length(changePoints1)>1){
							for(s in 2:length(changePoints1))
							{
								b = (log(INT1MOD[s-1]+1e-23) - log(INT1MOD_Minus[s]+1e-23))/(changePoints1[s] - changePoints1[s-1])
								a = INT1MOD[s-1]*exp(b* changePoints1[s-1])
								if(a >= 1e-10 & !is.infinite(a) & !is.infinite(b)){
									Group1_tempMOD[s-1] = - (a/b)*(exp(-b* changePoints1[s])- exp(-b* changePoints1[s-1]))
								}									
								b = (log(INT4[s-1]+1e-23) - log(INT4_Minus[s]+1e-23))/(changePoints1[s] - changePoints1[s-1])
								a = INT4[s-1]*exp(b* changePoints1[s-1])
								if(a>= 1e-10 & !is.infinite(a) & !is.infinite(b)){
									Group4_temp[s-1] = - (a/b)*(exp(-b* changePoints1[s])- exp(-b* changePoints1[s-1]))
								}
							}
						}
						for(k in 1:length(times))
						{	
							t = times[k]	
							changePoints = sort(unique(c(0,Basehaz13[,1], Basehaz24[,1], Basehaz14[,1],max(t - Basehaz34[,1],0),t)))
							changePoints = changePoints[changePoints<=t]
							INT2 = integrate2_short(changePoints)
							INT2_Minus = integrate2_short(pmax(changePoints-1e-10,0))
							
							Group1 = Group2 = Group4 = 0
							if(length(changePoints)>1){
								for(s in 2:length(changePoints))
								{
									b = (log(INT2[s-1]+1e-23) - log(INT2_Minus[s]+1e-23))/(changePoints[s] - changePoints[s-1])
									a = INT2[s-1]*exp(b*changePoints[s-1])
									if(a >= 1e-10 & !is.infinite(a) & !is.infinite(b)){Group2 = Group2 - (a/b)*(exp(-b*changePoints[s])- exp(-b*changePoints[s-1]))}
								}
							}
							
							K = max(which(changePoints1<=t))
							if(changePoints1[K]==0){
								Group1 = 0
								Group4 = 0
							}else if(changePoints1[K]==t & t > 0){
								Group1 = cumsum(Group1_tempMOD)[K-1] - Group2
								Group4 = cumsum(Group4_temp)[K-1]
							}else{
								Group1 = cumsum(Group1_tempMOD)[K-1] - Group2 + cubature::adaptIntegrate(Vectorize(integrate1_shortMOD), lowerLimit = changePoints1[K], upperLimit = t)$integral 
								Group4 = cumsum(Group4_temp)[K-1] + cubature::adaptIntegrate(Vectorize(integrate4_short), lowerLimit = changePoints1[K], upperLimit = t)$integral 
							}
							
							Group3 = exp(-Baseline_Hazard(t, Basehaz13 )*exp(XB_beta13))*exp(-Baseline_Hazard(t, Basehaz14 )*exp(XB_beta14))
						
							Prob1Save[k] = Group1*prob_Noncure
							Prob2Save[k] = Group2*prob_Noncure
							Prob3Save[k] = Group3*prob_Noncure + (exp(-Baseline_Hazard(t, Basehaz24 )*exp(XB_beta24))* (1-prob_Noncure))
							Prob4Save[k] = Group4*prob_Noncure +  ((1-exp(-Baseline_Hazard(t, Basehaz24 )*exp(XB_beta24)))* (1-prob_Noncure))
						}						
						SUBSET = c(1:length(times))	
						#par(mar = c(5,5,5,6), xpd = TRUE)
						plot(c(),c(), main = paste0('State Occupancy Probabilities \n', plotname), ylim = c(0,1), xlim = c(0,max(times)), xlab = 'Time from Baseline', 
							ylab = 'Proportion of Subjects in Group')
						
						lower = rep(0,length(SUBSET))
						upper = Prob2Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[8])		
						
						lower = upper
						upper = lower + Prob1Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])
							
						lower = upper
						upper = lower + Prob4Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[4])
						
						lower = upper
						upper = lower + Prob3Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
			
						legend(x='topleft', legend = c('Recurred and Alive', 'Recurred and Died', 'Died without Recurrence', 'Alive without Recurrence'), 
								fill = cbPalette[c(8,5,4,3)],
								bty = 'n', cex = 0.7)
					}else if(input$plottype == 'OS'){
						SUBSET = c(1:length(times))									
						for(k in 1:length(times))
						{	
							t = times[k]	
							changePoints = sort(unique(c(0,Basehaz13[,1], Basehaz24[,1], Basehaz14[,1],max(t - Basehaz34[,1],0),t)))
							changePoints = changePoints[changePoints<=t]
							INT2 = integrate2_short(changePoints)
							INT2_Minus = integrate2_short(pmax(changePoints-1e-10,0))									
							Group2 = 0
							if(length(changePoints)>1){
								for(s in 2:length(changePoints))
								{
									b = (log(INT2[s-1]+1e-23) - log(INT2_Minus[s]+1e-23))/(changePoints[s] - changePoints[s-1])
									a = INT2[s-1]*exp(b*changePoints[s-1])
									if(a >= 1e-10 & !is.infinite(a) & !is.infinite(b)){Group2 = Group2 - (a/b)*(exp(-b*changePoints[s])- exp(-b*changePoints[s-1]))}
								}
							}
							Group3 = exp(-Baseline_Hazard(t, Basehaz13 )*exp(XB_beta13))*exp(-Baseline_Hazard(t, Basehaz14 )*exp(XB_beta14))								
							Prob2Save[k] = Group2*prob_Noncure
							Prob3Save[k] = Group3*prob_Noncure + (exp(-Baseline_Hazard(t, Basehaz24 )*exp(XB_beta24))* (1-prob_Noncure))
						}
		
						plot(c(),c(), main = paste0('Overall Survival Probability \n', plotname), ylim = c(0,1), xlim = c(0,max(times)), xlab = 'Time from Baseline', 
							ylab = 'Survival Probability')
						
						lower = rep(0,length(SUBSET))
						upper = Prob2Save[SUBSET] + Prob3Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])		
						
						lower = upper
						upper = rep(1, length(lower))
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
			
						legend(x='bottomleft', legend = c('Alive', 'Died'), 
								fill = cbPalette[c(5,3)],
								bty = 'n', cex = 0.7)				
					}else if(input$plottype == 'EF'){						
						for(k in 1:length(times))
						{	
							t = times[k]	
							Group3 = exp(-Baseline_Hazard(t, Basehaz13 )*exp(XB_beta13))*exp(-Baseline_Hazard(t, Basehaz14 )*exp(XB_beta14))								
							Prob3Save[k] = Group3*prob_Noncure + (exp(-Baseline_Hazard(t, Basehaz24 )*exp(XB_beta24))* (1-prob_Noncure))
						}						
						SUBSET = c(1:length(times))	
						#par(mar = c(5,5,5,6), xpd = TRUE)
						plot(c(),c(), main = paste0('Event-Free Probability \n', plotname), ylim = c(0,1), xlim = c(0,max(times)), xlab = 'Time from Baseline', 
							ylab = 'Event-Free Probability')
						
						lower = rep(0,length(SUBSET))
						upper = Prob3Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])		
						
						lower = upper
						upper = rep(1, length(lower)) 
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
			
						legend(x='bottomleft', legend = c('Alive without Recurrence', 'Died and/or Recurred'), 
								fill = cbPalette[c(5,3)],
								bty = 'n', cex = 0.7)			
					}else{
						
						H13 = sapply(times , Baseline_Hazard, Basehaz13)*exp(XB_beta13)
						H14 = sapply(times , Baseline_Hazard, Basehaz14)*exp(XB_beta14)
						H24 = sapply(times , Baseline_Hazard, Basehaz24)*exp(XB_beta24)
						H34 = sapply(times , Baseline_Hazard, Basehaz34)*exp(XB_beta34)

						plot(c(),c(), main = paste0('Transition Event Hazards \n', plotname), ylim = c(0,max(c(H13, H14, H24, H34))), xlim = c(0,max(times)), xlab = 'Time (with Clock Reset)', 
							ylab = 'Cumul. Hazard')
						lines(times, H13, col = cbPalette[8], lwd = 2)
						lines(times,H24, col = cbPalette[5], lwd = 2)
						lines(times,H14, col = cbPalette[4], lwd = 2)
						lines(times,H34, col = cbPalette[3], lwd = 2)
			
						legend(x='topleft', legend = c('Not Cured -> Recurrence', 'Cured -> Death', 'Not Cured -> Death', 'Recurrence -> Death'), 
								fill = cbPalette[c(8,5,4,3)],
								bty = 'n', cex = 0.7)	
					}
					if(input$SavePlot == TRUE){
					grDevices::dev.print(pdf, paste(plotname, '.pdf', sep = ''))
					}					
				}, height = 280, width = 250)#end of renderplot
			})#end of local			
		}#end of for loop
		
		output$plots <- shiny::renderUI({print(plotInput())})	

	}
	)
	shiny::runApp(app)
}


			