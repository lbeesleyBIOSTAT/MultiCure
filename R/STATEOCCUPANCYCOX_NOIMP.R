

#' STATEOCCUPANCYCOX_NOIMP
#' @description This function uses RShiny to plot predicted state occupancy probabilities, overall survival probabilities, event-free probabilities, and the transition cumulative hazards for a particular set of covariate values. These probabilities are estimated based on a Multistate cure model fit with COX baseline hazards fit using the EM algorithm. The estimated probabilities can incorporate only baseline information OR can incorporate baseline information and some followup information after baseline. This function cannot be applied for prediction when recurrence time is included in the model from recurrence to death. 

#'
#' @param times number of iterations for the EM or MCEM algorithm
#' @param newCov A dataframe with columns corresponding to newCovariates used in the model fit. The rows should correspond to the newCovariate values for which we will make our predictions.
#' @param TransCov a list with elements: Trans13, Trans24, Trans14, Trans34, PNonCure. Each list element is a vector containing the names of the variables in newCov to be used in the model for the corresponding transition. 13 is NonCured -> Recurrence, 24 is Cured -> Death, 14 is NonCured -> Death, 34 is Recurrence -> Death. PNonCure contains the names of the newCovariates for the logistic regression for P(NonCure). 
#' @param beta Estimate from multistate cure model fit
#' @param alpha Estimate from multistate cure model fit
#' @param Haz_13 Estimate of the cumulative baseline hazard for the 1->3 transition. This is estimated using function BaselineHazard_NOIMP.
#' @param Haz_24 Estimate of the cumulative baseline hazard for the 2->4 transition. This is estimated using function BaselineHazard_NOIMP.
#' @param Haz_14 Estimate of the cumulative baseline hazard for the 1->4 transition. This is estimated using function BaselineHazard_NOIMP.
#' @param Haz_34 Estimate of the cumulative baseline hazard for the 3->4 transition. This is estimated using function BaselineHazard_NOIMP.
#' @param shape Estimate from multistate cure model fit
#'
#' @return NULL
#' @examples
#' attach(SimulateMultiCure(type = "NoMissingness"))
#' Cov = data.frame(X1,X2)
#' VARS = names(Cov)
#' TransCov = list(Trans13 = VARS, Trans24 = VARS, Trans14 = VARS, Trans34 = VARS, PNonCure = VARS)
#' datWIDE = data.frame( Y_R, Y_D, delta_R , delta_D, G)
#' fit = MultiCure(iternum = 100, datWIDE, Cov, ASSUME = "SameHazard", TransCov = TransCov, BASELINE = "cox")
#' Haz = BaselineHazard_NOIMP(datWIDE, Cov, beta = fit[[1]], alpha = fit[[2]], TransCov, ASSUME = "SameHazard", p = fit[[5]][,100])	
#' STATEOCCUPANCYCOX_NOIMP(times = seq(0,100,1), TransCov, newCov = data.frame(X1 = 0, X2 = 0), beta = fit[[1]], alpha = fit[[2]], Haz_13 = Haz[[1]], Haz_24 = Haz[[2]], Haz_14 = Haz[[3]], Haz_34 = Haz[[4]])  
#' @export


 
 
 
STATEOCCUPANCYCOX_NOIMP <- function(times, TransCov, newCov, beta, alpha, Haz_13, Haz_24, Haz_14, Haz_34) {

	app = list(
	ui = shiny::fluidPage(
		shiny::titlePanel("Plotting Options"),
		shiny::sidebarLayout(
			shiny::sidebarPanel(
				shiny::selectInput("plottype", "Plot Type:",
	                c("State Occupancy Probabilities" = "S",
	                  "Overall Survival Probability" = "OS",
	                  "Event Free Probability" = "EF",
	                  "Transition Hazards" = "HAZ")),
			   	shiny::conditionalPanel(condition="input. plottype !='HAZ'",
	             		shiny::sliderInput('CurTime', "Patient known to be alive at time: ", 0, min = 0, max = max(times)) ,
					shiny::checkboxInput("RecurEvent", label = "Patient had observed recurrence", value = FALSE) ),
			   	shiny::conditionalPanel(condition="input.RecurEvent==true",
	                  shiny::sliderInput('RecurTime', "When did the recurrence occur?: ", min = 0, max=max(times), value = 0) ),
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
					newCovTEMP = newCov[my_i,]				
					Prob1Save = rep(0,length(times))
					Prob2Save = rep(0,length(times))
					Prob3Save = rep(0,length(times))
					Prob4Save = rep(0,length(times))											
					integrate1_short = function(vr)
					{
						
						Cumhazard13_temp = exp(XB_beta13)*as.numeric(Haz_13(vr))
						Cumhazard14_temp = exp(XB_beta14)*as.numeric(Haz_14(vr))
						Cumhazard34_temp = exp(XB_beta34)*as.numeric(Haz_34(t-vr) )				
						Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
						Surv3_temp = exp(-Cumhazard34_temp)
						hazard13_temp = exp(XB_beta13)*as.numeric((Haz_13(vr+0.0001)-  Haz_13(vr-0.0001)    ))
						return(hazard13_temp*Surv1_temp*(1-Surv3_temp))   
					}
					integrate2_short = function(vr)
					{
						Cumhazard13_temp = exp(XB_beta13)*as.numeric(Haz_13(vr))
						Cumhazard14_temp = exp(XB_beta14)*as.numeric(Haz_14(vr))
						Cumhazard34_temp = exp(XB_beta34)*as.numeric(Haz_34(t-vr) )				
						Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
						Surv3_temp = exp(-Cumhazard34_temp)
						hazard13_temp = exp(XB_beta13)*as.numeric((Haz_13(vr+0.0001)-  Haz_13(vr-0.0001)    ))
						return(hazard13_temp*Surv1_temp*Surv3_temp)   
					}
					integrate4_short = function(vd)
					{	
						Cumhazard13_temp = exp(XB_beta13)*as.numeric(Haz_13(vd))
						Cumhazard14_temp = exp(XB_beta14)*as.numeric(Haz_14(vd))
						Surv1_temp = exp(-Cumhazard13_temp-Cumhazard14_temp)
						hazard14_temp = exp(XB_beta14)*as.numeric((Haz_14(vd+0.0001)-  Haz_14(vd-0.0001)    ))
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
												
					if(input$RecurTime > input$CurTime & input$RecurEvent == TRUE){
						stop('Recurrence time must be at or before current time')
					}
					recurevents_long = knots(Haz_13)
					deathevents_long = knots(Haz_14)
					if(input$CurTime == 0 & input$RecurEvent==FALSE){
						for(k in 1:length(times))
						{
							t = times[k]	
							recurevents = recurevents_long[recurevents_long <=t]
							deathevents = deathevents_long[deathevents_long <=t]
							Group1 = ifelse(length(recurevents)!=0,sum(integrate1_short(recurevents )), 0)
							Group2 =  ifelse(length(recurevents)!=0,sum(integrate2_short(recurevents )), 0)
							Group3 = exp(-Haz_13(t)*exp(XB_beta13))*exp(-Haz_14(t)*exp(XB_beta14))
							Group4 = ifelse(length(deathevents)!=0,sum(integrate4_short(deathevents )),0)
							Prob1Save[k] = Group1*prob_Noncure
							Prob2Save[k] = Group2*prob_Noncure
							Prob3Save[k] = Group3*prob_Noncure + (exp(-Haz_24(t)*exp(XB_beta24))* (1-prob_Noncure))
							Prob4Save[k] = Group4*prob_Noncure +  ((1-exp(-Haz_24(t)*exp(XB_beta24)))* (1-prob_Noncure))
						}
					}else{
						if(input$RecurEvent == FALSE){ #no observed recurrence
							for(k in which(times>=input$CurTime))
							{
								t = times[k]		
								recurevents = recurevents_long[recurevents_long <=t]
								deathevents = deathevents_long[deathevents_long <=t]
								Group1 = ifelse(length(recurevents)!=0,sum(integrate1_short(recurevents )), 0)
								Group2 =  ifelse(length(recurevents)!=0,sum(integrate2_short(recurevents )), 0)
								Group3 = exp(-Haz_13(t)*exp(XB_beta13))*exp(-Haz_14(t)*exp(XB_beta14))
								Group4 = ifelse(length(deathevents)!=0,sum(integrate4_short(deathevents )),0)
								recurevents = recurevents_long[recurevents_long <=input$CurTime ]
								deathevents = deathevents_long[deathevents_long <=input$CurTime ]
								Group1STAR = ifelse(length(recurevents)!=0,sum(integrate1_short(recurevents )), 0)
								Group2STAR = ifelse(length(recurevents)!=0,sum(integrate2_short(recurevents )), 0)
								Group3STAR = exp(-Haz_13(input$CurTime)*exp(XB_beta13))*exp(-Haz_14(input$CurTime)*exp(XB_beta14))
								Group4STAR = ifelse(length(deathevents)!=0,sum(integrate4_short(deathevents )),0)
								prob_NoncureSTAR = Group3STAR*prob_Noncure/(Group3STAR*prob_Noncure + exp(-Haz_24(input$CurTime)*exp(XB_beta24))*(1-prob_Noncure))
								Prob1Save[k] = ((Group1-Group1STAR)/Group3STAR)* prob_NoncureSTAR
								Prob2Save[k] = ((Group2-Group2STAR)/Group3STAR)* prob_NoncureSTAR
								Prob3Save[k] = (exp(-Haz_24(t)*exp(XB_beta24))/exp(-Haz_24(input$CurTime)*exp(XB_beta24)))*(1-prob_NoncureSTAR)+
												(Group3/Group3STAR)*prob_NoncureSTAR
								Prob4Save[k] = (1-(exp(-Haz_24(t)*exp(XB_beta24))/exp(-Haz_24(input$CurTime)*exp(XB_beta24))))*(1-prob_NoncureSTAR)+
												((Group4-Group4STAR)/Group3STAR)*prob_NoncureSTAR
							}
						}else{ #observed recurrence
							for(k in which(times>=input$CurTime))
							{
								t = times[k]		
								Prob1Save[k] =  1-(exp(-Haz_34(t-input$RecurTime)*exp(XB_beta34))/ exp(-Haz_34(input$CurTime-input$RecurTime)*exp(XB_beta34)))
								Prob2Save[k] =  exp(-Haz_34(t-input$RecurTime)*exp(XB_beta34))/ exp(-Haz_34(input$CurTime-input$RecurTime)*exp(XB_beta34))
								Prob3Save[k] = 0
								Prob4Save[k] = 0
							}					
						}
						
					}
				#2, recurred and alive; 1, recurred and died; 4, died without recurrence; 3, alive without recurrence
		
  					cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
					if(input$plottype=='S'){
						SUBSET = which(times>=input$CurTime)
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
						if(input$CurTime>0){abline(v=min(times[SUBSET]), lwd = 3)}
			
						legend(x='topleft', legend = c('Recurred and Alive', 'Recurred and Died', 'Died without Recurrence', 'Alive without Recurrence'), 
								fill = cbPalette[c(8,5,4,3)],
								bty = 'n', cex = 0.7)
					}else if(input$plottype == 'OS'){
									SUBSET = which(times>=input$CurTime)
						#par(mar = c(5,5,5,6), xpd = TRUE)
						plot(c(),c(), main = paste0('Overall Survival Probability \n', plotname), ylim = c(0,1), xlim = c(0,max(times)), xlab = 'Time from Baseline', 
							ylab = 'Survival Probability')
						
						lower = rep(0,length(SUBSET))
						upper = Prob2Save[SUBSET] + Prob3Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])		
						
						lower = upper
						upper = lower + Prob1Save[SUBSET] + Prob4Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
						if(input$CurTime>0){abline(v=min(times[SUBSET]), lwd = 3)}
			
						legend(x='bottomleft', legend = c('Alive', 'Died'), 
								fill = cbPalette[c(5,3)],
								bty = 'n', cex = 0.7)				
					}else if(input$plottype == 'EF'){
						SUBSET = which(times>=input$CurTime)
						#par(mar = c(5,5,5,6), xpd = TRUE)
						plot(c(),c(), main = paste0('Event-Free Probability \n', plotname), ylim = c(0,1), xlim = c(0,max(times)), xlab = 'Time from Baseline', 
							ylab = 'Event-Free Probability')
						
						lower = rep(0,length(SUBSET))
						upper = Prob3Save[SUBSET]
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])		
						
						lower = upper
						upper = lower + Prob1Save[SUBSET] + Prob4Save[SUBSET] + Prob2Save[SUBSET] 
						lines(times[SUBSET], upper  , lwd = 3)
						polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
						if(input$CurTime>0){abline(v=min(times[SUBSET]), lwd = 3)}
			
						legend(x='bottomleft', legend = c('Alive without Recurrence', 'Died and/or Recurred'), 
								fill = cbPalette[c(5,3)],
								bty = 'n', cex = 0.7)			
					}else{
						H13 = sapply(times ,Haz_13)*exp(XB_beta13)
						H14 = sapply(times ,Haz_14)*exp(XB_beta14)
						H24 = sapply(times ,Haz_24)*exp(XB_beta24)
						H34 = sapply(times ,Haz_34)*exp(XB_beta34)			
								
						plot(c(),c(), main = paste0('Cumulative Hazard \n', plotname), ylim = c(0,max(c(H13, H14, H24, H34))), xlim = c(0,max(times)), xlab = 'Time (with Clock Reset)', 
							ylab = 'Cumulative Hazard')
						lines(times,H13, col = cbPalette[8], lwd = 2)
						lines(times,H24, col = cbPalette[5], lwd = 2)
						lines(times,H14, col = cbPalette[4], lwd = 2)
						lines(times,H34, col = cbPalette[3], lwd = 2)
			
						legend(x='bottomright', legend = c('Not Cured -> Recurrence', 'Cured -> Death', 'Not Cured -> Death', 'Recurrence -> Death'), 
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


			