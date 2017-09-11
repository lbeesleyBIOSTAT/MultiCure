
#' VISUALIZEDATA
#' @description This function uses RShiny to plot Kaplan-Meier plots for time to recurrence and time to death. This function can also be used to visualize the occurrence of unequal follow-up and when the recurrence event/censorings and death event/censorings occur.
#'
#' @param datWIDE A data frame with the following columns: 
#' \itemize{
#' \item Y_R, the recurrence event/censoring time
#' \item delta_R, the recurrence event/censoring indicator
#'\item Y_D, the death event/censoring time 
#' \item delta_D, the death event/censoring indicator
#'}
#'
#' @return NULL
#'
#' @author Lauren J Beesley, \email{lbeesley@umich.edu}
#' @export


VISUALIZEDATA <- function(datWIDE) {

	app = list(
			ui = shiny::fluidPage(
		
		  # Application title
		  shiny::titlePanel("Plotting Options"),
		
		  # Sidebar with a slider input for the number of bins
		  shiny::sidebarLayout(
		    shiny::sidebarPanel(width = 4,
		      shiny::textInput("title", "Title:", "Follow-up for Recurrence and Death"),
		     	 shiny::checkboxInput("RecurEvent", label = "Plot Recurrence Events", value = FALSE),
		     	 shiny::checkboxInput("RecurCens", label = "Plot Recurrence Censorings", value = FALSE),
		      	shiny::checkboxInput("DeathEvent", label = "Plot Death Events", value = FALSE),
		      	      	shiny::checkboxInput("DeathCens", label = "Plot Death Censorings", value = FALSE),
		
		      # selectInput("RecurEvent", "Plot Recurrence Events:", 
		                  # choices = c("No", "Yes")),
		      
		       # selectInput("RecurCens", "Plot Recurrence Censorings:", 
		                  # choices = c("No", "Yes")),
		       # selectInput("DeathEvent", "Plot Death Events:", 
		                  # choices = c("No", "Yes")),
		                  
		        # selectInput("DeathCens", "Plot Death Censorings:", 
		                  # choices = c("No", "Yes")),
		                 
		        shiny::selectInput("UnequalCol", "Subjects to Plot:", 
		                  choices = c("All Subjects", "Only Unequal Follow-up", "Only Equal Follow-up"))		                  
		                  
		    ),

		    shiny::mainPanel(
		      shiny::plotOutput("distPlot"), 
		      
		     shiny::column(width = 3,style='padding:0px;height:50px' ,shiny::plotOutput("legend")),
		   	shiny::column(width = 9,style='padding:0px;height:150px', shiny::plotOutput("KMPlots"))
		   	)
		  )

		),#height:100px

	server = function(input, output ) {
		  	output$distPlot <- shiny::renderPlot({
		  	cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
		  	COL = cbPalette[c(3,7,8,4)]
		  	UnequalCens = ifelse(datWIDE$Y_R < datWIDE$Y_D & datWIDE$delta_R == 0, 1, 0)
		
		  	Nobs = length(datWIDE[,1])
			plot(c(),c(), ylim = c(0,max(datWIDE$Y_D)), xlim = c(0, Nobs), xlab = 'Subjects', ylab = 'Time from Baseline', main = input$title)
			
			if(input$UnequalCol=='All Subjects'){
				PlotSubject = rep(1,Nobs)
			}else if(input$UnequalCol=='Only Unequal Follow-up'){
				PlotSubject = UnequalCens
			}else if(input$UnequalCol=='Only Equal Follow-up'){
				PlotSubject = ifelse(UnequalCens==1,0,1)
			}
			segments(x0= c(1:Nobs)[PlotSubject==1], x1 =c(1:Nobs)[PlotSubject==1], 
				y0 = rep(0,Nobs)[PlotSubject==1], y1 = datWIDE$Y_D[PlotSubject==1])
			
			if(input$RecurEvent == TRUE){
				points(c(1:Nobs)[datWIDE$delta_R == 1 & PlotSubject==1], datWIDE$Y_R[datWIDE$delta_R == 1& PlotSubject==1], col = COL[1], pch = 16)
			}
			if(input$RecurCens == TRUE){
				points(c(1:Nobs)[datWIDE$delta_R == 0& PlotSubject==1], datWIDE$Y_R[datWIDE$delta_R == 0& PlotSubject==1], col = COL[2], pch = 16)
			}
			if(input$DeathEvent == TRUE){
				points(c(1:Nobs)[datWIDE$delta_D == 1& PlotSubject==1], datWIDE$Y_D[datWIDE$delta_D == 1& PlotSubject==1], col = COL[3], pch = 16)
			}
			if(input$DeathCens == TRUE){
				points(c(1:Nobs)[datWIDE$delta_D == 0& PlotSubject==1], datWIDE$Y_D[datWIDE$delta_D == 0& PlotSubject==1], col = COL[4], pch = 16)
			}
		
		  })
		  
		  output$legend <- shiny::renderPlot({
				cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
		  		COL = cbPalette[c(3,7,8,4)]
				plot(1,type="n", axes=FALSE, yaxt = 'n', xlab = '', ylab = '')
				par(xpd = T)
				legend(x='center', fill = c('black',COL), legend = c('Follow-up for Death', 'Recurrence Event', 'Recurrence Censoring', 'Death Event', 'Death Censoring' ), 
				horiz=F, cex =0.8)
				par(xpd = F)
		
		  })
		  
		  output$KMPlots <- shiny::renderPlot({
				par(mfrow=c(1,2))
				plot(survival::survfit(survival::Surv(datWIDE$Y_R, datWIDE$delta_R)~1),mark.time = T, xlab = 'Time from Baseline', ylab = 'Event-Free Probability', main = 'KM Plot for Recurrence \n (All Subjects)')
				plot(survival::survfit(survival::Surv(datWIDE$Y_D, datWIDE$delta_D)~1),mark.time = T, xlab = 'Time from Baseline', ylab = 'Event-Free Probability', main = 'KM Plot for Death \n (All Subjects)')
		
		
		  })
		  
		  
	}
	)
shiny::runApp(app)

}
