
#' COVIMPUTEINITIALIZE
#' @description The function COVIMPUTEINITIALIZE will initialize the missing covariate values. An example function is included in the package, but the user must specify their own function when applying MultiCure to datasets with covariate missingness. This function must have input/output as described below.
#' @param Cov matrix of covariates, same as in MultiCure 
#' @param CovMissing: A matrix indicating which elements of Cov are missing. 
#'
#' @return Cov a initialized version of the covariate matrix
#' @details The example code included in the package initializes missing covariate X2 in the Multistate cure model example
#' @export

COVIMPUTEINITIALIZE = function(Cov, CovMissing){
	### Code for initializing missing covariate X2 in the example
	Cov[CovMissing[,'X2'] == T,'X2'] = sample(x=Cov[CovMissing[,'X2']==F,'X2'], size = sum(as.numeric(CovMissing[,'X2'])), replace = T)	
    	return(Cov)
}

