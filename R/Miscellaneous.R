#' @export

mSample = function(Y){
	return(sample(x=c(0,1), size = 1, prob = c(1-Y,Y)))
}

 #' @export

expit = function(x){exp(x)/(1+exp(x))}


 #' @export

mHPropose = function(Y){
	minimum = Y[1]
	maximum = Y[2]
	return(runif(n=1,min = minimum, max = maximum))
}




 #' @export

#Outputs means and variance
RubinMe <- function(means,vars,impNum)
{
	N = length(means[,1])
    M1 = means
    V1 = vars
    
    Q1 = matrix(NA,N,1)
    Q1 = (1/impNum)*apply(M1,1,sum)
    
    U1 = matrix(NA,N,1)
    U1 = (1/impNum)*apply(V1,1,sum)
    
    B1 = matrix(0,N,N)
    for(i in 1:impNum)
    {
        B1 = B1 + (M1[,i]-Q1)%*%t(M1[,i]-Q1)
    }
    B1 = (1/(impNum-1)) * B1
    
    CovMat1 = U1 + (1 + (1/impNum)) * as.vector(diag(B1))
    r = (1+(1/impNum))*diag(B1)/U1
	v = (impNum-1)*((1+(1/r))^2)

    return(data.frame(mean = Q1 , var = CovMat1, v = v))
}

