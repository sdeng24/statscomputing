#' Title
#'
#' @param X input matrix, of dimension nobs x nvars; each row is an observation vector
#' @param y response variable
#' @param nlambda The number oflambdavalues - default is 100
#' @param alpha  The elasticnet mixing paramete
#' @param thresh  Convergence threshold for coordinate descent.
#' @param maxiter
#'
#' @return length(lambda)matrix of coefficients
#' @export
#'
#' @examples
elasso<-function(X,y,nlambda=100,alpha,thresh=1e-05,maxiter=10000)
{  n = nrow(X)
   p = ncol(X)

   x1<-scale(X,scale=F)
   sdxinv = 1/sqrt(colSums(x1^2)/(n - 1))
   xx = scale(x1)

   ym = mean(Y)
   yy = Y - ym
   if(n>=p){
    ratio=0.001
    }
   else{
     ratio=0.01
    }
   if(alpha>0){
   lambda.max = max(abs(crossprod(xx, yy/n)))/alpha}
   else{
   max(abs(crossprod(xx, yy/n)))/0.01
   }

   lambda.min.value=lambda.max*ratio
   lambda.s= exp(seq(log(lambda.max), log(lambda.min.value),
                  length = nlambda))
   beta_list = list()
   beta0<-rep(0,p)

   for (i in 1:nlambda)
    {
     denominator<-1+lambda.s[i]*(1-alpha)
     beta<-beta0+0.01
     iter=0
    while(max(abs(beta-beta0))>thresh & iter<maxiter){
      iter=iter+1
      beta<-beta0
      for (j in 1:p){
        const = mean(yy - xx %*% beta0)
        r = yy - xx %*% beta0 - const
        term1 = 2*(sum(r*xx[,j])/n + beta0[j])
        beta0[j] = sthreshold(term1, lambda.s[i]* alpha)/(2*denominator)
       }
     }
    beta_list[[i]] = beta0
    }
   beta_e = sapply(beta_list, function(x){x *sdxinv})

   return(list(beta = beta_e, lambda = lambda.s))
}

sthreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  if (x>lambda){
    return (x-lambda)}
  else {
    if (x< (-lambda)){
      return (x+lambda)}
    else {
      return (0); }
  }
}
