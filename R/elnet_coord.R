#' Title
#'
#'Solve the an elastic net problem via coordinate descent
#'
#' @param X input matrix, of dimension nobs x nvars; each row is an observation vector
#' @param y response variable
#' @param lambda  default is null
#' @param nlambda The number oflambda values - default is 100
#' @param alpha  The elasticnet mixing paramete
#' @param thresh  Convergence threshold for coordinate descent.
#' @param maxiter the maximal number of iterations in each coordinate descent
#'
#' @return length(lambda)matrix of coefficients
#' @export
#'
#' @examples
#'n=50
#'p=20
#'beta<-c(2,0,-2,0,1,0,-1,0,rep(0,12))
#'X=matrix(rnorm(n*p),n,p)
#'Y<-X%*%beta+rnorm(n)
#'fit1<-elnet_coord(X,Y,alpha=0.3)
#'beta_cd<-t(as.matrix(fit1$beta))
#'pct <- rowSums(abs(beta_cd))
#'matplot(pct,beta_cd,type="l",lty=1,
#'        xlab="|beta|_1",ylab="Coefficients")
#'text(max(pct)+0.2,beta_cd[dim(beta_cd)[1],],1:p,col=1:p,cex=0.7)
#'
#'
elnet_coord<-function(X,y,lambda = NULL,nlambda=100,alpha,thresh=1e-05,maxiter=10000)
{  n = nrow(X)
   p = ncol(X)

   x1<-scale(X,scale=F)
   sdxinv = 1/sqrt(colSums(x1^2)/(n - 1))
   xx = scale(x1)

   ym = mean(Y)
   yy = Y - ym

   if (alpha > 1) {
     warning("alpha >1; set to 1")
     alpha = 1
   }
   if (alpha < 0) {
     warning("alpha<0; set to 0")
     alpha = 0
   }

   if(is.null(lambda))
   {
     if(n>=p){
      ratio=0.001
      }
     else{
      ratio=0.01
      }
     if(alpha>0){
     lambda.max = max(abs(crossprod(xx, yy/n)))/alpha*1.1}
     else{
     lambda.max=max(abs(crossprod(xx, yy/n)))/0.01*1.1
     }
     lambda.min.value=lambda.max*ratio
     lambda.s= exp(seq(log(lambda.max), log(lambda.min.value),
                       length = nlambda))
   }else{
     nlambda=1
     lambda.s=c(lambda)
    }


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
        term1 = (sum(r*xx[,j])/n + beta0[j])
        beta0[j] = sthreshold(term1, lambda.s[i]* alpha)/denominator
       }
     }
    beta_list[[i]] = beta0
    }
   beta_e = sapply(beta_list, function(x){x *sdxinv})

   return(list(beta = beta_e, lambda = lambda.s))
}

#' plot the elnet_coord
#'
#' @param fit the fitted object via elnet_coord
#'
#' @return a plot
#' @export
#'
#' @examples
#'
plot_el<-function(fit){

  beta_cd<-t(as.matrix(fit$beta))
  pct <- rowSums(abs(beta_cd))

  matplot(pct,beta_cd,type="l",lty=1,
          xlab="|beta|_1",ylab="Coefficients")
  text(max(pct)+0.2,beta_cd[dim(beta_cd)[1],],1:p,col=1:p,cex=0.7)
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
