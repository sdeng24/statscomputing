#' Weighted leveraging Regression
#'
#' @param y response variable
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector
#' @param r number of subsamples
#' @param method "weighted" or "unif"
#'
#' @return a fitted linear regression
#' @export
#'
#' @examples
#' n=500
#' x<-rt(500,df=6)
#' noise<-rnorm(500)
#' y<--x+noise
#' algo_leverage(y,x,r=100,method="weighted")
#'
algo_leverage<-function(y,x,r,method){

 n=length(y)
 if(method=="weighted")
  {
    X = cbind(rep(1, n), x)
    H = X %*% solve(t(X) %*% X) %*% t(X)
    h1<-diag(H)
    h<-h1/sum(h1)
    subsamples<-sample(1:n,r,replace=T,prob=h)
    lmfit = lm(y[subsamples] ~ x[subsamples], weights=1/h[subsamples])
  }
  else{
    if (method=='unif') {
      subsamples = sample(1:n, r)
      lmfit = lm(y[subsamples] ~ x[subsamples])
     } else {
      stop("Method Input Wrong")
     }
  }
    return(lmfit)
}

