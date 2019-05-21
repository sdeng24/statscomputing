
#' Title
#'
#'Apply Gauss-Seidel or Jacobi method to solve the linear equations
#'
#' @param A  the linear matrix of coefficient
#' @param v  the underlying true solution
#' @param method  "G-S" stands for Gauss-Seidel method and "Jacobi" uses Jacobi method
#' @param ncores  1 stands for non-parallel computing and usually ncores<=4
#' @param iter.max the maximal number of iteration
#'
#' @return  relative error: the relative error at each iteration
#' @return  solution
#' @import  foreach
#' @import  doParallel
#'
#' @export
#'
#' @examples
#' n=100
#' alpha=3
#' A<-HWGen(n,alpha)
#' v<-rep(c(1,0),50)
#' solve_ols(A,v,method="Jacobi",ncores=1)
#'

solve_ols<-function(A,v,method,ncores=1,iter.max=10000)
{
  n=dim(A)[2]
  b<-A%*%v
  x<-rep(0,n)
  iter<-0
  relative<-c()
  if(method=="G-S")
  {
    while(iter<iter.max)
    {
      for(i in 1:n)
      {
        x[i] =(b[i]-sum( A[i, -i] * x[-i]))/A[i,i]
      }
      iter=iter+1
      relative<-c(relative,sum(abs(x-v))/sum(v))
    }
  }

  else{
   if(ncores==1){
      while(iter<iter.max)
     {
      x0<-x
      for(i in 1:n)
      {
        x[i] = (b[i]- A[i, -i] %*% x0[-i])/A[i,i]
      }
      iter=iter+1
      relative<-c(relative,sum(abs(x-v))/sum(v))
     }
    }
    else{
      doParallel::registerDoParallel(ncores)
       for(iter in 1:iter.max)
        {
          x0<-x
          x<- foreach(i=1:n,.combine = 'c') %dopar% {
          x[i] = (b[i] - sum(A[i, -i] * x0[-i]))/A[i,i]
        }
         relative<-c(relative,sum(abs(x-v))/sum(v))
        }
      }
  }
  return(list(error=relative,solution=x))
}



#' Title
#' Generating 3 banded matrix to solve the hw1
#'
#' @param n  dimension of the square matrix
#' @param alpha  the diagonal value
#'
#' @return 3 banded matrix
#' @import Matrix
#' @export
#'
#' @examples
HWGen<-function(n,alpha)
{
 row<-matrix(c(rep(alpha,n),rep(-1,n)),nrow=n)
 A<-Matrix::bandSparse(n,k=c(0,1),diag=row,symm=T)
 A<-as.matrix(A)
 return(A)
}

