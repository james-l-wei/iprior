# # Example data 3
# 
set.seed(1234)
A2<-Sigma_theta
V02<-(norm(A2,type="F")^-2)*A2
V02<-0.5*Sigma_inverse+0.5*(norm(A2,type="F")^-2)*A2

# Don't work
# V02<-abs(Sigma_inverse)
# V02<-diag(diag(Sigma_inverse))
# V02<-solve(Sigma_theta)+rnorm(n,mean=0,sd=0.02)
# V02<-Sigma_inverse
# V02<-A2
# V02<-Sigma_inverse/(n)
# V02<-(norm(Sigma_inverse,type="F")^-2)*Sigma_inverse
# work
# V02<-(norm(A2,type="F")^-2)*A2

A3<-matrix(1,100,100)
for(i in 1:100){
  for (j in 1:(i-1)){
    A3[i,j]<- -1
  }
  A3[i,i]<-100+i
}
tol<-10^-10
V03<-(norm(A,type="F")^-2)*A3

start.time<-Sys.time()
Results<-Approx.Inverse1(A2,"Exact",tol,V02)
V<-Results[["V"]]
err.path<-Results[["err.path"]]
end.time<-Sys.time()
AI.time.taken3<-end.time-start.time

#############################################################################################

#############################################################################################

#############################################################################################

#############################################################################################

## Schulz Algorithm ##

Approx.Inverse1<-function(A,inv.steps,inv.tol,V0){
  
  # Initialisation
  n<-dim(A)[1]
  V<-V0
  err.path<-c()
  err1<-norm(diag(n)-A%*%V,type="F")
  err<-err1
  
  if(inv.steps=="Exact"){
    # Ensures that if the number of steps is set to "Exact", a minimum of m1 loops are performed
    # so that "Exact" runtime is always greater than any other runtime.
    m1<-ceiling(n^0.5)
    for(i in 1:m1){ 
      if(err>err1*inv.tol){
        V<-V%*%(2*diag(n)-A%*%V)
        err<-norm(diag(n)-A%*%V,type="F")
        err.path<-c(err.path,err)
      }
    }
    while(err>err1*inv.tol){
      V<-V%*%(2*diag(n)-A%*%V)
      err<-norm(diag(n)-A%*%V,type="F")
      err.path<-c(err.path,err)
    }
  }else{
    # If number of steps is not set to "Exact", simply repeat the update step the required number
    # of times, with empty loops when the tolerance is reached
    for(i in 1:inv.steps){
      if(err>err1*inv.tol){
        V<-V%*%(2*diag(n)-A%*%V)
        err<-norm(diag(n)-A%*%V,type="F")
        err.path<-c(err.path,err)
      }
    }
  }
  return(list("V"=V,"err.path"=err.path))
}


#############################################################################################

#############################################################################################

#############################################################################################

#############################################################################################


## Modified Root-Finding Algorithm ##

Approx.Inverse2<-function(A,inv.steps,inv.tol,V0){
  
  # Initialisation
  n<-dim(A)[1]
  V<-V0
  err.path<-c()
  err1<-norm(diag(n)-A%*%V,type="F")
  err<-err1
  
  if(inv.steps=="Exact"){
    # Ensures that if the number of steps is set to "Exact", a minimum of m1 loops are performed
    # so that "Exact" runtime is always greater than any other runtime.
    m1<-ceiling(n^0.5)
    for(i in 1:m1){ 
      if(err>err1*inv.tol){
        V<-0.5*V%*%(9*diag(n)-A%*%V%*%(16*diag(n)-A%*%V%*%(14*diag(n)-A%*%V%*%(6*diag(n)-A%*%V))))
        err<-norm(diag(n)-A%*%V,type="F")
        err.path<-c(err.path,err)
      }
    }
    while(err>err1*inv.tol){
      V<-0.5*V%*%(9*diag(n)-A%*%V%*%(16*diag(n)-A%*%V%*%(14*diag(n)-A%*%V%*%(6*diag(n)-A%*%V))))
      err<-norm(diag(n)-A%*%V,type="F")
      err.path<-c(err.path,err)
    }
  }else{
    # If number of steps is not set to "Exact", simply repeat the update step the required number
    # of times, with empty loops when the tolerance is reached
    for(i in 1:inv.steps){
      if(err>err1*inv.tol){
        V<-0.5*V%*%(9*diag(n)-A%*%V%*%(16*diag(n)-A%*%V%*%(14*diag(n)-A%*%V%*%(6*diag(n)-A%*%V))))
        err<-norm(diag(n)-A%*%V,type="F")
        err.path<-c(err.path,err)
      }
    }
  }
  return(list("V"=V,"err.path"=err.path))
}



#############################################################################################

#############################################################################################

#############################################################################################

#############################################################################################

## Hyperpower Iterative Algorithm ##

Approx.Inverse3<-function(A,inv.steps,inv.tol,V0){
  
  # Initialisation
  n<-dim(A)[1]
  V<-V0
  err.path<-c()
  err1<-norm(diag(n)-A%*%V,type="F")
  err<-err1
  
  #Constants
  c1<-(1/4)*(sqrt(27-2*sqrt(93))+1)
  c2<-(1/4)*(1-sqrt(27-2*sqrt(93)))
  c3<-(1/496)*(5*sqrt(93)-93)
  d1<-(1/496)*(-93-5*sqrt(93))
  d2<- -sqrt(93)/4
  mu<-3/8
  psi<-321/1984
  
  R<-diag(n)-A%*%V
  
  if(inv.steps=="Exact"){
    # Ensures that if the number of steps is set to "Exact", a minimum of m1 loops are performed
    # so that "Exact" runtime is always greater than any other runtime.
    m1<-ceiling(n^0.5)
    for(i in 1:m1){ 
      if(err>err1*inv.tol){
        R2<-R%*%R
        R4<-R2%*%R2
        M<-(diag(n)+c1*R2+R4)%*%(diag(n)+c2*R2+R4)
        T.<-M+c3*R2
        S<-M+d1*R2+d2*R4
        V<-V%*%((diag(n)+R)%*%(T.%*%S)+mu*R2+psi*R4)
        #V<-V%*%A%*%V
        
        R<-diag(n)-A%*%V
        err<-norm(R,type="F")
        err.path<-c(err.path,err)
      }
    }
    while(err>err1*inv.tol){
      R2<-R%*%R
      R4<-R2%*%R2
      M<-(diag(n)+c1*R2+R4)%*%(diag(n)+c2*R2+R4)
      T.<-M+c3*R2
      S<-M+d1*R2+d2*R4
      V<-V%*%((diag(n)+R)%*%(T.%*%S)+mu*R2+psi*R4)
      #V<-V%*%A%*%V
      
      R<-diag(n)-A%*%V
      err<-norm(R,type="F")
      err.path<-c(err.path,err)
    }
  }else{
    # If number of steps is not set to "Exact", simply repeat the update step the required number
    # of times, with empty loops when the tolerance is reached
    for(i in 1:inv.steps){
      if(err>err1*inv.tol){
        R2<-R%*%R
        R4<-R2%*%R2
        M<-(diag(n)+c1*R2+R4)%*%(diag(n)+c2*R2+R4)
        T.<-M+c3*R2
        S<-M+d1*R2+d2*R4
        V<-V%*%((diag(n)+R)%*%(T.%*%S)+mu*R2+psi*R4)
        #V<-V%*%A%*%V
        
        R<-diag(n)-A%*%V
        err<-norm(R,type="F")
        err.path<-c(err.path,err)
      }
    }
  }
  return(list("V"=V,"err.path"=err.path))
}

