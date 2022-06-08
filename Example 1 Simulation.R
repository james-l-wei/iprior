setwd("C:/! Files/Seraph Research/Knowledge Base/Statistics and Programming/MSc Statistics/ST499 Dissertation/Code/New code")


# library(plotly)

seed<-1 # 1,4,19
set.seed(seed)

n<-200
x_dot1<-runif(n,-1,1)
x_dot2<-runif(n,0,24)

# Example 1: lambda_1=10, lambda_2=0.1, psi=0.1
# Example 2: lambda_1=50, lambda_2=0.1, psi=0.1
# Example 3: lambda_1=10, lambda_2=0.5, psi=0.1
# Example 4: lambda_1=10, lambda_2=0.1, psi=0.5

lambda_1<-10
lambda_2<-0.1
psi<-0.1

h_lambda<-function(x_i1,x_j1,x_i2,x_j2,lambda1,lambda2){
  1+0.5*lambda1*((1/2)*(x_i1^2+1)+(1/2)*(x_j1^2+1)-abs(x_i1-x_j1)-(2/3))+
    0.5*lambda2*((1/24)*(x_i2^2-24*x_i2+288)+(1/24)*(x_j2^2-24*x_j2+288)-
                   abs(x_i2-x_j2)-8)
}

w<-rnorm(n,mean=0,sd=sqrt(psi))
f<-rep(0,n)
for(i in 1:n){
  for(j in 1:n)
    f[i]<-f[i]+w[j]*h_lambda(x_dot1[i],x_dot1[j],x_dot2[i],x_dot2[j],
                             lambda_1,lambda_2)
}
y<-f+rnorm(n,mean=0,sd=1/sqrt(psi))

Data.Ex1<-data.frame(y=y,x1=x_dot1,x2=x_dot2)

# plot_ly(x=x_dot1,y=x_dot2,z=y,type="scatter3d",mode="markers")%>%
#   layout(
#     title="",
#     scene=list(
#       xaxis=list(title = "x_i1"),
#       yaxis=list(title = "x_i2"),
#       zaxis=list(title = "y")
#     ))

n<-length(y)

set.seed(seed)
mod<-kernL(y~x1+x2,Data.Ex1,kernel="fbm",est.hurst=FALSE)

