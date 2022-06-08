library(beepr)


# Experiment on memory weights; set inv.method 1, 2, 3

mem_weights<-seq(0,1,0.05)

Results.Matrix<-matrix(NA,20,length(mem_weights))
Results.Matrix2<-matrix(NA,20,length(mem_weights))
for(i in 1:20){ 
  for(j in 1:length(mem_weights)){
    rm(mod.fit)
    set.seed(seed)
    try(mod.fit<-iprior3(mod,method="em",inv.steps=i,inv.tol=10^-10,inv.method=2,
                         Mem.Var=mem_weights[j],stop.crit=0.001))
    try(Results.Matrix[i,j]<-mod.fit$niter)
    try(Results.Matrix2[i,j]<-summary(mod.fit)$train.rmse)
  }
}
write.csv(Results.Matrix,file="Results.Matrix.csv")
write.csv(Results.Matrix2,file="Results.Matrix2.csv")
beep("fanfare")

#############################################################################################

#############################################################################################

#############################################################################################

#############################################################################################

# Error tolerances

inv_tol.vector<-rep(1,30)
for(i in 1:length(inv_tol.vector)){
  inv_tol.vector[i]<-inv_tol.vector[i]*2^(-i)
}

Summary<-list()

nx<-10000
ErrorPath<-rep(NA,nx)
LogLikPath<-rep(NA,nx)

for(i in 1:30){ 
  rm(mod.fit)
  rm(err.path2)
  rm(LogLikPath2)
  rm(n1)
  rm(n2)
  set.seed(seed)
  start.time<-Sys.time()
  try(mod.fit<-iprior3(mod,method="em",inv.steps="Exact",inv.tol=inv_tol.vector[i],
                       inv.method=1,Mem.Var=0.4,stop.crit=0.001))
  end.time<-Sys.time()
  try(Summary<-rbind(Summary,list("Test"=i,"lambda1"=mod.fit$coef[1][[1]],
                                  "lambda2"=mod.fit$coef[2][[1]],
                                  "psi"=mod.fit$coef[3][[1]],
                                  "rmse"=summary(mod.fit)$train.rmse,
                                  "Runtime"=end.time-start.time,
                                  "Steps"=mod.fit$niter)))
  try(err.path2<-mod.fit$err.path)
  try(length(err.path2)<-nx-1)
  try(ErrorPath<-cbind(ErrorPath,c(i,err.path2)))
  
  try(LogLikPath2<-mod.fit$LogLikPath)
  try(length(LogLikPath2)<-nx-1)
  try(LogLikPath<-cbind(LogLikPath,c(i,LogLikPath2)))
}
colnames(ErrorPath)<-NULL
colnames(LogLikPath)<-NULL
write.csv(Summary,file="Summary.csv")
write.csv(ErrorPath,file="ErrorPath.csv")
write.csv(LogLikPath,file="LogLikPath.csv")
beep("fanfare")




#############################################################################################

#############################################################################################

#############################################################################################

#############################################################################################

# Initialisations


lambdapath<-c()

for(i in 1:30){ 
  rm(mod.fit)
  start.time<-Sys.time()
  try(mod.fit<-iprior3(mod,method="em",inv.steps="Exact",inv.tol=2^(-3),
                       inv.method=1,Mem.Var=0.4,stop.crit=0.001,
                       theta0<-c(rnorm(1,0,100),rnorm(1,0,100),rnorm(1))))
  end.time<-Sys.time()
  try(Summary<-rbind(Summary,list("Test"=i,"lambda1"=mod.fit$coef[1][[1]],
                                  "lambda2"=mod.fit$coef[2][[1]],
                                  "psi"=mod.fit$coef[3][[1]],
                                  "rmse"=summary(mod.fit)$train.rmse,
                                  "Runtime"=end.time-start.time,
                                  "Steps"=mod.fit$niter)))
  try(lambdapath<-cbind(lambdapath,mod.fit$lambda.path))
}
write.csv(lambdapath,file="lambdapath.csv")
beep("fanfare")