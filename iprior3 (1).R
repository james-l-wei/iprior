library(iprior)
library(pacman)

pack.name <- "iprior"

hidden <- setdiff(p_funs(pack.name, TRUE, character.only = TRUE), p_funs(pack.name, character.only = TRUE))

invisible(lapply(hidden, function(x) {
  
  a <- strtrim(x, 1) == "%" 
  b <- substring(x, nchar(x)) == "%"
  
  if (a && b) {
    x2 <- paste0("`", x, "`")
  } else {
    x2 <- x
  }
  
  try(assign(x, eval(parse(text=paste(pack.name,":::",x2,sep=""))), 
             envir = .GlobalEnv))
}))



####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################


iprior3<-function (y, ..., kernel = "linear", method = "direct", 
                   control = list(), interactions = NULL, est.lambda = TRUE, 
                   est.hurst = FALSE, est.lengthscale = FALSE, est.offset = FALSE, 
                   est.psi = TRUE, fixed.hyp = NULL, lambda = 1, psi = 1, nystrom = FALSE, 
                   nys.seed = NULL, model = list(), train.samp, test.samp,inv.steps,
                   inv.tol,inv.method,Mem.Var,stop.crit,theta0=NULL) 
{
  if (is.ipriorKernel(y)) {
    mod <- y
  }
  else {
    mod <- kernL(y = y, ..., kernel = kernel, interactions = interactions, 
                 est.lambda = est.lambda, est.hurst = est.hurst, est.lengthscale = est.lengthscale, 
                 est.offset = est.offset, est.psi = est.psi, fixed.hyp = fixed.hyp, 
                 lambda = lambda, psi = psi, nystrom = nystrom, nys.seed = nys.seed, 
                 model = model, train.samp = train.samp, test.samp = test.samp)
  }
  if (is.categorical(mod)) {
    warning("Categorical responses loaded. Consider using iprobit package.", 
            call. = FALSE)
  }
  method <- tolower(method)
  method <- match.arg(method, c("direct", "em", 
                                "fixed", "canonical", "mixed"))
  control_ <- list(maxit = 500, em.maxit = 20, par.maxit = 5, 
                   stop.crit = stop.crit, theta0 = NULL, silent = FALSE, report = 10, 
                   psi.reg = FALSE, restarts = 0, no.cores = parallel::detectCores(), 
                   optim.method = "L-BFGS", omega = 0)
  control <- update_control(control, control_)
  control.optim <- list(fnscale = -2, trace = ifelse(isTRUE(control$silent), 
                                                     0, 1), maxit = max(0, control$maxit - 1), REPORT = control$report, 
                        factr = control$stop.crit/.Machine$double.eps)
  # if (is.null(control$theta0)) {
  #   theta0 <- rnorm(mod$thetal$n.theta)
  # }
  # else {
  #   theta0 <- control$theta0
  #   if (length(theta0) != mod$thetal$n.theta) {
  #     stop(paste("Incorrect number of parameters specified. Should be", 
  #                mod$thetal$n.theta))
  #   }
  # }
  if (as.numeric(control$restarts) >= 1) {
    res <- iprior_parallel(mod, method, control)
    res$est.method <- paste0(gsub("\\.", "", 
                                  res$est.method), " with random restarts.")
  }
  else {
    est.method <- iprior_method_checker(mod, method)
    if (est.method["fixed"]) {
      res <- iprior_fixed(mod)
      res$est.method <- "Estimation with fixed hyperparameters."
      res$est.conv <- "Convergence not assessed."
    }
    else if (est.method["canonical"]) {
      res <- iprior_canonical(mod, theta0, control.optim)
      res$est.method <- "Efficient canonical method."
    }
    else {
      if (est.method["em.closed"]) {
        res <- iprior_em_closed3(mod, control$maxit, control$stop.crit, 
                                 control$silent, theta0, omega = control$omega,
                                 inv.steps=inv.steps,inv.tol=inv.tol,
                                 inv.method=inv.method,Mem.Var=Mem.Var)
        res$est.method <- "Closed-form EM algorithm."
      }
      if (est.method["em.reg"]) {
        res <- iprior_em_reg(mod, control$maxit, control$stop.crit, 
                             control$silent, theta0)
        res$est.method <- "Regular EM algorithm."
      }
      if (est.method["direct"]) {
        res <- iprior_direct3(mod, loglik_iprior, theta0, 
                              control.optim, control$optim.method)
        res$est.method <- "Direct optimisation method."
      }
      if (est.method["nystrom"]) {
        res <- iprior_direct3(mod, loglik_nystrom, theta0, 
                              control.optim, control$optim.method)
        res$est.method <- "Nystrom approximated optimisation."
      }
      if (est.method["mixed"]) {
        res <- iprior_mixed3(mod, theta0, control$em.maxit, 
                             control$stop.crit, control$silent, control.optim, 
                             control$optim.method, inv.steps=inv.steps,inv.tol=inv.tol,
                             inv.method=inv.method,Mem.Var=Mem.Var)
        res$est.method <- paste0("EM algorithm (", 
                                 control$em.maxit, " steps) + direct minimisation.")
      }
      if (res$conv == 0) 
        res$est.conv <- paste0("Converged to within ", 
                               control$stop.crit, " tolerance.")
      else if (res$conv == 1) 
        res$est.conv <- "Convergence criterion not met."
      else res$est.conv <- res$message
    }
  }
  res$ipriorKernel <- mod
  
  
  res$coefficients <- reduce_theta(res$param.full, mod$estl)$theta.reduced
  tmp <- mse_iprior(mod$y, y.hat = res$y.hat)
  res$fitted.values <- tmp$y + get_intercept(mod)
  names(res$fitted.values) <- attr(mod$y, "dimnames")[[1]]
  res$residuals <- tmp$resid
  res$train.error <- tmp$train.error
  if (!is.null(mod$y.test) & !is.null(mod$Xl.test)) {
    res$test <- structure(predict_iprior(res, mod$Xl.test, 
                                         mod$y.test), class = "ipriorPredict")
  }
  res$method <- method
  res$control <- control
  res$call <- fix_call_default(match.call(), "iprior")
  res$ipriorKernel$call <- fix_call_default(match.call(), "kernL")
  res$LogLik
  class(res) <- "ipriorMod"
  res
}



####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################



iprior_mixed3<-function (mod, theta0 = NULL, em.maxit = 20, stop.crit = 1e-05, 
                         silent = FALSE, control.optim = list(), optim.method =
                           "L-BFGS",inv.steps,inv.tol,inv.method,Mem.Var) 
{
  control.optim_ <- list(fnscale = -2, trace = ifelse(isTRUE(silent), 
                                                      0, 1), maxit = 500, REPORT = 10, factr = 1e+07)
  control.optim <- update_control(control.optim, control.optim_)
  if (!isTRUE(silent)) 
    cat(paste0("Running ", em.maxit, " initial EM iterations\n"))
  start.time <- Sys.time()
  em.method <- iprior_method_checker(mod, "em")
  if (em.method["em.closed"]) {
    tmp <- iprior_em_closed3(mod, em.maxit, stop.crit, silent, 
                             theta0, mixed = TRUE, inv.steps=inv.steps,
                             inv.tol=inv.tol,inv.method=inv.method,Mem.Var=Mem.Var)
  }
  if (em.method["em.reg"]) {
    tmp <- iprior_em_reg(mod, em.maxit, stop.crit, silent, 
                         theta0, mixed = TRUE)
  }
  if (!isTRUE(silent)) 
    cat("Now switching to direct optimisation\n")
  
  tmp$theta<-abs(tmp$theta)
  
  res <- iprior_direct3(mod, loglik_iprior, tmp$theta, control.optim, 
                        optim.method)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)
  res$time <- time.taken
  res$start.time <- start.time
  res$end.time <- end.time
  res$loglik <- c(tmp$loglik, res$loglik)
}



####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################


iprior_em_closed3<-function (mod, maxit = 500, stop.crit = 1e-05, silent = FALSE, 
                             theta0 = NULL, lambda0 = NULL, psi0 = NULL, mixed = FALSE, 
                             omega = 0,inv.steps,inv.tol,inv.method, Mem.Var) 
{
  iprior.env <- environment()
  list2env(mod, iprior.env)
  list2env(BlockBStuff, iprior.env)
  environment(BlockB) <- iprior.env
  environment(em_loop_logical3) <- iprior.env
  maxit <- max(1, maxit)
  
  #######################################
  # theta0 <- c(rnorm(mod$thetal$n.theta)
  
  theta0 <- c(rnorm(2,0,100),rnorm(1))
  psi <- theta_to_psi(theta0, mod)
  lambda.new <- lambda <- theta_to_collapsed_param(theta0, 
                                                   mod)[seq_len(mod$p)]
  niter <- 0
  err.path<-c()
  LogLik.Path<-c()
  lambda.path<-matrix(rep(NA,1100),550,2)
  
  loglik <- rep(NA, maxit)
  Hl <- expand_Hl_and_lambda(Hl, rep(1, p), intr, intr.3plus)$Hl
  if (!silent) 
    pb <- txtProgressBar(min = 0, max = maxit, style = 1)
  start.time <- Sys.time()
  
  while (em_loop_logical3()) {
    Hlam <- get_Hlam(mod, lambda, theta.is.lambda = TRUE)
    eigen_Hlam(Hlam, iprior.env)
    z <- psi * u^2 + 1/psi
    zinv.Vt <- t(V)/z
    Vy.inv.y <- as.numeric(crossprod(y, V) %*% zinv.Vt)
    
    # Modifications begin
    Sigma_theta<-psi*Hlam%*%Hlam+(1/psi)*diag(n)
    if(niter==0){
      Sigma_inverse<-(norm(Sigma_theta,type="F")^-2)*Sigma_theta
    }else{
      Sigma_inverse<-Mem.Var*Sigma_inverse+(1-Mem.Var)*(norm(Sigma_theta,type="F")^-2)*Sigma_theta
    }
    if(inv.method==1){
      Sigma_inverse_list<-Approx.Inverse1(Sigma_theta,inv.steps=inv.steps,
                                          inv.tol=inv.tol,V0=Sigma_inverse)
    }else if(inv.method==2){
      Sigma_inverse_list<-Approx.Inverse2(Sigma_theta,inv.steps=inv.steps,
                                          inv.tol=inv.tol,V0=Sigma_inverse)
    }else if(inv.method==3){
      Sigma_inverse_list<-Approx.Inverse3(Sigma_theta,inv.steps=inv.steps,
                                          inv.tol=inv.tol,V0=Sigma_inverse)
    }
    else if(inv.method==4){
      Sigma_inverse_list<-list("V"=solve(Sigma_theta),"err.path"=0)
    }
    Sigma_inverse<-Sigma_inverse_list[["V"]]
    err.path<-c(err.path,Sigma_inverse_list[["err.path"]])
    
    
    w <- psi * Hlam %*% Sigma_inverse%*%y
    W <- Sigma_inverse + tcrossprod(w)
    
    # w <- psi * Hlam %*% Vy.inv.y
    # W <- V %*% zinv.Vt + tcrossprod(w)
    for (k in seq_len(p)) {
      lambda <- expand_Hl_and_lambda(lambda[1:p], lambda[1:p], 
                                     intr, NULL)$lambda
      BlockB(k)
      T1 <- sum(Psql[[k]] * W)
      T2 <- crossprod(y, Pl[[k]]) %*% w - sum(Sl[[k]] * 
                                                W)/2
      lambda.new[k] <- as.numeric(T2/T1)
    }
    lambda <- (1 + omega) * lambda.new - omega * lambda[1:p]
    Hlamsq <- V %*% (t(V) * u^2)
    T3 <- crossprod(y) + sum(Hlamsq * W) - 2 * crossprod(y, 
                                                         Hlam %*% w)
    psi.new <- sqrt(max(0, as.numeric(sum(diag(W))/T3)))
    psi <- (1 + omega) * psi.new - omega * psi
    logdet <- sum(log(z))
    loglik[niter + 1] <- -n/2 * log(2 * pi) - logdet/2 - 
      crossprod(y, Vy.inv.y)/2
    niter <- niter + 1
    
    lambda<-abs(lambda)
    lambda.path[niter,]<-t(lambda[1:2])
    
    print(lambda)
    # if (!silent){
    #   setTxtProgressBar(pb, niter)
    #   print(niter)
    #   print(loglik[niter])
    # }
  }
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)
  y.hat <- get_y.hat(u, V, w)
  if (!isTRUE(mixed)) {
    Vy.inv <- V %*% zinv.Vt
    dVy <- NULL
    for (i in seq_len(p)) {
      dVy[[i]] <- Vy.inv %*% (psi * (2 * lambda[i] * Psql[[i]] +
                                       Sl[[i]]))
    }
    dVy[[p + 1]] <- diag(1/psi, n) - (2/psi^2) * Vy.inv
    Fi <- matrix(0, nrow = p + 1, ncol = p + 1)
    for (i in seq_len(p + 1)) {
      for (j in seq_len(p + 1)) {
        Fi[i, j] <- sum(dVy[[i]] * dVy[[j]])/2
      }
    }
    se <- sqrt(diag(solve(Fi)))
    #se<-NULL
  }
  else {
    se <- NULL
  }
  convergence <- niter == maxit
  param <- kernel_to_param(kernels, lambda[1:p])
  theta <- param_to_theta(param, estl, logpsi = log(psi))$theta
  param.full <- theta_to_collapsed_param(theta, mod)
  if (!silent) {
    close(pb)
    if (isTRUE(mixed)) 
      cat("")
    else if (convergence) 
      cat("Convergence criterion not met.\n")
    else cat("Converged after", niter, "iterations.\n")
  }
  list(theta = theta, param.full = param.full, se = se, loglik = as.numeric(na.omit(loglik)), 
       w = as.numeric(w), y.hat = y.hat, niter = niter, start.time = start.time, 
       end.time = end.time, time = time.taken, convergence = convergence, 
       message = NULL,err.path=err.path,LogLik.Path=LogLik.Path,lambda.path=lambda.path)
}


####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################


iprior_direct3<-function (mod, lik.fn, theta0, control = list(), method = "L-BFGS") 
{
  control_ <- list(fnscale = -2, trace = 1, maxit = 500, REPORT = 10, 
                   factr = 10^-3)
  control <- update_control(control, control_)
  iprior.env <- environment()
  w <- loglik <- NULL
  start.time <- Sys.time()
  
  # mods begin
  control[["maxit"]]<-1
  iter<-0
  lik.old<-99999
  lik.chg<-1
  while(lik.chg>10^-3){
    res <- optim(abs(theta0), lik.fn, object = mod, env = iprior.env,
                 trace = TRUE, method = method, control = control, hessian = TRUE)
    iter<-iter+res$count[1]
    lik.new<-res$value
    lik.chg<-abs(lik.new-lik.old)
    lik.old<-lik.new
  }

  # mods end
  
  # res <- optim(theta0, lik.fn, object = mod, env = iprior.env, 
  #              trace = TRUE, method = method, control = control, hessian = TRUE)
  end.time <- Sys.time()
  time.taken <- as.time(end.time - start.time)
  w <- get_w(u, V, Vy.inv.y, psi)
  y.hat <- get_y.hat(u, V, w)
  tmp <- eigenCpp(-res$hessian)
  u <- tmp$val + 1e-09
  V <- tmp$vec
  Fi.inv <- V %*% (t(V)/u)
  se <- sqrt(diag(Fi.inv))
  se <- convert_se(se, res$par, mod)
  loglik <- as.numeric(na.omit(loglik))
  param.full <- theta_to_collapsed_param(res$par, mod)
  list(theta = res$par, param.full = param.full, se = se, loglik = loglik, 
       w = w, y.hat = y.hat, niter = iter, start.time = start.time, 
       end.time = end.time, time = time.taken, convergence = res$convergence, 
       message = res$message)
}


####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

####################################################################################

em_loop_logical3<-function () 
{
  ll.diff <- loglik[niter] - loglik[niter - 1]
  crit1 <- (niter != maxit)
  crit2 <- (abs(ll.diff) > stop.crit)
  if (niter == 0) {
    return(TRUE)
  }
  else if (niter == 1) {
    return(crit1)
  }
  else {
    if (ll.diff < 0) {
      warning(paste0("Log-likelihood decreased at iteration ", 
                     niter), call. = FALSE)
    }
    return(crit1 & crit2)
  }
}
