
spdmar <- function(p, n) {
  stopifnot(p > 0, p == as.integer(p), n>0, n== as.integer(n))
  
  pars <- lapply(seq_len(p+1), function(i) {
    # A <- matrix(rnorm(n * n), nrow = n, ncol = n)
    # crossprod(A) + diag(n)
    A <- diag(n)
    
  })
  
  #convention pars$[[1]] is the intercept W_0, pars$[[2]] is W_1 (weight for most recent lag), etc... 
  
  structure(
    list(p= p, n=n, pars=pars),
    class = "spdmar"
  )
  
}

predict.spdmar <- function(object, input_series){
  #this is just one-step forward prediction at the moment. 
  
  #input_series shall be pxnxn (the required p input SPD matrices)
  p <- object$p
  n <- object$n
  pars <- object$pars
  
  dims <- dim(input_series)
  
  stopifnot(length(dims)== 3, dims[1]==p, dims[2]==n, dims[3]==n)
  
  # input_series <- time_series[t-p:t,,] 
  #input_series[1,,] contains the least recent matrix, whereas input_series[p,,] the most recent i.e. Sigma(t-1)
  prediction <- pars[[1]]%*%pars[[1]]
    
  for (i in 1:p){
    prediction <- prediction +  (pars[[p+2-i]]  %*% input_series[i,,]) %*% pars[[p+2-i]]
  }
    
  return (prediction)
}

grad.spdmar <- function(object, time_series){
  p <- object$p
  n <- object$n
  pars <- object$pars
  dims <- dim(time_series)
  stopifnot(length(dims) == 3, dims[1] >= p + 1, dims[2] == n,dims[3] == n)
  T_ <- dim(time_series)[1]
  
  cum_gradients <- lapply(seq_len(p + 1), function(i) {
    matrix(0, nrow = n, ncol = n)
  })
  
  for (t_ in (p+1):T_){
    #build lagged series
    input_series = time_series[(t_-p):(t_-1),,,drop = FALSE]
    target = time_series[t_,,]
    E_t = target-predict(object, input_series)
    
    #intercept gradient
    grad_int = -2*(E_t %*% pars[[1]] + pars[[1]] %*% E_t )
    cum_gradients[[1]] = cum_gradients[[1]] + grad_int
    
    for (j in (1:p)){
      Wj = pars[[j+1]]
      Sigma = time_series[t_-j,,]
      gradient = -2*(E_t %*% Wj %*% Sigma + Sigma %*% Wj %*% E_t)
      cum_gradients[[j+1]] = cum_gradients[[j+1]] + gradient
    }
    
  }
  return(lapply(cum_gradients, function(G) {
    G / (T_ - p)
  }))
  
}

frob_loss.spdmar <- function(object, time_series){
  #computes \|E\|_F^2
  
  #time_series shall be txnxn
  p <- object$p
  n <- object$n
  pars <- object$pars
  
  dims <- dim(time_series)
  
  stopifnot(length(dims) == 3, dims[1] >= p + 1, dims[2] == n,dims[3] == n)
  
  T_ <- nrow(time_series)
  
  cum_loss = 0
  for (t_ in (p+1):T_){
    #build lagged series
    input_series = time_series[(t_-p):(t_-1),,,drop = FALSE]
    target = time_series[t_,,]
    E = predict(object, input_series) - target
    frob_ob = sum(E^2)
    cum_loss = cum_loss + frob_ob
  }
  
  return(cum_loss / (T_ - p))
}



print.spdmar <- function(object, ...) {
  cat("SPD Matrix Autoregressive Model\n")
  cat("Order:", object$p, "\n")
  cat("Number of weight matrices including intercept (pars):", object$p+1, "\n")
  cat("Input Dimension:", object$n, "times", object$n, "\n")
}



