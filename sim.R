source("~/Desktop/Research Projects/SPDMAR/code/spdmar/utils.R")

n=5
T_ = 10000  #nobs the estimator sees
burn_in = 100
p=1
n_traj = 1
n_epochs = 500
l_rate = 0.009
n_experiments = 1
nsim = T_+burn_in #nobs generated per simulation




buildWeights <- function(n,p, target_rho = 0.9){

  pars <- lapply(seq_len(p+1), function(i) {
      A <- matrix(rnorm(n * n), nrow = n, ncol = n)
      A <- (A+t(A))*0.5
      A <- higham_approx(A)
      if (i > 1) {
        #spectral radius 
        rho <- max(abs(eigen(A, symmetric = TRUE)$values))
        
        #shrink 
        if (rho >= target_rho) {
          A <- A * target_rho / rho
        }
      }
      return(A)
  })
    
  return (pars)
  
}


for (exp_id in 1:n_experiments){
  #generate the ground truth SPD parameters (DGP)
  set.seed(exp_id)
  Ws <- buildWeights(n, p) #W_0,W_1,...,W_p
  Y_final <- array(0,dim = c(n_traj, T_, n, n))
  
  for(v in 1:n_traj){
    seed <- exp_id * n_traj + v
    set.seed(seed)
    
    first_p_obs <- lapply(seq_len(p), function(i) {
      A <- matrix(rnorm(n * n), nrow = n, ncol = n)
      A <- (A+t(A))*0.5
      A <- higham_approx(A)
    })
    
    Yt <- array(0,dim = c(nsim, n, n))
    for (j in 1:p) {
      Yt[j, , ] <- first_p_obs[[j]]
    }
    
  
    for (i in (p+1):nsim) {
      #generate the innovation
      Z <- matrix(rnorm(n * n), nrow = n, ncol = n) #E[ZZ^T]=n*I
      mu <- Ws[[1]] %*% Ws[[1]]
      E_t <- (Ws[[1]] %*% Z %*% t(Ws[[1]] %*% Z)) / n
      eps_t <- E_t - mu
      
      
      prediction <- Ws[[1]] %*% Ws[[1]]
      for (j in 1:p) {
        W <- Ws[[p + 2 - j]]
        prediction <- prediction + (W %*% Yt[i-j, , ]) %*% W
      }
      
      Yt[i,, ] <- prediction + eps_t
      
      Yt[i, , ] <- (Yt[i, , ] + t(Yt[i, , ])) / 2
    }
    
    train_set <- Yt[(burn_in+1):nsim,, , drop=FALSE]
    
    #discard burnin
    Y_final[v,,,] <- train_set
    
    
    
    ##PGD
    model <- spdmar(p, n)
  
    for (epoch in 1:n_epochs){
      
      #compute the gradient
      euclidean_gradients <- grad.spdmar(model, train_set)
      
      #update the weights
      for (i in 1:(p + 1)) {
        G <- euclidean_gradients[[i]]
        W <- model$pars[[i]]
        step_size <- l_rate * norm(W, "F") / norm(G, "F")
        model$pars[[i]] <- model$pars[[i]] - step_size * euclidean_gradients[[i]]
        
        #projection via higham approx
        model$pars[[i]] <- higham_approx(model$pars[[i]])
        #model$pars[[i]] <- model$pars[[i]]/ norm(model$pars[[i]], "F")
      }
      
      loss <- frob_loss.spdmar(model, train_set)
      cat("Epoch:", epoch, "Train Loss:", loss, "\n")
      
    }

    
    
    ##Riemannian optm
    model_r <- spdmar(p, n)
    

    for (epoch in 1:n_epochs){
      
      #compute the gradient
      euclidean_gradients <- grad.spdmar(model_r, train_set)
      
      #update the weights
      for (i in 1:(p + 1)) {
        G <- euclidean_gradients[[i]]
        W <- model_r$pars[[i]]
        
        sqrt_W <- matrix_sqrt(W)
        inv_sqrt_W <- matrix_inv_sqrt(W)
        
        riemannian_gradient <- sqrt_W %*% ((G+t(G))/2) %*% sqrt_W
        g_norm <- norm(riemannian_gradient, "F")
        
        step_size <- l_rate * norm(W, "F") / g_norm
        P <- step_size * riemannian_gradient 
        
        
        
        A <- -1*inv_sqrt_W %*% P %*% inv_sqrt_W
        
        W_new <- sqrt_W %*%
          expm(A) %*%
          sqrt_W
        
        
        W_new <- (W_new+t(W_new))/2
        model_r$pars[[i]] <- W_new
        
        
      }
      
      loss <- frob_loss.spdmar(model_r, train_set)
      cat("Epoch:", epoch, "Train Loss:", loss, "\n")
      
    }
    

    for (i in 1:(p + 1)) {
      
      cat("\n====================================\n")
      cat("Parameter W_", i - 1, "\n", sep = "")
      
      cat("\nGround truth:\n")
      print(Ws[[i]])
      
      cat("\nEstimated PGD (frob distance) :\n")
      print(sum((Ws[[i]]-model$pars[[i]])^2))
      
      cat("\nEstimated Riemannian :\n")
      print(sum((Ws[[i]]-model_r$pars[[i]])^2))
    }
    

    
    
  }
  
  #plot the data we got
  
  
  

}

