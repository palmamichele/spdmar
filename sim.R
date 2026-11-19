source("~/Desktop/Research Projects/SPDMAR/code/spdmar/utils.R")
source("~/Desktop/Research Projects/SPDMAR/code/spdmar/spdmar.R")
library("expm")


n=3
T_list = c(200,500,1000)  #nobs the estimator sees
#assume T_list is non decreasing, last element in T_list is the largest one

burn_in = 100
p=1
n_traj = 100
n_epochs = 10000
l_rate_pgd = 0.09
l_rate_rm = 0.09
n_methods=3

n_experiments = 1 #for generating the pars
N = n*(n+1)*(p+1)/2
n_cols_per_T <- 3 * n_methods 



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


vech_lt <- function(S) S[lower.tri(S, diag = TRUE)]

unvech_lt <- function(v, n) {
  S <- matrix(0, n, n)
  S[lower.tri(S, diag = TRUE)] <- v
  S + t(S) - diag(diag(S), n)
}

# SPD square root; eigenvalues below eps are clipped to eps first
# (only active when the argument is not positive definite)
spd_sqrt <- function(S, eps = 1e-8) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  lam <- pmax(e$values, eps)
  e$vectors %*% diag(sqrt(lam), length(lam)) %*% t(e$vectors)
}

# closed-form OLS of the VAR(p) on vech(Sigma_t), then recovery of W_0,...,W_p
# train_set: array T x n x n
fit_var_vech <- function(train_set, p) {
  T_ <- dim(train_set)[1]
  n  <- dim(train_set)[2]
  m  <- n * (n + 1) / 2
  
  V <- t(apply(train_set, 1, vech_lt))          # T x m, row t = vech(Sigma_t)
  
  Y <- V[(p + 1):T_, , drop = FALSE]            # responses vech(Sigma_t)
  X <- matrix(1, nrow = T_ - p, ncol = 1)       # intercept
  for (j in 1:p) {
    X <- cbind(X, V[(p + 1 - j):(T_ - j), , drop = FALSE])   # lag j
  }
  
  coef <- qr.coef(qr(X), Y)                     # (1 + p*m) x m,  Y = X coef
  
  c_h  <- coef[1, ]
  e_I  <- vech_lt(diag(n))                      # vech(I_n)
  
  pars <- vector("list", p + 1)
  pars[[1]] <- spd_sqrt(unvech_lt(c_h, n))      # W_0
  for (j in 1:p) {
    B_j <- t(coef[1 + (j - 1) * m + (1:m), , drop = FALSE])  # m x m
    pars[[j + 1]] <- spd_sqrt(unvech_lt(B_j %*% e_I, n))    # W_j
  }
  pars
}





for (exp_id in 1:n_experiments){
  #generate the ground truth SPD parameters (DGP)
  set.seed(exp_id)
  Ws <- buildWeights(n, p) #W_0,W_1,...,W_p
  
  vec_gt = c()
  for (i in 1:(p + 1)) {
    vec_gt = c(vec_gt, vech(Ws[[i]]))
  }
  

  stopifnot(length(vec_gt) == N)

  DF <- data.frame(matrix(NA, nrow = N, ncol = 1 + n_cols_per_T * length(T_list)))
  
  
  series_means <- array(0, dim=c(length(T_list),N,n_methods)) 
  series_rmse <- array(0, dim=c(length(T_list),N,n_methods)) 
  series_covg <- array(0, dim=c(length(T_list),N,n_methods)) 
  
  
  t_max = T_list[length(T_list)]
  pgd_estimates <- array(
    NA_real_,
    dim = c(n_traj, length(T_list), N)
  )
  
  pgd_time <- array(NA_real_, dim = c(n_traj, length(T_list)))
  
  rm_estimates <- array(
    NA_real_,
    dim = c(n_traj, length(T_list), N)
  )
  
  var_estimates <- array(
    NA_real_,
    dim = c(n_traj, length(T_list), N)
  )
  
  
  set.seed(exp_id * n_traj)
  
  first_p_obs <- lapply(seq_len(p), function(i) {
    A <- matrix(rnorm(n * n), nrow = n, ncol = n)
    A <- (A+t(A))*0.5
    A <- higham_approx(A)
  })
  
  
  for (v in 1:n_traj){
    nsim = t_max+burn_in #nobs generated per 
    Y_final <- array(0,dim = c(n_traj, t_max, n, n))
    cat("trajectory: ", v, "\n")
    
  
    seed <- exp_id * n_traj + v
    set.seed(seed)
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
    
    train_set_full <- Yt[(burn_in+1):nsim,, , drop=FALSE]
    
    #discard burnin
    Y_final[v,,,] <- train_set_full
    
    
    for(t_idx in 1:length(T_list)){

      T_ = T_list[t_idx]
      train_set <- train_set_full[1:T_,,, drop=FALSE]
      
      
      
      ##PGD
      model <- spdmar(p, n)
      model_r <- model
      
      pgd_time[v, t_idx] <- system.time({
      for (epoch in 1:n_epochs){
        
        #compute the gradient
        euclidean_gradients <- grad.spdmar(model, train_set)
        
        #update the weights
        for (i in 1:(p + 1)) {
          G <- euclidean_gradients[[i]]
          W <- model$pars[[i]]
          #step_size <- l_rate * norm(W, "F") / norm(G, "F")
          step_size <- l_rate_pgd
          model$pars[[i]] <- model$pars[[i]] - step_size * euclidean_gradients[[i]]
          
          #projection via higham approx
          model$pars[[i]] <- higham_approx(model$pars[[i]])
          #model$pars[[i]] <- model$pars[[i]]/ norm(model$pars[[i]], "F")
        }
        
        loss <- frob_loss.spdmar(model, train_set)
        #cat("Epoch:", epoch, "Train Loss:", loss, "\n")
        
      }
      })[["elapsed"]]
      
      
      
      
      ##Riemannian optm
      for (epoch in 1:n_epochs){
        
        #compute the gradient
        euclidean_gradients <- grad.spdmar(model_r, train_set)
        
        #update the weights
        for (i in 1:(p + 1)) {
          G <- euclidean_gradients[[i]]
          W <- model_r$pars[[i]]
          
          sqrt_W <- matrix_sqrt(W)
          inv_sqrt_W <- matrix_inv_sqrt(W)
          
          riemannian_gradient <- W %*% ((G+t(G))/2) %*% W
          g_norm <- norm(riemannian_gradient, "F")
          
          #step_size <- l_rate * norm(W, "F") / g_norm
          step_size <- l_rate_rm
          P <- step_size * riemannian_gradient 
          
          
          
          A <- -1*inv_sqrt_W %*% P %*% inv_sqrt_W
          
          W_new <- sqrt_W %*%
            expm(A) %*%
            sqrt_W
          
          
          W_new <- (W_new+t(W_new))/2
          model_r$pars[[i]] <- W_new
          
          
        }
        
        loss <- frob_loss.spdmar(model_r, train_set)
        #cat("Epoch:", epoch, "Train Loss:", loss, "\n")
        
      }
      
      ##VAR on vech (closed form) + recovery of the SPD weights
      var_pars <- fit_var_vech(train_set, p)
      
      
      vec_pgd = c()
      vec_rm = c()
      vec_var = c()
      for (i in 1:(p + 1)) {
        
        vec_pgd = c(vec_pgd, vech(model$pars[[i]]))
        
        vec_rm = c(vec_rm, vech(model_r$pars[[i]]))
        
        vec_var = c(vec_var, vech(var_pars[[i]]))
        
      }
      
      pgd_estimates[v,t_idx, ] <- vec_pgd
      rm_estimates[v, t_idx, ] <- vec_rm
      var_estimates[v, t_idx, ] <- vec_var
      
      

      
      
    }
    
    
 
    
  }
  
  for(t_idx in 1:length(T_list)){
    
    T_ = T_list[t_idx]
    
    pgd_t <- pgd_estimates[, t_idx, , drop = FALSE]
    rm_t  <- rm_estimates[, t_idx, , drop = FALSE]
    var_t <- var_estimates[, t_idx, , drop = FALSE]
    
    dim(pgd_t) <- c(n_traj, N)
    dim(rm_t)  <- c(n_traj, N)
    dim(var_t) <- c(n_traj, N)
    
    gt_mat <- matrix(
      vec_gt,
      nrow = n_traj,
      ncol = N,
      byrow = TRUE
    )
    
    pgd_rmse <- sqrt(colMeans((pgd_t - gt_mat)^2))
    rm_rmse  <- sqrt(colMeans((rm_t  - gt_mat)^2))
    var_rmse <- sqrt(colMeans((var_t - gt_mat)^2))
    
    rm_means <- colMeans(rm_t)
    pgd_means <- colMeans(pgd_t)
    var_means <- colMeans(var_t)
    
    series_means[t_idx,,1] <- rm_means
    series_means[t_idx,,2] <- pgd_means
    series_means[t_idx,,3] <- var_means
    

    
    #pgd_se <- apply(pgd_estimates, 2, sd) / sqrt(n_traj)
    #rm_se  <- apply(rm_estimates, 2, sd) / sqrt(n_traj)
    
    #pgd_lower <- pgd_means - 1.96 * pgd_se
    #pgd_upper <- pgd_means + 1.96 * pgd_se
    
    #rm_lower <- rm_means - 1.96 * rm_se
    #rm_upper <- rm_means + 1.96 * rm_se
    
    
    
    
    series_rmse[t_idx,,1] <- rm_rmse
    series_rmse[t_idx,,2] <- pgd_rmse
    series_rmse[t_idx,,3] <- var_rmse
    
    series_covg[t_idx,,1]<- series_rmse[t_idx,,1] #percentage of replications in which the true parameter falls withing the estimate 95% confidence interval
    series_covg[t_idx,,2] <- series_rmse[t_idx,,1]
    series_covg[t_idx,,3] <- series_rmse[t_idx,,1]
    
  }


  

  
  
  #fill the CSV table

  #col names
  t_names = c()
  for(t_idx in 1:length(T_list)){
    T_ = T_list[t_idx]
    #first 3 cols are for the RM, remaining 3 cols for PGD
    t_names= c(t_names, 
               paste(T_, "hatq_ij RM", sep=" "),  paste(T_, "RMSE RM", sep=" "),  paste(T_, "COVG RM", sep=" "), 
               paste(T_, "hatq_ij PGD", sep=" "),  paste(T_, "RMSE PGD", sep=" "),  paste(T_, "COVG PGD", sep=" "),
               paste(T_, "hatq_ij VAR", sep=" "), paste(T_, "RMSE VAR", sep=" "), paste(T_, "COVG VAR", sep=" ")
               )
  }
  
  colnames(DF) <- c("q_ij", t_names)
  

  for (idx in 1:N){
    
    row_tb = c()
    for(t_idx in 1:length(T_list)){
      T_ = T_list[t_idx]
      #first 3 cols are for the RM, remaining 3 cols for PGD
      for (k in 1:n_methods) {
        row_tb = c(row_tb, series_means[t_idx,idx,k], series_rmse[t_idx,idx,k], series_covg[t_idx,idx,k])
      }
    }
    
    
    DF[idx, ] <- c(vec_gt[idx], row_tb)
  }
  
  
  write.csv(DF, "table_csv.csv")
  print("finished")

  print(pgd_time) #With n_traj = 1 and T_list = c(200, 500, 1000), pgd_time is a 1x3 matrix. Its entries are the PGD times (in seconds) on 200, 500 and 1000 observations.
  }

