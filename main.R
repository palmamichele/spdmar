#install.packages("vars")   install.package("expm")

###ARGS
data_path="~/Desktop/Research Projects/SPDMAR/code/spdmar/data/adjusted_RCOV50_COVDJkernel.csv" 
p=3 #order of autoregressive models
set.seed(123)
n_epochs = 100
l_rate = 0.009
train_thr = 0.8 
###

source("~/Desktop/Research Projects/SPDMAR/code/spdmar/utils.R")
source("~/Desktop/Research Projects/SPDMAR/code/spdmar/spdmar.R")
library(vars)
library(expm)

vect_data <- read.csv(data_path, sep = ",")
vect_data$X <- NULL
T_or = nrow(vect_data)
N = nrow(build_matrix_from_vech(vect_data[1,]))
time_series_matrix <- array(0, dim=c(T_or,N,N))
time_series_vec <- array(0, dim=c(T_or, N*N))


#build time series of spd matrice, t=1 index of least recent obs, T_or the index for most recent obs
for (t_ in 1:T_or){
  time_series_matrix[t_,,] <- build_matrix_from_vech(vect_data[t_,]) 
  time_series_vec[t_,] <- vec(time_series_matrix[t_,,])
}


train_last_idx = floor(train_thr * T_or)
train_set = time_series_matrix[1:train_last_idx,,,drop = FALSE]
test_set = time_series_matrix[(train_last_idx+1):T_or,,,drop = FALSE]

cat("number of obs: ", T_or, "\n")
cat("number of training obs: ", train_last_idx, "\n")

#PGD update
model <- spdmar(p, N)
print(model)
cat("Train loss before training ", frob_loss.spdmar(model,time_series_matrix), "\n")

#train test split? rolling-econometric-window?


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


cum_test_loss = 0
cum_rw_loss = 0
for (t_ in (p+1):(nrow(test_set))){
  #build lagged series
  input_series = test_set[(t_-p):(t_-1),,,drop = FALSE]
  target = test_set[t_,,]
  rwalk = test_set[t_-1,,]
  prediction = predict.spdmar(model, input_series)
  E = prediction - target
  frob_ob = sum(E^2)
  cum_test_loss = cum_test_loss + frob_ob
  cum_rw_loss = cum_rw_loss + sum( (target - rwalk)^2)
}

cat("Model Test loss: ", cum_test_loss, "\n")
cat("RW Test loss: ", cum_rw_loss, "\n")


#Riemannian optm
model_r <- spdmar(p, N)
print(model_r)

#train test split? rolling-econometric-window?


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


cum_test_loss = 0
cum_rw_loss = 0
for (t_ in (p+1):(nrow(test_set))){
  #build lagged series
  input_series = test_set[(t_-p):(t_-1),,,drop = FALSE]
  target = test_set[t_,,]
  rwalk = test_set[t_-1,,]
  prediction = predict.spdmar(model_r, input_series)
  E = prediction - target
  frob_ob = sum(E^2)
  cum_test_loss = cum_test_loss + frob_ob
  cum_rw_loss = cum_rw_loss + sum( (target - rwalk)^2)
}

cat("Riemannian Model Test loss: ", cum_test_loss, "\n")
cat("RW Test loss: ", cum_rw_loss, "\n")

#p=1 RW test loss 148320.8 , model test loss : 122230.1
#p=5       147146.2       , 644356.4 

#estimate VAR on time_series_vec
# Y <- t(time_series_vec)
# var_fit <- VAR(
#   Y,
#   p = order,
#   type = "const"
# )
# 
# summary(var_fit)
# 
# A <- Acoef(var_fit)

#get W