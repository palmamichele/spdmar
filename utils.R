build_matrix_from_vech <- function(data){
  #convert dataframe to numeric vector and return symmetric matrix from it
  x <- as.numeric(data)
  
  m <- length(x)
  
  n <- (-1+sqrt(1+8*m))/2
  n <- as.integer(n)
  
  M <- matrix(0, n, n)
  
  for (j in 1:n) {
    
    start <- 1 + (j - 1) * n - (j - 1) * j / 2
    end <- start + n - j
    
    els <- x[start:end]
    
    M[j, j:n] <- els
    M[j:n, j] <- els
  }
  
  return (M)
}


vec <- function(matrix){
  n <- nrow(matrix)
  m <- ncol(matrix)
  
  v <- replicate(n*m, 0)
  for (i in 1:n){
    start <- (i-1)*m+1
    end <- i*m
    v[start:end] <- matrix[i,]
  }
  return (v)
} 

# vech <- function(matrix){
#   
# }


recover_from_kronecker <- function(C){
  n <- sqrt(nrow(C))
  
  W <- matrix(0,n,n) 
  
  i<-1
  # for (i in 1:n){
  start <- i*n+1
  end <- (i+1)*n
  C_i_block <- C[start:end, start:end]
    
    # if (sum(diag(C_i_block)!=0){
  W <- C_i_block*sqrt(sum(diag(C)))/sum(diag(C_i_block))
  return (W)
    # }
    
  # }
  
}



higham_approx <- function(A, eps=1e-6){
  #symmetric matrix A 
  
  #numerical stab
  A = (A + t(A)) / 2 
  
  eig <- eigen(A, symmetric=TRUE)
  eigvals<- eig$values
  eigvecs <- eig$vectors
  
  #reeig
  eigvals <- pmax(eigvals, eps)
  X <- eigvecs %*% diag(eigvals) %*% t(eigvecs)
  
  X<- (X + t(X)) / 2
  return (X)
}



#example usage
W = matrix(data=c(1,2,3,4),nrow=2,ncol=2, byrow=TRUE)
C = kronecker(W,W)
W_prime = recover_from_kronecker(C)



