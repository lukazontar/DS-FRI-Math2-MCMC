#### Assignment: 
####  Implement Metropolis-Hastings with Multivariate Normal proposal 
####  (mean 0, covariance matrix is a parameter).

# p - target distribution
# m - number of samples
# x_0 - starting state in MCMC algorithm
# cov_mtx - multivariate normal covariance matrix
# ix_chain - index of chain
metropolis_hastings_mv_norm <- function(p, m, x_0, cov_mtx, ix_chain) {
  # Start measuring time 
  start_time <- Sys.time()
  
  # Repeat 0s per dimension
  mu <- rep(0, sqrt(length(cov_mtx)))
  
  
  # Initialize current state as starting state
  x_i <- x_0
  
  # Samples list
  x <-  data.frame()
  
  # Generate m (possibly dependent) samples from density p
  for (i in 1:m) {
    repeat {
      x_star <- x_i + mvrnorm(n=1, mu=mu, Sigma=cov_mtx)
      alpha <- min(
        1, 
        p(x_star)/p(x_i)
      )
      sample_uniform <- runif(n=1, min=0, max=1)
      
      # If condition met, accept sample, else reject
      if (sample_uniform <= alpha) {
        x_i <- x_star
        x <- rbind(x, data.frame(matrix(data=x_i, nrow=1)))
        break()
      }
    }
  }
  
  # Add index and chain data to x
  x$Chain <- as.character(rep(ix_chain, m))
  x$Index <- 1:m
  
  # End measuring time 
  total_time_spent <- difftime(Sys.time(), start_time, units="secs")
  
  return_value <- c()
  return_value$x = x
  return_value$total_time_spent = total_time_spent
  
  # Return independent samples list x with  with total time spent
  return(return_value)
}