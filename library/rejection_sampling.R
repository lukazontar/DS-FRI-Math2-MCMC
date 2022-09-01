#### Assignment: 
####  Implement rejection sampling with your own choice of 2D envelope. 
####  This needn't be done for (4), because the dimensionality is too high.

# f - multiplier of target density function - f(x) = C * p(x)
# g_sample - 2D envelope
# g_density - calculate density of proposal density function 
# M - positive constant such that f(x) <= M*g(x)
# m - number of samples
# ix_chain - Index of chain
rejection_sampling <- function(f, g_sample, g_density, M, m, ix_chain) {
  # Start measuring time 
  start_time <- Sys.time()
  # Samples list
  x <-  data.frame()
  
  # Generate m independent samples from density p
  for (i in 1:m) {
    repeat {
      # Sample y from g
      y <- g_sample()
      
      # Sample u from Unif(0, 1)
      u <- runif(n=1, min=0, max=M*g_density(y))  
      if (u > f(y)) {
        break()
      }
    }
    x <- rbind(x, data.frame(matrix(data=y, nrow=1)))
    
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
