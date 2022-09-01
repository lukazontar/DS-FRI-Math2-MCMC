#### Assignment: 
####  Implement HMC (step-size and number of steps are parameters; 
####  you can use unit mass matrix, but if you want to impress, 
####  you can tune the diagonal elements as well)

# minus_log_f_grad - gradient of function that is proportional to target density
# q_start - starting value of position
# p_start - starting value of momentum
# L - number of steps
# epsilon - step size
run_leapfrog_steps <- function(minus_log_f_grad, q_start, p_start, L, epsilon) {
  # Initialize momentum and position
  q_star <- q_start
  p_star <- p_start
  
  
  # Half step of momentum
  p_star <- p_star - minus_log_f_grad(q_star) * epsilon / 2
  
  # L full steps step of momentum and position
  for(i in 1:L) {
    # Full step - position
    q_star <- q_star + epsilon * p_star
    
    # Full step - momentum - we make one less momentum step, because we make 2 halves before and after additionally! 
    if(i < L) {
      p_star <- p_star - minus_log_f_grad(q_star) * epsilon
    }
  }
  
  # Half step of momentum
  p_star <- p_star - epsilon/2 * minus_log_f_grad(q_star)
  
  return_value <- c()
  return_value$p_star <- p_star
  return_value$q_star <- q_star
  
  # Return new p* and q*
  return(return_value)
}

# minus_log_f - function proportional to our target density
# minus_log_f_grad - gradient of -log(f(x))
# q_0 - starting value
# epsilon - step size (strictly larger than 0)
# L - number of steps
# m - number of samples
hamiltonian_monte_carlo <- function(minus_log_f, minus_log_f_grad, q_0, epsilon, L, m, ix_chain) {
  # Start measuring time 
  start_time <- Sys.time()
  # Samples list
  x <-  data.frame()
  
  q_i <- q_0
  
  # Generate m (possibly dependent) samples from density p
  for (i in 1:m) {
    # Sample new momentum
    p <- rnorm(n=length(q_i), mean=0, sd=1)
    
    # Run L Leapfrog steps with size epsilon from (q_i, p) to get p_star and q_star
    result <- run_leapfrog_steps(minus_log_f_grad=minus_log_f_grad, q_start=q_i, p_start=p, L=L, epsilon=epsilon)
    p_star <- result$p_star
    q_star <- result$q_star
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    p_star=-p_star
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U = minus_log_f(q_i)
    current_K = sum(p^2) / 2
    proposed_U = minus_log_f(q_star)
    proposed_K = sum(p_star^2) / 2
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
    {
      # Accept sample
      q_i <- q_star
    }
    x <- rbind(x, data.frame(matrix(data=q_i, nrow=1)))
    
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