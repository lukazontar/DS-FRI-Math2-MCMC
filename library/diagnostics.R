# variables - vector
# samples - vector of generated samples
# total_execution_time - total execution time in seconds
standard_mcmc_diagnostics <-
  function(variables, samples, total_execution_time, ground_truth) {
    ess_list <- c()
    ess_per_sec_list <- c()
    estimate_mean_list <- c()
    estimate_error_list <- c()
    
    plots <- list()
    
    ix <- 1
    iy <- 1
    for (variable in variables) {
      # ggplot trace plot
      g_traceplot <-
        ggplot(data = samples) + geom_line(aes_string(x = "Index", y = variable, color =
                                                        "Chain")) + ggtitle(paste("Sampling trace plot - ", variable, sep = "")) + theme_bw()
      
      plots[[ix]] <- g_traceplot
      
      # ggplot autocorrelation plot
      lag_max = 100
      autocov_data <-
        data.frame(Correlation = as.vector(acf(
          samples[variable], lag.max = lag_max, plot = FALSE
        )$acf),
        Lag = 1:(lag_max + 1))
      g_autocorrelation <-
        ggplot(data = autocov_data) + geom_bar(aes(x = Lag, y = Correlation), stat =
                                                 'identity') + ggtitle(paste("Lag-k autocorrelation plot - ", variable, sep =
                                                                               "")) + theme_bw()
      
      plots[[ix + 1]] <- g_autocorrelation
      
      ## Overall statistics
      # ESS
      ess_variable <- ess(samples[variable])
      ess_list <- c(ess_list, ess_variable)
      
      # ESS / second
      ess_variable_per_sec <- ess_variable / total_execution_time
      ess_per_sec_list <- c(ess_per_sec_list, ess_variable_per_sec)
      
      # Estimate mean, error
      estimate_mcse <- mcse(as.vector(unlist((
        sqrt((samples[variable] - ground_truth[iy])^2)
        ))))
      estimate_mean_list <- c(estimate_mean_list, estimate_mcse$est)
      estimate_error_list <- c(estimate_error_list, estimate_mcse$se)
      
      iy <- iy + 1
      ix <- ix + 2
    }
    df_stats <-
      data.frame(
        Variable = variables,
        ESS = ess_list,
        ESS_per_Second = ess_per_sec_list,
        Est_Mean = estimate_mean_list,
        Est_StdErr = estimate_error_list
      )
    
    return_value <- c()
    return_value$plots <- plots
    return_value$stats <- df_stats
    
    return(return_value)
  }

# variables - list of parameters we are trying to estimate
# hmc_results - results from Hamiltonian Monte Carlo algorithm
# rej_samp_results - results from rejection sampling algorithm
# mh_results - results from Metropolis-Hastings algorithm
# ground_truth - ground_truth estimates
summarize_results <-
  function(variables,
           hmc_results,
           rej_samp_results,
           mh_results, include_rejection_sampling=TRUE) {
    annotated_plots <- list()
    mh_plot <-
      ggarrange(
        plotlist = mh_results$plots,
        ncol = 2,
        nrow = length(variables),
        common.legend = TRUE,
        legend = "bottom"
      )
    
    annotated_plots[[1]] <-
      annotate_figure(p = mh_plot,
                      top = text_grob("Metropolis Hastings", face = "bold", size = 14))
    
    hmc_plot <-
      ggarrange(
        plotlist = hmc_results$plots,
        ncol = 2,
        nrow = length(variables),
        common.legend = TRUE,
        legend = "bottom"
      )
    
    annotated_plots[[2]] <-
      annotate_figure(p = hmc_plot,
                      top = text_grob("Hamiltonian Monte Carlo", face = "bold", size = 14))

    if (include_rejection_sampling) {
      rej_samp_plot <-
        ggarrange(
          plotlist = rej_samp_results$plots,
          ncol = 2,
          nrow = length(variables),
          common.legend = TRUE,
          legend = "bottom"
        )
      
      annotated_plots[[3]] <-
        annotate_figure(p = rej_samp_plot,
                        top = text_grob("Rejection sampling", face = "bold", size = 14))
      
    plot <-
      ggarrange(
        plotlist = annotated_plots,
        ncol = 1,
        nrow = 3,
        common.legend = TRUE,
        legend = "bottom"
      )
        
    } else {
      plot <-
        ggarrange(
          plotlist = annotated_plots,
          ncol = 2,
          nrow = 1,
          common.legend = TRUE,
          legend = "bottom"
        )
    }
  
    
    hmc_full_string <- "Hamiltonian Monte Carlo"
    mh_full_string <- "Metropolis Hastings"
    rej_samp_full_string <- "Rejection sampling"
    for (i in 1:length(variables)) {
      hmc_var_string <-
        sprintf(
          "%s: ESS=%#.2f, ESS/s=%#.2f, mean(est)=%#.2f, se(est)=%#.2f",
          hmc_results$stats[i, "Variable"],
          hmc_results$stats[i, "ESS"],
          hmc_results$stats[i, "ESS_per_Second"],
          hmc_results$stats[i, "Est_Mean"],
          hmc_results$stats[i, "Est_StdErr"]
        )
      mh_var_string <-
        sprintf(
          "%s: ESS=%#.2f, ESS/s=%#.2f, mean(est)=%#.2f, se(est)=%#.2f",
          mh_results$stats[i, "Variable"],
          mh_results$stats[i, "ESS"],
          mh_results$stats[i, "ESS_per_Second"],
          mh_results$stats[i, "Est_Mean"],
          mh_results$stats[i, "Est_StdErr"]
        )
      if (include_rejection_sampling) {
      rej_samp_var_string <-
        sprintf(
          "%s: ESS=%#.2f, ESS/s=%#.2f, mean(est)=%#.2f, se(est)=%#.2f",
          rej_samp_results$stats[i, "Variable"],
          rej_samp_results$stats[i, "ESS"],
          rej_samp_results$stats[i, "ESS_per_Second"],
          rej_samp_results$stats[i, "Est_Mean"],
          rej_samp_results$stats[i, "Est_StdErr"]
        )
        rej_samp_full_string <-
          paste(rej_samp_full_string, rej_samp_var_string, sep = "\n")
      }
      
      hmc_full_string <-
        paste(hmc_full_string, hmc_var_string, sep = "\n")
      mh_full_string <- paste(mh_full_string, mh_var_string, sep = "\n")
    
    }
    
    stats <-
      paste(hmc_full_string, mh_full_string, rej_samp_full_string, sep = "\n")
    
    return_value <- c()
    return_value$summary_plot <- plot
    return_value$summary_stats <- stats
    return(return_value)
  }

generate_results <-
  function(minus_log_f,
           minus_log_f_grad,
           starting_point,
           epsilon,
           L,
           m,
           f,
           g_sample,
           g_density,
           M,
           p,
           cov_mtx,
           ground_truth, 
           include_rejection_sampling=TRUE) {
    # Initialize global variables for Rejection sampling
    rej_samp_samples <- data.frame()
    rej_samp_total_execution_times <- c()
    
    # Initialize global variables for Metropolis Hastings
    mh_samples <- data.frame()
    mh_total_execution_times <- c()
    
    hmc_samples <- data.frame()
    hmc_total_execution_times <- c()
    
    for (i in 1:5) {
      print(i)
      # Hamiltonian Monte Carlo - Scenario 1
      hmc_result_single_chain <-
        hamiltonian_monte_carlo(
          minus_log_f = minus_log_f,
          minus_log_f_grad = minus_log_f_grad,
          q_0 = starting_point,
          epsilon = epsilon,
          L = L,
          m = m,
          ix_chain = i
        )
      hmc_samples <- rbind(hmc_samples, hmc_result_single_chain$x)
      hmc_total_execution_times <-
        c(hmc_total_execution_times,
          hmc_result_single_chain$total_time_spent)
      print("HMC done")
      if (include_rejection_sampling) {
        # Rejection sampling - Scenario 1
        rej_samp_result_single_chain <-
          rejection_sampling(
            f = f,
            g_sample = g_sample,
            g_density = g_density,
            M = M,
            m = m,
            ix_chain = i
          )
        rej_samp_samples <-
          rbind(rej_samp_samples, rej_samp_result_single_chain$x)
        rej_samp_total_execution_times <-
          c(
            rej_samp_total_execution_times,
            rej_samp_result_single_chain$total_time_spent
          )
        print("Rejection sampling done")  
      }
      
      # Metropolis Hastings - Scenario 1
      mh_result_single_chain <-
        metropolis_hastings_mv_norm(
          p = p,
          m = m,
          x_0 = starting_point,
          cov_mtx = cov_mtx,
          ix_chain = i
        )
      mh_samples <- rbind(mh_samples, mh_result_single_chain$x)
      mh_total_execution_times <-
        c(mh_total_execution_times,
          mh_result_single_chain$total_time_spent)
      print("MH done")
    }

    
    # Evaluating for Hamiltonian Monte Carlo
    hmc_total_execution_time_mean <- mean(hmc_total_execution_times)
    hmc_results <-
      standard_mcmc_diagnostics(
        variables = variables,
        samples = hmc_samples,
        total_execution_time = hmc_total_execution_time_mean,
        ground_truth = ground_truth
      )
    
    if (include_rejection_sampling) {
      # Evaluating for Rejection sampling
      rej_samp_total_execution_time_mean <-
        mean(rej_samp_total_execution_times)
      rej_samp_results <-
        standard_mcmc_diagnostics(
          variables = variables,
          samples = rej_samp_samples,
          total_execution_time = rej_samp_total_execution_time_mean,
          ground_truth = ground_truth
        )
    }
    
    # Evaluating for Metropolis Hastings
    mh_total_execution_time_mean <- mean(mh_total_execution_times)
    mh_results <-
      standard_mcmc_diagnostics(
        variables = variables,
        samples = mh_samples,
        total_execution_time = mh_total_execution_time_mean,
        ground_truth = ground_truth
      )
    
    return_value <- c()
    return_value$hmc_total_execution_time_mean = hmc_total_execution_time_mean
    return_value$hmc_results = hmc_results
    if (include_rejection_sampling) {
      return_value$rej_samp_total_execution_time_mean = rej_samp_total_execution_time_mean
      return_value$rej_samp_results = rej_samp_results
    }
    return_value$mh_total_execution_time_mean = mh_total_execution_time_mean
    return_value$mh_results = mh_results
    return(return_value)
  }
