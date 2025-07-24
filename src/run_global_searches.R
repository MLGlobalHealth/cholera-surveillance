#' Run global search on a pomp object
#' 
#' Starting points are passed to the function together with the search parameters.
#' @import dplyr, pomp
#' @return A pomp object
#' 

run_global_search <- function(
    model,
    guess_points,
    rw.sd,
    Np,
    Nmif,
    Nreps,
    Nreps_eval,
    out_filename
) {
  start <- Sys.time()
  model_params <- coef(model)
  fixed_params <- model_params[!(names(model_params) %in% names(guess_points))]
  
  
  m_global <- foreach (guess=iter(guess_points, "row"), .combine=c) %dopar% {
    (
      model
      %>% mif2(
        Np=Np,
        Nmif=Nmif, 
        params=c(fixed_params, unlist(guess)),
        cooling.fraction.50 = 0.5,
        rw.sd = rw.sd
      )
      %>% mif2(cooling.fraction.50 = 0.1)
    )
  }
  
  loglik_global <- (
    foreach(
      i=1:Nreps,
      .combine=rbind
    ) %dopar% {
      logmeanexp(
        replicate(
          Nreps_eval,
          logLik(pfilter(model, params=coef(m_global[[i]]), Np=1000))
        ),
        se=TRUE
      )
    }
  )
  
  end <- Sys.time()
  cat("Run time of global iterated particle filter: ", as.integer(round(difftime(end, start, units='mins'))), " minutes")
}


#' Collect and save results of global search
#' 
#' @import dplyr, pomp
#' @return A pomp object
#' 

collect_and_save_loglik <- function(
    model_search,
    loglik_values,
    out_filename,
    reset_data = FALSE
) {
  
  df_results <- (
    data.frame(
      t(sapply(model_search, coef)),
      log_lik=loglik_values[,1],
      log_lik_se=loglik_values[,2]
    )
  )
  
  if (reset_data) {
    write.table(
      df_results, 
      file = out_filename, 
      append = FALSE,
      col.names = TRUE,
      row.names = FALSE
    )
  } else {
    append <- FALSE
    col.names <- TRUE
    row.names <- FALSE
    if (file.exists(out_filename)) {
      append <- TRUE
      col.names <- FALSE
    }
    (
      write.table(
        df_results, 
        file = out_filename, 
        append = append,
        col.names = col.names,
        row.names = row.names
      )
    )
  }
  df_results
}

#' Simulate from best model
#' 
#' @import dplyr, pomp
#' @return A pomp object
#' 
#' 

simulate_best_model <- function(
    model_search,
    loglik_values,
    nsim
) {
  sims <- (
    model_search[[which.max(loglik_values[,1])]]
    %>% simulate(
      nsim = nsim,
      format = "data.frame",
      include.data = FALSE
    )
  )
  sims
}


#' Aggregate simulation results
#' 
calculate_sim_stats <- function(sims) {
  sims_stats <- (
    sims
    %>% filter(.id != "data")
    %>% select(.id, date, y)
    %>% group_by(date)
    %>% summarise(
      q05 = quantile(y, 0.025, na.rm = TRUE),
      q50 = quantile(y, 0.5, na.rm = TRUE),
      mean = mean(y, na.rm = TRUE),
      q95 = quantile(y, 0.975, na.rm = TRUE),
    )
    %>% ungroup()
  )
  sims_stats
}