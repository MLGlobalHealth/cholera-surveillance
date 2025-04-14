#' Builder for pomp models of the Uvira, DRC cholera data
#' 
#' Generate a pomp object for fitting the data with several levels of complexity. 
#' The simplest is a SEIAR model with a constant force of infection. The following features can be added:
#' 1. Switch between suspected cholera to detected cholera cases
#' 2. Time-varying force of infection:
#'   - Using pre-calculated combined covariate
#'   - Adding a linear regression model for the top n covariates
#' 3. Separate compartments for people recoverying from an asymptomatic and a symptomatic infection
#' 4. A hospitalized compartment that is a fraction of the infected compartment and includes people who are not infective because they are isolated
#' 5. Influx and outflux of internally displaced people
#' 6. Effect of natural disasters
#' 7. Two compartments for people vaccinated during the vaccination campaigns
#' 
#' 
#' @import dplyr, pomp
#' @return A pomp object

build_pomp_model <- function(
    data,
    use_confirmed = FALSE,
    use_cov = FALSE,
    use_cov_regression = FALSE,
    scale_infectious = FALSE,
    scale_hospitalized = FALSE,
    use_separate_recovery = FALSE,
    add_hospitalized = FALSE,
    add_idps = FALSE,
    add_disasters = FALSE,
    add_vaccination = FALSE
) {
  
  data_cols <- names(data)
  
  # Data to model
  if (use_confirmed) {
    y = c("rdt_confirmed")
  } else {
    y = c("reports")
  }
  df_model <- (
    data
    %>% select(
      date, 
      y = y
    )
    %>% filter(date >= 0)
  )
  
  # Create covariate table
  covar_names = c("birthrate", "deathrate")
  if (use_cov) {
    if (use_cov_regression) {
      cov_reg_names <- data_cols[grepl("cov_", data_cols)]
      covar_names <-  c(covar_names, cov_reg_names)
    } else {
      covar_names <- c(covar_names, "cov")
    }
  }
  if (add_idps) {
    covar_names <- c(covar_names, "idps")
  }
  if (add_disasters) {
    covar_names <- c(covar_names, "dis")
  }
  if (add_vaccination) {
    covar_names <- c(covar_names, "vacc_1_prop", "vacc_2_prop")
  }
  
  df_covar <- data %>% select(date, all_of(covar_names))
  covar <- covariate_table(
    df_covar,
    times = "date"
  )
  
  # Observations
  obs_names <- c("y")
  
  # Accumulation variables
  accum_vars <- c("C", "W")
  
  # States
  states_names <- c("S", "E", "A", "I", "C", "W", "N")
  ivp_names <- c("S", "E", "A", "I")
  ivp_values <- c("s_0", "e_0", "a_0", "i_0")
  if (add_hospitalized) {
    states_names <- c(states_names, "H")
    ivp_names <- c(ivp_names, "H")
    ivp_values <- c(ivp_values, "h_0")
  }
  if (use_separate_recovery) {
    states_names <- c(states_names, "RA", "RI")
    ivp_names <- c(ivp_names, "RA", "RI")
    ivp_values <- c(ivp_values, "ra_0", "ri_0")
  } else {
    states_names <- c(states_names, "R")
    ivp_names <- c(ivp_names, "R")
    ivp_values <- c(ivp_values, "r_0")
  }
  if (add_vaccination) {
    states_names <- c(states_names, "V1", "V2")
  }
  
  # Transformation groups
  
  rp_log_names <- c(
    "sigmaSE",
    "tau"
  )
  rp_logit_names <- c(
    "f_ai",
    "eps_A",
    "rho"
  )
  
  if (use_cov) {
    rp_log_names <- c(rp_log_names, "a_cov")
    if (use_cov_regression) {
      for (i in 1:length(cov_reg_names)) {
        rp_log_names <- c(rp_log_names, paste0("b_cov_", i))
      }
    } else {
      rp_log_names <- c(rp_log_names, "b_cov")
    }
  } else {
    rp_log_names <- c(rp_log_names, "Beta")
  }
  
  if (add_hospitalized) {
    rp_logit_names <- c(rp_logit_names, "f_h")
  }
  if (add_idps) {
    rp_logit_names <- c(rp_logit_names, "b_idps")
  }
  if (add_disasters) {
    rp_log_names <- c(rp_log_names, "b_dis")
  }
  if (add_vaccination) {
    rp_log_names <- c(rp_log_names, "mu_vs1", "mu_vs2")
    rp_logit_names <- c(rp_logit_names, "f_v1", "f_v2")
  }
  if (use_separate_recovery) {
    rp_log_names <- c(rp_log_names, "mu_ras", "mu_ris")
  } else {
    rp_log_names <- c(rp_log_names, "mu_rs")
  }
  
  # Fixed parameters 
  N_0 <- 243763
  mu_latend <- 1 / ((5/7) / 52) # inverse of latent period (5 days transformed in weeks)
  mu_ar <- 1 / ((10/7) / 52) # inverse of asymptomatic period (10 days transformed in weeks)
  mu_ir <- 1 / ((10/7) / 52)  # inverse of infectious period (10 days transformed in weeks)
  
  fixed_params <- c(
    N_0 = N_0,
    mu_latend = mu_latend,
    mu_ar = mu_ar,
    mu_ir = mu_ir
  )
  if (add_hospitalized) {
    mu_hr <- 1 / ((10/7) / 52)  # inverse of infectious period (10 days transformed in weeks)
    fixed_params <- c(fixed_params, mu_hr = mu_hr)
  }
  
  par_trans <- parameter_trans(
    log = rp_log_names,
    logit = rp_logit_names,
    barycentric = ivp_values
  )
  
  param_names <- c(rp_log_names, rp_logit_names, ivp_values, names(fixed_params))
  
  ##################################
  # Model functions
  
  # Build observation functions
  dmeas <- Csnippet("
    double f;
    if (tau > 0.0) {
      f = dnbinom_mu(nearbyint(y), 1.0/tau, rho*C, give_log);
    }
    else {
      f = dpois(nearbyint(y), rho*C, give_log);
    }
    lik = (give_log) ? f : exp(f);
  ")
  
  rmeas <- Csnippet("
    if (tau > 0.0) {
      y = rnbinom_mu(1.0/tau, rho*C);
    }
    else {
      y = rpois(rho*C);
    }
  ")
  
  # Build rinit
  scaling_factor <- paste0("double m = N_0/(", paste(ivp_values, collapse = "+"), ");\n")
  rinit_str <- paste0(
    scaling_factor,
    paste(paste0(ivp_names, " = nearbyint(m*", ivp_values, ");"), collapse = "\n"),
    "\n",
    if (add_vaccination) "V1 = 0;\nV2 = 0;\n",
    "C = 0;\nW = 0;\n",
    paste0("N = ", paste(ivp_names, collapse = "+"), ";")
  )
  rinit <- Csnippet(rinit_str)
  
  # Build rprocess
  
  ## Parameters
  num_rates <- 11
  if (add_hospitalized) {
    num_rates <- num_rates + 3
  } 
  if (use_separate_recovery) {
    num_rates <- num_rates + 2
  }
  if (add_vaccination) {
    num_rates <- num_rates + 6
  }
  
  num_idps_rates <- 4
  if (add_hospitalized) {
    num_idps_rates <- num_idps_rates + 1
  }
  if (add_vaccination) {
    num_idps_rates <- num_idps_rates + 2
  }
  
  ## Universal rates
  
  births <- "births = rpois(birthrate * N * dt);\n\n"
  noise <- "dw = rgammawn(sigmaSE, dt);\n\n"
  
  variable_definitions <- paste0(
    "double foi, births, dw;\n",
    if (add_hospitalized) "double hosp;\n",
    paste0("double rate[", num_rates, "];\n"),
    paste0("double trans[", num_rates, "];\n"),
    if (add_idps) {
      paste0(
        "double idps_scaled, s;\n",
        paste0("double rate_idps[", num_idps_rates, "];\n"),
        paste0("int trans_idps[", num_idps_rates, "];\n")
      )
    },
    if (add_vaccination) {
      if (use_separate_recovery) {
        "int SV1, RAV1, RIV1, AV1, V1V2;\n\n"
      } else {
        "int SV1, RV1, AV1, V1V2;\n\n"
      }
    },
    "\n"
  )
  
  ## IDPS
  
  if (add_hospitalized) {
    ref_comp <- "H"
    
  } else {
    ref_comp <- "I"
  }
  ref_pop <- paste0("(N - ", ref_comp, ")")
  no_trans_comp <- c("C", "W", "N")
  trans_comp <- states_names[!(states_names %in% no_trans_comp)]
  out_comp <- states_names[!(states_names %in% c(ref_comp, no_trans_comp))]
  vacc_comp <- states_names[startsWith(states_names, "V")]
  in_comp <- states_names[!(states_names %in% c(ref_comp, no_trans_comp, vacc_comp))]
  idps <- paste0(
    "if (idps < 0) {\n",
    paste0("  if (N > ", ref_comp, ") {\n"),
    paste0("    idps_scaled = floor(fmin((-idps * b_idps), ", ref_pop, "));\n"),
    paste0("    s = ", ref_pop, " / (", paste(out_comp, collapse = "+"), ");\n"),
    paste0(paste0("    rate_idps[", seq(0, num_idps_rates - 1, 1), "] = s*", out_comp, "/", ref_pop, ";", collapse = "\n")),
    paste0("\n    rmultinom(idps_scaled, &rate_idps[0], ", num_idps_rates, ", &trans_idps[0]);\n    "),
    paste0(out_comp, " -= fmin(trans_idps[", seq(0, num_idps_rates - 1, 1), "], ", out_comp, ");", collapse = "\n    "),
    "\n  }\n",
    "}\n",
    "if (idps > 0) {\n",
    paste0("  idps_scaled = floor(fmin((idps * b_idps), ", ref_pop, "));\n"),
    paste0(paste0("  ", in_comp, " += nearbyint(idps_scaled/", length(in_comp), ");"), collapse = "\n"),
    "\n}\n",
    paste0("\nN = ", paste(trans_comp, collapse = "+"), ";\n\n")
  )
  
  ## Vaccination
  
  if (use_separate_recovery) {
    vacc <- paste0(
      "if (vacc_1_prop > 0) {\n",
      "  SV1 = rbinom(S, vacc_1_prop);\n",
      "  RAV1 = rbinom(RA, vacc_1_prop);\n",
      "  RIV1 = rbinom(RI, vacc_1_prop);\n",
      "  AV1 = rbinom(A, vacc_1_prop);\n",
      "  S -= SV1;\n",
      "  RA -= RAV1;\n",
      "  RI -= RIV1;\n",
      "  A -= AV1;\n",
      "  V1 += SV1 + RAV1 + RIV1 + AV1;\n",
      "}\n",
      "if (vacc_2_prop > 0) {\n",
      "  SV1 = rbinom(S, vacc_2_prop);\n",
      "  RAV1 = rbinom(RA, vacc_2_prop);\n",
      "  RIV1 = rbinom(RI, vacc_2_prop);\n",
      "  AV1 = rbinom(A, vacc_2_prop);\n",
      "  V1V2 = rbinom(V1, vacc_2_prop);\n",
      "  S -= SV1;\n",
      "  RA -= RAV1;\n",
      "  RI -= RIV1;\n",
      "  A -= AV1;\n",
      "  V1 += SV1 + RV1 + AV1 - V1V2;\n",
      "  V2 += V1V2;\n",
      "}\n"
    )
  } else {
    vacc <- paste0(
      "if (vacc_1_prop > 0) {\n",
      "  SV1 = rbinom(S, vacc_1_prop);\n",
      "  RV1 = rbinom(R, vacc_1_prop);\n",
      "  AV1 = rbinom(A, vacc_1_prop);\n",
      "  S -= SV1;\n",
      "  R -= RV1;\n",
      "  A -= AV1;\n",
      "  V1 += SV1 + RV1 + AV1;\n",
      "}\n",
      "if (vacc_2_prop > 0) {\n",
      "  SV1 = rbinom(S, vacc_2_prop);\n",
      "  RV1 = rbinom(R, vacc_2_prop);\n",
      "  AV1 = rbinom(A, vacc_2_prop);\n",
      "  V1V2 = rbinom(V1, vacc_2_prop);\n",
      "  S -= SV1;\n",
      "  R -= RV1;\n",
      "  A -= AV1;\n",
      "  V1 += SV1 + RV1 + AV1 - V1V2;\n",
      "  V2 += V1V2;\n",
      "}\n\n"
    )
  }
  
  ## Force of infection
  if (use_cov) {
    if (use_cov_regression) {
      Beta <- paste0(
        "(a_cov + ",
        paste0("b_cov_", seq(1, length(cov_reg_names), 1), " * ", cov_reg_names, collapse = " + ")
      )
    } else {
      Beta <- "(a_cov + b_cov * cov"
    }
    if (add_disasters) {
      Beta <- paste0(Beta, " + b_dis * dis)")
    } else {
      Beta <- paste0(Beta, ")")
    }
  } else {
    Beta <- "Beta"
  }
  inf <- "(I + eps_A * A)"
  if (scale_infectious) { 
    inf <- paste0("pow(", inf, ", alpha)")
  }
  foi <- paste0(
    "foi = ",
    Beta,
    "*",
    inf,
    "/N;\n\n"
  )
  
  ## Hospitalization rate
  
  if (scale_hospitalized) {
    hosp <- "hosp = f_h * (1/(1 + H/N));\n\n"
  } else {
    hosp <- "hosp = f_h;\n\n"
  }
  
  ## Transition rates
  
  rates <- c(
    "foi * dw/dt",
    "deathrate",
    "mu_latend * f_ai"
  )
  if (add_hospitalized) {
    rates <- c(
      rates, 
      "mu_latend * (1 - f_ai) * (1 - hosp)", 
      "mu_latend * (1 - f_ai) * hosp"
    )
  } else {
    rates <- c(
      rates, 
      "mu_latend * (1 - f_ai)"
    )
  }
  rates <- c(
    rates, 
    "deathrate",
    "mu_ar",
    "deathrate",
    "mu_ir",
    "deathrate"
  )
  if (add_hospitalized) {
    rates <- c(
      rates, 
      "mu_hr",
      "deathrate"
    )
  }
  if (use_separate_recovery) {
    rates <- c(
      rates,
      "mu_ras",
      "deathrate",
      "mu_ris", 
      "deathrate"
    )
  } else {
    rates <- c(
      rates,
      "mu_rs",
      "deathrate"
    )
  }
  if (add_vaccination) {
    rates <- c(
      rates,
      "f_v1 * foi",
      "mu_vs1",
      "deathrate",
      "f_v2 * foi",
      "mu_vs2",
      "deathrate"
    )
  }
  
  rates <- paste0(paste0("rate[", seq(0, length(rates) - 1, 1), "] = ", rates, collapse = ";\n"), ";\n\n")
  
  ## Transitions
  
  if (add_hospitalized) { 
    trans_num <- c(2, 4, 2, 2, 2)
  } else {
    trans_num <- c(2, 3, 2, 2)
  }
  if (use_separate_recovery) {
    trans_num <- c(trans_num, 2, 2)
  } else {
    trans_num <- c(trans_num, 2)
  }
  if (add_vaccination) {
    trans_num <- c(trans_num, 3, 3)
  }
  trans_starts <- cumsum(c(0, trans_num[1:length(trans_num) - 1]))
  trans_ends <- cumsum(trans_num)
  
  
  trans <- paste0(
    paste0(
      "reulermultinom(",
      trans_num,
      ", ",
      trans_comp,
      ", &rate[",
      trans_starts,
      "], dt, &trans[",
      trans_starts,
      "]);",
      collapse = "\n"
    ),
    "\n\n"
  )
  
  ## Compartment updates
  
  ### Outflows
  
  outflow_trans <- c()
  for (i in 1:length(trans_starts)) {
    outflow_idxs <- seq(trans_starts[i], trans_ends[i] - 1, 1)
    outflow_trans[[i]] <- paste0(
      "trans[",
      outflow_idxs,
      "]",
      collapse = "+"
    )
  }
  
  out_updates <- paste0(
    paste0(
      trans_comp,
      " -= (",
      outflow_trans,
      ");",
      collapse = "\n"
    ),
    "\n\n"
  )
  
  ### in_updatess
  
  rec_names <- trans_comp[which(startsWith(trans_comp, "R"))]
  rec_idxs <- which(startsWith(trans_comp, "R"))
  exp_idxs <- which(startsWith(trans_comp, "E"))
  asymp_idxs <- which(startsWith(trans_comp, "A"))
  symp_idxs <- which(startsWith(trans_comp, "I"))
  hosp_idxs <- which(startsWith(trans_comp, "H"))
  vacc_idxs <- which(startsWith(trans_comp, "V"))

  if (add_vaccination) {
    in_updates_s_idxs <- c(
      trans_starts[rec_idxs],
      trans_starts[vacc_idxs[[1]]] + 1
    )
  } else {
    in_updates_s_idxs <- c(
      trans_starts[rec_idxs]
    )
  }
  in_updates_s <- paste0(
    "S += births+",
    paste0(
      "trans[",
      in_updates_s_idxs,
      "]",
      collapse = "+"
    ),
    ";\n"
  )
  in_updates_e <- "E += trans[0];\n"
  if (add_vaccination) {
    in_updates_a_idxs <- c(
      trans_starts[exp_idxs],
      trans_starts[vacc_idxs]
    )
  } else {
    in_updates_a_idxs <- c(
      trans_starts[exp_idxs]
    )
  }
  in_updates_a <- paste0(
    "A += ",
    paste0(
      "trans[",
      in_updates_a_idxs,
      "]",
      collapse = "+"
    ),
    ";\n"
  )
  in_updates_i <- paste0(
    "I += trans[",
    trans_starts[exp_idxs] + 1,
    "];\n"
  )
  
  if (add_hospitalized) {
    in_updates_h <- paste0(
      "H += trans[",
      trans_starts[exp_idxs] + 2,
      "];\n"
    )
  }
  
  
  if (use_separate_recovery) {
    if (add_hospitalized) {
      in_updates_r_idxs <- list(trans_starts[asymp_idxs], c(trans_starts[symp_idxs], trans_starts[hosp_idxs]))
    } else {
      in_updates_r_idxs <- list(trans_starts[asymp_idxs], trans_starts[symp_idxs])
    }
    in_updates_r_parts <- c()
    for (i in 1:length(rec_names)) {
      in_updates_r_parts[[i]] <- paste0(
        rec_names[i],
        " += ",
        paste0(
          "trans[",
          unlist(in_updates_r_idxs[i]),
          "]",
          collapse = "+"
        )
      )
    }
    in_updates_r <- paste0(
      paste0(in_updates_r_parts, collapse = ";\n"),
      ";\n"
    )
  } else {
    if (add_hospitalized) {
      in_updates_r_idxs <- c(trans_starts[asymp_idxs], trans_starts[symp_idxs], trans_starts[hosp_idxs])
    } else {
      in_updates_r_idxs <- c(trans_starts[asymp_idxs], trans_starts[symp_idxs])
    }
    in_updates_r <- paste0(
      rec_names,
      " += ",
      paste0(
        "trans[",
        in_updates_r_idxs,
        "]",
        collapse = "+"
      ),
      ";\n"
    )
  }
  
  
  in_updates_v <- paste0("V1 += trans[", trans_starts[vacc_idxs[2]] + 1, "];\n")
  
  in_updates <- paste0(
    in_updates_s,
    in_updates_e,
    in_updates_a,
    in_updates_i,
    if (add_hospitalized) in_updates_h,
    in_updates_r,
    if (add_vaccination) in_updates_v,
    "\n"
  )
  
  if (add_hospitalized) {
    meas_state <- paste0(
      "C += trans[",
      trans_starts[exp_idxs] + 2,
      "];\n"
    )
  } else {
    meas_state <- paste0(
      "C += trans[",
      trans_starts[exp_idxs] + 1,
      "];\n"
    )
  }
  
  rprocess_str <- paste0(
    variable_definitions,
    if (add_idps) idps,
    if (add_vaccination) vacc,
    births,
    noise,
    foi,
    if (add_hospitalized) hosp,
    rates,
    trans,
    out_updates,
    in_updates,
    meas_state,
    "\n",
    "W += (dw - dt)/sigmaSE;\n\n",
    paste0("N = ", paste(trans_comp, collapse = "+"), ";\n")
  )
  
  rprocess <- Csnippet(rprocess_str)
  
  # Create pomp object
  
  m <- (
    df_model
    %>% pomp(
      times = "date",
      t0 = df_covar$date[1],
      rprocess = euler(
        step.fun = rprocess,
        delta.t = 1/52/7
      ),
      rinit = rinit,
      rmeasure = rmeas,
      dmeasure = dmeas,
      obsnames = obs_names,
      covar = covar,
      covarnames = covar_names,
      accumvars = accum_vars,
      statenames = states_names,
      paramnames = param_names,
      partrans = par_trans,
      params = fixed_params
    )
  )

  m
  # rprocess_str
}


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
    out_filename,
    reset_data = FALSE
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
  
  df_results <- (
    data.frame(
      t(sapply(m_global, coef)),
      log_lik=loglik_global[,1],
      log_lik_se=loglik_global[,2]
    )
  )
  
  end <- Sys.time()
  
  
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
  
  cat("Run time of global iterated particle filter: ", as.integer(round(difftime(end, start, units='mins'))), " minutes")
  
  df_results
  
}