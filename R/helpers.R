output_to_df <- function(output, sim_dates, sim_name){

  output_df <- NULL
  if(length(dim(output)) == 2){
    output_new <- array(0, dim=c(1,dim(output)))
    output_new[1, , ] <- output
    output <- output_new
  }
  ens_members <- dim(output)[2]

  for(ens in 1:ens_members){
    df <- as.data.frame(matrix(output[,ens, ], ncol = 12))
    names(df) <- c("leaves","wood","som", "lai", "gpp", "nee", "ra", "npp_l", "npp_w", "rh", "litterfall","mortality")
    df <- bind_cols(tibble(datetime = sim_dates, ensemble = ens), df)
    output_df <- bind_rows(output_df, df)
  }
  output_df <- output_df |>
    pivot_longer(-c("datetime","ensemble"), names_to = "variable", values_to = "prediction") |>
    mutate(model_id = sim_name)

  return(output_df)
}


parameters_to_df <- function(fit_params, sim_dates, sim_name, param_names){

  output_df <- NULL
  ens_members <- dim(fit_params)[2]

  for(ens in 1:ens_members){
    df <- as.data.frame(fit_params[,ens, ])
    names(df) <- param_names
    df <- bind_cols(tibble(datetime = sim_dates, ensemble = ens), df)
    output_df <- bind_rows(output_df, df)
  }
  output_df <- output_df |>
    pivot_longer(-c("datetime","ensemble"), names_to = "variable", values_to = "prediction") |>
    mutate(model_id = sim_name)
}



combine_model_obs <- function(out, obs, variables, sds){

  combined <- out |>
    left_join(obs, by = c("datetime", "variable")) |>
    filter(variable %in% variables) |>
    na.omit() |>
    left_join(tibble(variable = variables,
                     sds = sds), by = "variable")

  return(combined)
}

assign_met_ensembles <- function(inputs, ens_members, var_order = c("temp", "PAR")){

  met_ensembles <- unique(inputs$parameter)

  dates <- unique(inputs$datetime)

  inputs <- inputs |>
    pivot_wider(names_from = variable, values_from = prediction) |>
    select(all_of(c("datetime", "parameter", var_order))) |>
    arrange(datetime, parameter)

  var_names <- names(inputs |> select(-datetime, -parameter))

  inputs_ensemble <- array(0, dim = c(length(dates),ens_members, length(var_names)))
  for(ens in 1:ens_members){
    ensemble_index <- sample(met_ensembles, 1)
    curr_inputs <- inputs |>
      filter(parameter == ensemble_index) |>
      select(-datetime, -parameter)

    inputs_ensemble[, ens, ] <- as.matrix(curr_inputs)
  }
  dimnames(inputs_ensemble)[[3]] <- var_names
  return(inputs_ensemble)
}


