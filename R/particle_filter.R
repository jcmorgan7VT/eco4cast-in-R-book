particle_filter <- function(t, forecast, obs, sim_dates, wt, fit_params){

  curr_obs <- obs |>
    filter(datetime == sim_dates[t],
           variable %in% c("lai", "wood", "som", "nee"))

  if(nrow(curr_obs) > 0){

    forecast_df <- output_to_df(forecast[t, , ], sim_dates[t], sim_name = "particle_filter")


    combined_output_obs <- combine_model_obs(forecast_df,
                                             obs = curr_obs,
                                             variables = c("lai", "wood", "som", "nee"),
                                             sds = c(0.1, 1, 20, 0.005))

    likelihood <- rep(0, ens_members)
    for(ens in 1:ens_members){

      curr_ens <- combined_output_obs |>
        filter(ensemble == ens)

      likelihood[ens] <- exp(sum(dnorm(x =  curr_ens$observation,
                                       mean = curr_ens$prediction,
                                       sd = curr_ens$sds, log = TRUE)))
    }

    wt[t, ] <- likelihood * wt[t-1, ]

    wt_norm <-  wt[t, ]/sum(wt[t, ])
    Neff = 1/sum(wt_norm^2)

    if(Neff < ens_members/2){
      ## resample ensemble members in proportion to their weight
      resample_index <- sample(1:ens_members, ens_members, replace = TRUE, prob = wt_norm )
      forecast[t, , ] <- as.matrix(forecast[t, resample_index, 1:12])  ## update state
      fit_params[t, , ] <- as.matrix(fit_params[t, resample_index, ])
      wt[t, ] <- rep(1, ens_members)
    }
  }
  return(list(forecast = forecast, fit_params = fit_params, wt = wt, sim_dates = sim_dates))
}
