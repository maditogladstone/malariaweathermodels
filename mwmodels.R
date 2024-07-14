# This script provides the code used in the study comparing different approaches of modelling the effects of temperature and rainfall in high and low transmission settings

# packages and themes ####
## libraries
library(pacman)
p_load(
  readxl,
  tidyverse,
  deSolve,
  scales
)

## plots theme
theme_manuscript <- function(){
  theme_minimal() + 
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 14, color = "black"), 
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title = element_text(size = 14, colour = 'black'),
          legend.text = element_text(size = 14, colour = 'black'),
          legend.key.height = unit(1, "cm"),
          legend.position = "bottom"
    )
}

# data ####
## initial conditions 
model_states <- function(){
  
  initial_states <- read_excel("vector_states_parameters.xlsx", sheet = "States")
  
  temp_states <- initial_states %>% 
    select(!description) %>% 
    pivot_wider(names_from = "state", values_from = "value") %>% 
    {
      list(
        aquatic = as.list(select(., ends_with("_a"))),
        adults = as.list(select(., ends_with("_m"))),
        hosts = as.list(select(., ends_with("_h"))),
        counts = as.list(select(., starts_with("C")))
      )
    }
  
  return(temp_states)
}

## base model parameters
model_parameters <- function(){
  
  base_parameters <- read_excel("vector_states_parameters.xlsx", sheet = "Parameters")
  
  temp_parameters <- base_parameters %>% 
    select(!description) %>% 
    pivot_wider(names_from = "parameter", values_from = "value") %>% 
    {
      as.list(.)
    }
  
  return(temp_parameters)
}

## temperature and rainfall  
## datasets were obtained from: https://cckpapi.worldbank.org/cckp/v1/era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_mean/ZAF?_format=json
cckp_weather <- function(){
  
  temps_values <- read_excel("era5-x0.25_timeseries_tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_mean.xlsx")
  
  rains_values <- read_excel("era5-x0.25_timeseries_pr_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_mean.xlsx")
  
  temp_weather <- function(values){
    values %>%
      select(!c("code", "name")) %>% 
      pivot_longer(cols = everything(), names_to = "time", values_to = "value") %>% 
      mutate(time = 31*seq(1, nrow(.), 1))
  } 
  
  temps_t <- approxfun(temp_weather(temps_values)$time, temp_weather(temps_values)$value, rule = 2:1)
  rains_t <- approxfun(temp_weather(rains_values)$time, temp_weather(rains_values)$value, rule = 2:1)
  
  return(
    list(
      temps = temps_t, rains = rains_t
    )
  )
}

## parameter functions 
### eggs laid
n_e_func <- function(times, temps){
  
  temp_n_e <- -0.61411*(temps(times) + 2)^3 + 38.93*(temps(times) + 2)^2 - 801.27*(temps(times) + 2) + 5391.4
  
  return(
    pmax(120, temp_n_e)
  )
}

### egg deposition rate
theta_func <- function(times, temps){
  
  temp_theta <- 0.00054*temps(times)^3 - 0.038*temps(times)^2 + 0.88*temps(times) + 1
  
  return(
    pmax(0, temp_theta)
  )
}

### biting rates
weather_a_func <- function(times, temps){
  
  temp_a <- 0.000203*(temps(times)^2 - 11.7*temps(times))*sqrt(42.3 - temps(times))
  
  return(
    pmax(0, temp_a)
  )
}

weather_fit_a_func <- function(times, temps){
  
  return(
    1/(107.204 - 13.3523*temps(times) + 0.677509*temps(times)^2 - 0.0159732*temps(times)^3 + 0.000144876*temps(times)^4)
  )
}

### adult mosquito mortality rates
thermo_mu_m_func <- function(times, temps){
  
  temp_mu_m <- -0.000091*temps(times)^3 + 0.059*temps(times)^2 + 1.3*temps(times) + 9.9
  
  return(
    pmax(1/21, temp_mu_m)
  )
}

hydro_mu_m_func <- function(times, temps){
  
  return(
    0.07/exp(-4.40 + 1.31*temps(times) - 0.03*temps(times)^2)
  )
}

weather_fit_mu_m_func <- function(times, temps){
  
  return(
    (3.04 + 29.564*exp(-(temps(times) + 273.15 - 278)/2.7035))/30.4
  )
}

### egg clearance rate
thermo_mu_e_func <- function(times, temps){
  
  temp_mu_e <- 0.0033*(temps(times) + 2)^3 - 0.23*(temps(times) + 2)^2 + 5.3*(temps(times) + 2) - 40
  
  return(
    pmax(0.5, temp_mu_e)    
  )
}

### larval mortality rates
thermo_mu_l_func <- function(times, temps){
  
  temp_mu_l <- 0.00081*(temps(times) + 2)^3 - 0.056*(temps(times) + 2)^2 + 1.3*(temps(times) + 2) - 8.6
  
  return(
    pmax(0.4, temp_mu_l)    
  )
}

hydro_mu_l_func <- function(times, temps){
  
  return(
    0.0025*temps(times)^2 - 0.094*temps(times) + 1.0257
  )
}

weather_fit_mu_l_func <- function(times, temps){
  
  return(
    1/(-4.4 + 1.31*temps(times) - 0.03*temps(times)^2)
  )
}

### pupae clearance rate
thermo_mu_p_func <- function(times, temps){
  
  temp_mu_p <- 0.0034*(temps(times) + 2)^3 - 0.22*(temps(times) + 2)^2 - 4.9*(temps(times) + 2) - 34
  
  return(
    pmax(0.3, temp_mu_p)    
  )
}

### parasite development rate in mosquitoes
gamma_m_func <- function(times, temps){
  
  temp_gamma_m <- (temps(times) - 16)/111
  
  return(
    case_when(temps(times) < 16 ~ 0, .default = temp_gamma_m)
  )
}

### egg environmental carrying capacity
K_e_func <- function(times, variables, rains){
  
  temp_K_e <- (variables$parameters$constants[["P_A"]]*0.225/0.01)*rains(times)
  
  return(
    case_when(temp_K_e == 0 ~ 1e8, .default = temp_K_e)
  )
}

### egg development rate
thermo_kappa_e_func <- function(times, temps){
  
  temp_kappa_e <- 0.012*(temps(times) + 2)^3 - 0.81*(temps(times) + 2)^2 + 18*(temps(times) + 2) - 135.93
  
  return(
    pmax(1/2, temp_kappa_e)
  )
}

### larvale development rate
thermo_kappa_l_func <- function(times, temps){
  
  temp_kappa_l <- -0.002*(temps(times) + 2)^3 + 0.14*(temps(times) + 2)^2 - 3*(temps(times) + 2) + 22
  
  return(
    pmax(1/3, temp_kappa_l)
  )
}

### pupae development rates
thermo_kappa_p_func <- function(times, temps){
  
  temp_kappa_p <- 0.0034*(temps(times) + 2)^3 - 0.22*(temps(times) + 2)^2 - 4.9*(temps(times) + 2) - 34
  
  return(
    pmax(1/4, temp_kappa_p)    
  )
}

hydro_kappa_p_func <- function(times, variables, temps, rains){
  
  # larval development period
  n_e <- variables$parameters$constants[["n_e"]]
  t_e <- variables$parameters$constants[["t_e"]]
  t_l <- 1/(0.0557*(temps(times) + 2) - 0.06737)
  t_p <- variables$parameters$constants[["t_p"]]
  
  # probability of survival
  R_l <- variables$parameters$constants[["R_l"]]
  max_p_e <- variables$parameters$constants[["max_p_e"]]
  max_p_l <- variables$parameters$constants[["max_p_l"]]
  max_p_p <- variables$parameters$constants[["max_p_p"]]
  p_e <- min(max(0, (4*max_p_e/R_l^2)*rains(times)*(R_l - rains(times))), max_p_e)
  p_l <- exp(-1/t_l)*min(max(0, (4*max_p_l/R_l^2)*rains(times)*(R_l - rains(times))), max_p_l)
  p_p <- min(max(0, (4*max_p_p/R_l^2)*rains(times)*(R_l - rains(times))), max_p_p)
  
  return(
    (n_e*p_e*p_l*p_p)/(t_e + t_l + t_p)
  )
}

parameters_func <- list(
  n_e_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = rep(variables$parameters$constants[["n_e"]], length(times)),
        Approach_B = n_e_func(times, temps),
        Approach_A = rep(variables$parameters$constants[["n_e"]], length(times))
      )
    )
  },
  theta_wdm = function(times, temps){
    
    return(
      list(
        Approach_C = rep(variables$parameters$constants[["theta"]], length(times)),
        Approach_B = theta_func(times, temps),
        Approach_A = rep(variables$parameters$constants[["theta"]], length(times))/3
      )
    )
  },
  a_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = 10*weather_a_func(times, temps),
        Approach_B = 500*weather_a_func(times, temps),
        Approach_A = 10*weather_fit_a_func(times, temps)
      )
    )
  },
  mu_m_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = hydro_mu_m_func(times, temps),
        Approach_B = thermo_mu_m_func(times, temps),
        Approach_A = weather_fit_mu_m_func(times, temps)
      )
    )
  },
  mu_e_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = rep(variables$parameters$constants[["mu_e"]], length(times)),
        Approach_B = thermo_mu_e_func(times, temps),
        Approach_A = rep(variables$parameters$constants[["mu_e"]], length(times))
      )
    )
  },
  mu_l_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = hydro_mu_l_func(times, temps),
        Approach_B = thermo_mu_l_func(times, temps),
        Approach_A = weather_fit_mu_l_func(times, temps)
      )
    )
  },
  mu_p_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = rep(variables$parameters$constants[["mu_p"]], length(times)),
        Approach_B = thermo_mu_p_func(times, temps),
        Approach_A = rep(variables$parameters$constants[["mu_p"]], length(times))      )
    )
  },
  gamma_m_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = gamma_m_func(times, temps)/100,
        Approach_B = 100*rep(variables$parameters$constants[["gamma_m"]], length(times)),
        Approach_A = rep(variables$parameters$constants[["gamma_m"]], length(times))/2
      )
    )
  },
  K_e_wdm = function(times, variables, rains){
    
    return(
      list(
        Approach_C = K_e_func(times, variables, rains),
        Approach_B = rep(variables$parameters$constants[["K_e"]], length(times)),
        Approach_A = rep(variables$parameters$constants[["K_e"]], length(times))
      )
    )
  },
  kappa_e_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = 10*rep(variables$parameters$constants[["kappa_e"]], length(times)),
        Approach_B = thermo_kappa_e_func(times, temps),
        Approach_A = rep(variables$parameters$constants[["kappa_e"]], length(times))
      )
    )
  },
  kappa_l_wdm = function(times, variables, temps){
    
    return(
      list(
        Approach_C = rep(variables$parameters$constants[["kappa_l"]], length(times)),
        Approach_B = thermo_kappa_l_func(times, temps),
        Approach_A = rep(variables$parameters$constants[["kappa_l"]], length(times))
      )
    )
  },
  kappa_p_wdm = function(times, variables, temps, rains){
    
    return(
      list(
        Approach_C = hydro_kappa_p_func(times, variables, temps, rains),
        Approach_B = thermo_kappa_p_func(times, temps),
        Approach_A = rep(variables$parameters$constants[["kappa_p"]], length(times))/3
      )
    )
  }
)

## vector control intervention
vc_efficacy = list(
  "1" = 0, 
  "2" = 0.25,
  "3" = 0.50,
  "4" = 0.75,
  "5" = 0.90
)

scenarios = list(
  "1" = "0% efficacy",
  "2" = "25% efficacy",
  "3" = "50% efficacy",
  "4" = "75% efficacy",
  "5" = "90% efficacy"
)

# model ####
variables <- list(
  statesNames = list(
    # total populations
    tot_pop = list(
      "tot_a",
      "tot_m",
      "tot_h"
    ),
    # vectors 
    vectors = list(
      ## aquatic 
      aquatic = list(
        "E_a",
        "L_a",
        "P_a"
      ),
      ## adults 
      adults = list(
        "S_m",
        "E_m",
        "I_m"
      )
    ),
    # hosts 
    hosts = list(
      "S_h",
      "E_h",
      "A_h",
      "Iu_h",
      "Is_h",
      "Tu_h",
      "Ts_h",
      "R_h"
    ),
    # counts 
    counts = list(
      "CInc",
      "CSev",
      "CTrt_u",
      "CTrt_s",
      "CDth",
      "CMov", 
      "CPrv",
      "CVec"
    )
  ),
  states = model_states(),
  parameters = list(
    constants = model_parameters(),
    functions = parameters_func
  ),
  effects = list(
    weather = list(
      "cckp" = cckp_weather()
    )
  ),
  interventions  = list(
    vector_control = vc_efficacy
  )
)

# Simulator ####
simulate <- function(){
  
  # model
  mod <- function(times, state, parameters, vars, effects_func, interventions_func) {
    with(as.list(c(state, parameters, effects_func, interventions_func)), {
      
      # total number of eggs laid 
      n_e <- effects_func[["n_e"]][times + 1]
      
      # egg deposition rate
      theta <- effects_func[["theta"]][times + 1]
      
      # carrying capacity
      K_e <- effects_func[["K_e"]][times + 1]
      
      # death rates
      mu_e <- effects_func[["mu_e"]][times + 1]
      mu_l <- effects_func[["mu_l"]][times + 1]
      mu_p <- effects_func[["mu_p"]][times + 1]
      mu_m <- effects_func[["mu_m"]][times + 1]
      
      # development rate of mosquitoes
      kappa_e <- effects_func[["kappa_e"]][times + 1]
      kappa_l <- effects_func[["kappa_l"]][times + 1]
      kappa_p <- effects_func[["kappa_p"]][times + 1]
      
      # mosquito progression rate
      gamma_m <- effects_func[["gamma_m"]][times + 1]
      
      # biting rate
      a <- effects_func[["a"]][times + 1] 
      
      # population control 
      eta <- 0 # movement
      mu_b <- mu_h # birth
      mu_s <- 0 # disease induced death 
      
      # total populations 
      tot_a = E_a + L_a + P_a
      tot_m = S_m + E_m + I_m
      tot_h = S_h + E_h + A_h + Iu_h + Is_h + Tu_h + Ts_h + R_h
      
      # vector control
      zeta = interventions_func[["zeta"]][times + 1]
      sigma = interventions_func[["sigma"]][times + 1]
      eff_vc_cov = 1/(1 + interventions_func[["vc_eff"]][times + 1]*CVec)
      
      # infectious humans 
      Infectious <- zeta_a*A_h + Iu_h + Is_h + zeta_u*Tu_h + zeta_s*Ts_h
      
      # force of infection 
      lambda_m <- eff_vc_cov*a*b*(Infectious/tot_h)
      lambda_h <- eff_vc_cov*a*c*(tot_m/tot_h)*(I_m/tot_m)
      
      ## aquatic mosquitoes
      dE_a = n_e*theta*(1 - tot_m/K_e)*tot_m - kappa_e*E_a - mu_e*E_a
      dL_a = kappa_e*E_a - kappa_l*L_a - mu_l*L_a
      dP_a = kappa_l*L_a - kappa_p*P_a - mu_p*P_a
      
      ## adult mosquitoes
      dS_m = kappa_p*P_a - lambda_m*S_m - mu_m*S_m
      dE_m = lambda_m*S_m - gamma_m*E_m - mu_m*E_m
      dI_m = gamma_m*E_m - mu_m*I_m
      
      ## host population
      dS_h = mu_b*tot_h - lambda_h*S_h - mu_h*S_h + rho*R_h
      dE_h = lambda_h*S_h - gamma_h*E_h - mu_h*E_h
      dA_h = pa*gamma_h*E_h + omega*Iu_h - delta_r*A_h - mu_h*A_h
      dIu_h = eta + (1 - pa)*gamma_h*E_h - delta_r*Iu_h - omega*Iu_h - nu*Iu_h - tau_u*Iu_h - mu_h*Iu_h + alpha*Is_h
      dIs_h = nu*Iu_h - tau_s*Is_h - alpha*Is_h - mu_h*Is_h - mu_s*Is_h
      dTu_h = tau_u*Iu_h - delta_u*Tu_h - mu_h*Tu_h
      dTs_h = tau_s*Is_h - delta_s*Ts_h - mu_h*Ts_h
      dR_h = delta_u*Tu_h + delta_s*Ts_h + delta_r*A_h + delta_r*Iu_h - rho*R_h - mu_h*R_h
      
      ## counts 
      dCInc = lambda_h*S_h
      dCSev = nu*Iu_h
      dCTrt_u = tau_u*Iu_h
      dCTrt_s = tau_s*Is_h
      dCDth = mu_s*Is_h
      dCMov = eta
      dCPrv = A_h + Iu_h + Is_h + Tu_h + Ts_h
      dCVec = zeta - sigma*CVec
      
      # output 
      output <- c(
        dE_a, dL_a, dP_a, 
        dS_m, dE_m, dI_m, 
        dS_h, dE_h, dA_h, dIu_h, dIs_h, dTu_h, dTs_h, dR_h, 
        dCInc, dCSev, dCTrt_u, dCTrt_s, dCDth, dCMov, dCPrv,
        dCVec
      )
      
      return(list(output)) 
    })
  }
  
  # time
  no_yrs <- 1
  times <- seq(0, 365*no_yrs, 1)
  
  # parameters 
  parms <- variables$parameters$constants
  
  # states
  initstates <- variables$states %>%
    list_c() %>%
    as_vector()
  
  # effects 
  model_effects <- function(times, variables, mod_type){
    
    temps <- variables$effects$weather[["cckp"]][["temps"]]
    
    rains <- variables$effects$weather[["cckp"]][["rains"]]
    
    return(
      list(
        n_e = variables$parameters$functions$n_e_wdm(times, variables, temps)[[mod_type]],
        theta = variables$parameters$functions$theta_wdm(times, temps)[[mod_type]],
        a = variables$parameters$functions$a_wdm(times, variables, temps)[[mod_type]],
        mu_m = variables$parameters$functions$mu_m_wdm(times, variables, temps)[[mod_type]],
        mu_e = variables$parameters$functions$mu_e_wdm(times, variables, temps)[[mod_type]],
        mu_l = variables$parameters$functions$mu_l_wdm(times, variables, temps)[[mod_type]],
        mu_p = variables$parameters$functions$mu_p_wdm(times, variables, temps)[[mod_type]],
        gamma_m = variables$parameters$functions$gamma_m_wdm(times, variables, temps)[[mod_type]],
        K_e = variables$parameters$functions$K_e_wdm(times, variables, rains)[[mod_type]],
        kappa_e = variables$parameters$functions$kappa_e_wdm(times, variables, temps)[[mod_type]],
        kappa_l = variables$parameters$functions$kappa_l_wdm(times, variables, temps)[[mod_type]],
        kappa_p = variables$parameters$functions$kappa_p_wdm(times, variables, temps, rains)[[mod_type]]
      )
    )
  }
  
  # interventions
  vec_control <- "itn"
  
  model_interventions <- function(times, variables, no_scenario, vector_control){
    
    intervention_scenario <- switch(
      vec_control,
      "irs" = list(
        zeta = rep(1/365, length(times)), # annual spraying
        sigma = rep(2/365, length(times)), # 6 month operational effectiveness
        vc_eff = rep(variables$interventions$vector_control[[no_scenario]], length(times)) # efficacy
      ),
      "itn" = list(
        zeta = rep(1/365, length(times)), # annual deployment
        sigma = rep(1/(365*3), length(times)), # 3 year operational effectiveness
        vc_eff = rep(variables$interventions$vector_control[[no_scenario]], length(times)) # efficacy
      )
    )
    
    return(
      intervention_scenario
    )
  }
  
  # simulate
  int_type <- 1:5
  
  sim_model <- function(int_type, mod_type){
    
    model_output <- ode(
      times = times,
      y = initstates,
      parms = parms,
      vars = variables,
      effects = model_effects(times, variables, mod_type),
      interventions = model_interventions(times, variables, int_type),
      func = mod
    )
    
    return(model_output)
  }
  
  result <- list(
    Approach_A = sapply(int_type, sim_model, "Approach_A", simplify = "array"),
    Approach_B = sapply(int_type, sim_model, "Approach_B", simplify = "array"),
    Approach_C = sapply(int_type, sim_model, "Approach_C", simplify = "array")
  )
  
  return(result)
}

# Postprocessor ####
postproc <- function(result){
  
  convert_to_longer <- function(result, mod_type){
    
    temp_result <- result[[mod_type]] %>% 
      as.data.frame() %>%
      mutate(
        # total populations
        tot_a.1 = select(., ends_with("_a.1")) %>% apply(., 1, sum),
        tot_a.2 = select(., ends_with("_a.2")) %>% apply(., 1, sum),
        tot_a.3 = select(., ends_with("_a.3")) %>% apply(., 1, sum),
        tot_a.4 = select(., ends_with("_a.4")) %>% apply(., 1, sum),
        tot_a.5 = select(., ends_with("_a.5")) %>% apply(., 1, sum),
        
        tot_m.1 = select(., ends_with("_m.1")) %>% apply(., 1, sum),
        tot_m.2 = select(., ends_with("_m.2")) %>% apply(., 1, sum),
        tot_m.3 = select(., ends_with("_m.3")) %>% apply(., 1, sum),
        tot_m.4 = select(., ends_with("_m.4")) %>% apply(., 1, sum),
        tot_m.5 = select(., ends_with("_m.5")) %>% apply(., 1, sum),
        
        tot_h.1 = select(., ends_with("_h.1")) %>% apply(., 1, sum),
        tot_h.2 = select(., ends_with("_h.2")) %>% apply(., 1, sum),
        tot_h.3 = select(., ends_with("_h.3")) %>% apply(., 1, sum),
        tot_h.4 = select(., ends_with("_h.4")) %>% apply(., 1, sum),
        tot_h.5 = select(., ends_with("_h.5")) %>% apply(., 1, sum),
        
        # counts
        Inc.1 = c(0, diff(CInc.1)), # incidence
        Inc.2 = c(0, diff(CInc.2)),
        Inc.3 = c(0, diff(CInc.3)),
        Inc.4 = c(0, diff(CInc.4)),
        Inc.5 = c(0, diff(CInc.5)),
        
        Sev.1 = c(0, diff(CSev.1)), # severe cases
        Sev.2 = c(0, diff(CSev.2)),
        Sev.3 = c(0, diff(CSev.3)),
        Sev.4 = c(0, diff(CSev.4)),
        Sev.5 = c(0, diff(CSev.5)),
        
        Trt_u.1 = c(0, diff(CTrt_u.1)), # treated uncomplicated cases
        Trt_u.2 = c(0, diff(CTrt_u.2)),
        Trt_u.3 = c(0, diff(CTrt_u.3)),
        Trt_u.4 = c(0, diff(CTrt_u.4)),
        Trt_u.5 = c(0, diff(CTrt_u.5)),
        
        Trt_s.1 = c(0, diff(CTrt_s.1)), # treated severe cases
        Trt_s.2 = c(0, diff(CTrt_s.2)),
        Trt_s.3 = c(0, diff(CTrt_s.3)),
        Trt_s.4 = c(0, diff(CTrt_s.4)),
        Trt_s.5 = c(0, diff(CTrt_s.5)),
        
        Dth.1 = c(0, diff(CDth.1)), # disease induced death
        Dth.2 = c(0, diff(CDth.2)),
        Dth.3 = c(0, diff(CDth.3)),
        Dth.4 = c(0, diff(CDth.4)),
        Dth.5 = c(0, diff(CDth.5)),
        
        Mov.1 = c(0, diff(CMov.1)), # incidence of infections acquired from outside the community
        Mov.2 = c(0, diff(CMov.2)),
        Mov.3 = c(0, diff(CMov.3)),
        Mov.4 = c(0, diff(CMov.4)),
        Mov.5 = c(0, diff(CMov.5)),
        
        Prv.1 = c(0, diff(CPrv.1)), # prevalence
        Prv.2 = c(0, diff(CPrv.2)),
        Prv.3 = c(0, diff(CPrv.3)),
        Prv.4 = c(0, diff(CPrv.4)),
        Prv.5 = c(0, diff(CPrv.5)),
        
        time = time.1 # time
      ) %>%
      select(!starts_with("time.")) %>% 
      pivot_longer(!time, names_to = "variable", values_to = "value") %>% 
      mutate(
        model = mod_type,
        scenario = case_when(
          substring(variable, nchar(variable), nchar(variable)) == "1" ~ as_factor(scenarios[["1"]]),
          substring(variable, nchar(variable), nchar(variable)) == "2" ~ as_factor(scenarios[["2"]]),
          substring(variable, nchar(variable), nchar(variable)) == "3" ~ as_factor(scenarios[["3"]]),
          substring(variable, nchar(variable), nchar(variable)) == "4" ~ as_factor(scenarios[["4"]]),
          substring(variable, nchar(variable), nchar(variable)) == "5" ~ as_factor(scenarios[["5"]])
        ),
        state_var = substring(variable, 1, nchar(variable) - 2)
      )
    
    return(temp_result)
  }
  
  # post processing
  thermal_fit_pproc <- convert_to_longer(result, "Approach_A")
  
  thermo_pproc <- convert_to_longer(result, "Approach_B")
  
  hydro_pproc <- convert_to_longer(result, "Approach_C")
  
  pproc <- bind_rows(
    thermal_fit_pproc, thermo_pproc, hydro_pproc
  )
  
  return(pproc)
}

# Plotter ####
plotter <- function(pproc){

  skip_yrs <- 0 # start plot after skip years
  
  # inputs
  temps_plotter <- function(pproc, plot_type){
    
    temps_plot <- ggplot() +
      geom_line(
        aes(
          x = unique(filter(pproc, time >= 365*skip_yrs)$time), # convert to years 
          y = variables$effects$weather[[plot_type]]$temps(unique(filter(pproc, time >= 365*skip_yrs)$time))
        ),
        colour = "red",
        linewidth = 0.5
      ) +
      geom_point(
        aes(
          x = seq(365*skip_yrs, max(pproc$time), by = 31), # convert to months 
          y = variables$effects$weather[[plot_type]]$temps(seq(365*skip_yrs, max(pproc$time), by = 31))
        ),
        colour = "red",
        size = 2
      ) +
      scale_x_continuous(
        breaks = seq(0, max(pproc$time), by = 365), 
        labels = 0:(max(pproc$time) %/% 365)
      ) +
      ylim(0, NA) +
      labs(
        title = "Mean Monthly Temperature in South Africa",
        x = "Time (years)",
        y = expression("Temperature ("*degree*"C)")
      ) +
      theme_manuscript()
    
    
    return(temps_plot)
  }
  
  rains_plotter <- function(pproc, plot_type){
    
    rains_plot <- ggplot() +
      geom_line(
        aes(
          x = unique(filter(pproc, time >= 365*skip_yrs)$time), # convert to years 
          y = variables$effects$weather[[plot_type]]$rains(unique(filter(pproc, time >= 365*skip_yrs)$time))
        ),
        colour = "blue",
        linewidth = 0.5
      ) +
      geom_point(
        aes(
          x = seq(365*skip_yrs, max(pproc$time), by = 31), # convert to months 
          y = variables$effects$weather[[plot_type]]$rains(seq(365*skip_yrs, max(pproc$time), by = 31))
        ),
        colour = "blue",
        size = 2
      ) +
      scale_x_continuous(
        breaks = seq(0, max(pproc$time), by = 365), 
        labels = 0:(max(pproc$time) %/% 365)
      ) +
      ylim(0, NA) +
      labs(
        title = "Mean Monthly Rainfall in South Africa",
        x = "Time (years)",
        y = "Rainfall (mm)"
      ) +
      theme_manuscript()
    
    
    return(rains_plot)
  }
  
  # parameters
  thermo_time <- seq(0, 42, 1) # 0 to 42 degrees Celcius
  
  thermal <- function(times){
    
    temp = times     
    
    return(
      temp
    )
  }
  
  a_plotter <- function(thermo_time){
    
    a_plot <- ggplot() +
      geom_line(
        aes(
          x = thermo_time, 
          y = 10*weather_a_func(thermo_time, thermal)
        ),
        colour = "red",
        linewidth = 1
      ) +
      xlim(0, NA) +
      ylim(0, NA) +
      labs(
        title = "Mosquito Biting Rate (Entomological Inoculation Rate)",
        x = expression("temperature ("*degree*"C)"),
        y = "Bites/day"
      ) +
      theme_manuscript()
    
    
    return(a_plot)
  }
  
  pdr_plotter <- function(thermo_time){
    
    gamma_m_plot <- ggplot() +
      geom_line(
        aes(
          x = thermo_time,  
          y = gamma_m_func(thermo_time, thermal)
        ),
        colour = "red",
        linewidth = 1
      ) +
      xlim(0, NA) +
      ylim(0, NA) +
      labs(
        title = "Parasite Development Rate",
        x = expression("Temperature ("*degree*"C)"),
        y = "Infectious mosquitoes/day"
      ) +
      theme_manuscript()
    
    return(gamma_m_plot)
  }
  
  lmr_plotter <- function(thermo_time){
    
    mu_l_plot <- ggplot() +
      geom_line(
        aes(
          x = thermo_time, 
          y = hydro_mu_l_func(thermo_time, thermal)
        ),
        colour = "red",
        linewidth = 1
      ) +
      xlim(0, NA) +
      ylim(0, NA) +
      labs(
        title = "Larval Mortality Rate",
        x = expression("Temperature ("*degree*"C)"),
        y = "Larvae/day"
      ) +
      theme_manuscript()
    
    
    return(mu_l_plot)
  }
  
  # projections  
  inc_plotter <- function(pproc, plot_type){
    
    inc_plot <- pproc %>% 
      filter(
        time >= 365*skip_yrs,
        state_var %in% "Inc"
      ) %>% 
      group_by(state_var, model) %>% 
      ggplot() +
      geom_line(
        aes(
          x = time, y = value, colour = as_factor(scenario)#, linetype = scenario
        ),
        linewidth = 1
      ) +
      facet_wrap(~as_factor(model), scales = "free", ncol = 1) +
      scale_y_continuous(
        limits = c(0, NA),
        labels = label_number(suffix = "", scale = 1e-3)
      ) +
      scale_x_continuous(
        breaks = seq(0, max(pproc$time), by = 365), 
        labels = 0:(max(pproc$time) %/% 365)
      ) +
      labs(
        title = "Incidence",
        x = "Time (years)",
        y = "Cases per 1 000",
        colour = "Scenario"
      ) +
      theme_manuscript()
    
    return(inc_plot)
  }
  
  uncomplicated_plotter <- function(pproc, plot_type){
    
    uncomplicated_plot <- pproc %>% 
      filter(
        time >= 365*skip_yrs,
        state_var %in% "Iu_h"
      ) %>% 
      group_by(state_var, model) %>% 
      ggplot() +
      geom_line(
        aes(x = time, y = value, colour = as_factor(model)),
        linewidth = 1
      ) +
      facet_wrap(~as_factor(scenario), scales = "free") +
      scale_y_continuous(
        limits = c(0, NA),
        labels = label_number(suffix = "", scale = 1e-3)
      ) +
      scale_x_continuous(
        breaks = seq(0, max(pproc$time), by = 365), 
        labels = 0:(max(pproc$time) %/% 365)
      ) +
      labs(
        title = "Uncomplicated Infections",
        x = "Time (years)",
        y = "Cases per 1 000",
        colour = "Model"
      ) +
      theme_manuscript()
    
    return(uncomplicated_plot)
  }
  
  prv_plotter <- function(pproc, plot_type){
    
    prv_plot <- pproc %>% 
      filter(
        time >= 365*skip_yrs,
        state_var %in% "Prv" 
      ) %>% 
      group_by(state_var, model) %>% 
      ggplot() +
      geom_line(
        aes(
          x = time, y = value, colour = as_factor(scenario)#, linetype = scenario
        ),
        linewidth = 1
      ) +
      facet_wrap(~as_factor(model), scales = "free", ncol = 1) +
      scale_y_continuous(
        limits = c(0, NA),
        labels = label_number(suffix = "", scale = 1e-3)
      ) +
      scale_x_continuous(
        breaks = seq(0, max(pproc$time), by = 365), 
        labels = 0:(max(pproc$time) %/% 365)
      ) +
      labs(
        title = "Prevalence",
        x = "Time (years)",
        y = "Cases per 1 000",
        colour = "Scenario"
      ) +
      theme_manuscript()
    
    return(prv_plot)
  }
  
  # cases averted
  averted_func <- function(approach, pproc){
    
    temp_averted <- cbind(
      "25% efficacy" = sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "0% efficacy" & state_var == "Iu_h" & model == approach), value)) - sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "25% efficacy" & state_var == "Iu_h" & model == approach), value)),
      "50% efficacy" = sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "0% efficacy" & state_var == "Iu_h" & model == approach), value)) - sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "50% efficacy" & state_var == "Iu_h" & model == approach), value)), 
      "75% efficacy" = sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "0% efficacy" & state_var == "Iu_h" & model == approach), value)) - sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "75% efficacy" & state_var == "Iu_h" & model == approach), value)), 
      "90% efficacy" = sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "0% efficacy" & state_var == "Iu_h" & model == approach), value)) - sum(select(filter(pproc, time >= 365*skip_yrs, scenario == "90% efficacy" & state_var == "Iu_h" & model == approach), value))
    )
    
    return(temp_averted)
  }
  
  mod_approaches <- c("Approach_A", "Approach_B", "Approach_C")
  
  cases_averted <- sapply(mod_approaches, averted_func, pproc = pproc, simplify = "matrix")
  
  dimnames(cases_averted)[[1]] <- scenarios[2:5]
  dimnames(cases_averted)[[2]] <- mod_approaches
  
  cases_averted_plotter <- function(pproc, plot_type){
    
    cases_averted_plot <- cases_averted %>% 
      as.data.frame() %>%
      mutate(scenario = rownames(.)) %>%
      pivot_longer(!scenario, names_to = "Model", values_to = "Cases") %>% 
      ggplot() +
      geom_col(aes(y = scenario, x = Cases, fill = Model), position = position_dodge()) +
      scale_x_continuous(
        limits = c(0, NA),
        labels = scales::label_number(suffix = "", scale = 1e-5)
      ) +
      scale_colour_manual(
        values = c("firebrick", "green4", "blue3"),
        aesthetics = c("colour", "fill")
      ) +
      labs(
        title = "Cases Averted",
        y = "Scenario",
        x = "Cases per 100 000",
        fill = "Model"
      ) +
      theme_manuscript()
    
    return(cases_averted_plot)
  }
  
  # plots
  plots <- list(  
    # number of bites by a mosquito per day
    "bites" = a_plotter(thermo_time),
    
    # daily larval mortality rate
    "lmr" = lmr_plotter(thermo_time),
    
    # daily parasite development rate
    "pdr" = pdr_plotter(thermo_time),
    
    # temps 
    "temps" = temps_plotter(pproc, "cckp"),
    
    # rains 
    "rains" = rains_plotter(pproc, "cckp"),

    # incidence
    "incidence" = inc_plotter(pproc, "inc_plot"),

    # uncomplicated
    "uncomplicated" = uncomplicated_plotter(pproc, "uncomplicated_plot"),

    # prevalence
    "prevalence" = prv_plotter(pproc, "prv_plot"),

    # cases averted
    "casesaverted" = cases_averted_plotter(pproc, "cases_averted_plot")
    
  )
  
  return(plots)
}

# model function ####
models_func <- list(
  variables = variables,
  simulator = simulate,
  processor = postproc,
  plotter = plotter
)

# simulate model
model_sim = models_func$simulator()
model_sim

# post-process model results
model_pproc = models_func$processor(model_sim)
model_pproc 

# plot model results
plot_list = models_func$plotter(model_pproc)
plot_list
