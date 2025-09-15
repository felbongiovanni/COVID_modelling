lognormal_mean <- function(mu, sigma) {
  return(log(mu) - sigma^2/2)
}



prep_stan_data_bp_model <- function(individual, dat, protein, ig) {
  # Filer data per individual and antigen
  df <- dat %>% filter(Record.ID == individual, IgClass %in% c(ig, NA)) %>% select(Spike, NP, RBD, Date.of.blood.sample.collection, When.did.you.receive.your.diagnosis.., Dose.1, Dose.2, Dose.3.., Dose.4, Dose.5, days, type, Event.Name)
  df <- df %>% arrange(Date.of.blood.sample.collection)
  
  # Extract all vaccination and diagnosis dates during follow-up
  tau <- as.numeric(as.Date(c(df$Dose.1, df$Dose.2, df$Dose.3.., df$Dose.4, df$Dose.5, df$When.did.you.receive.your.diagnosis..)) - min(df$Date.of.blood.sample.collection, na.rm = T))
  tau <- sort(tau[!is.na(tau) & tau >= 0])
  
  # Extract any diagnosis dates for natural infections before study enrollment
  dd <- as.numeric(as.Date(df$When.did.you.receive.your.diagnosis..) - min(df$Date.of.blood.sample.collection, na.rm = T))
  dd <- abs(dd[!is.na(dd) & dd < 0])
  
  # Extract all natural infections that occurred during follow-up
  nat_inf <- as.numeric(as.Date(c(df$When.did.you.receive.your.diagnosis..)) - min(df$Date.of.blood.sample.collection, na.rm = T))
  nat_inf <- sort(nat_inf[!is.na(nat_inf) & nat_inf > 0]) 
  
  # Drop any NAs
  df <- df %>% drop_na(protein) %>% filter(!is.na(Date.of.blood.sample.collection))
  
  if (unique(df$type) == "recovered") { # if individual is in the recovered cohort
    nat_inf <- c(0, nat_inf + dd) # day 0 is from first exposure before study enrollment
    time <- df$days + dd
    tau <- c(0, tau + dd)
    prior_inf_data <- which(time < tau[2]) # want only samples before any natural infections during follow up
    if (length(tau) < 2) { # if no natural infections occurred 
      prior_inf_data <- 1:length(time)
    }
    if (length(nat_inf) > 1) { # if a natural infection occurred, filter for data only before that
      non_reinf_data <- which(time < nth(nat_inf,2))
    } else {
      non_reinf_data <- 1:length(time)
    }
  } else { # if individual is in the naive cohort
    nat_inf <- nat_inf
    time <- df$days 
    tau <- tau
    prior_inf_data <- 0
    if (length(nat_inf) > 0) { # if a natural infection occurred during follow-up, filter for data before that
      non_reinf_data <- which(time < min(nat_inf, na.rm = T))
    } else {
      non_reinf_data <- 1:length(time)
    }
  }
  
  # Extract data for samples before any natural infections during follow up (if they occurred)
  prior_inf_date <- df$When.did.you.receive.your.diagnosis..
  new_time <- time[non_reinf_data]
  new_tau <- tau[tau <= max(new_time)]
  antibodies <- df[[protein]]
  new_ab <- antibodies[non_reinf_data]
  
  # Search for any asymptomatic infections that may not have been captured
  NP_ab <- df$NP[non_reinf_data]
  if (protein == "Spike") {
    other_ab <- df$RBD[non_reinf_data]
  } else {
    other_ab <- df$Spike[non_reinf_data]
  }
  
  if (unique(df$type) == "recovered") {
    asymp_inf <- NP_ab > lag(NP_ab) & (new_ab > lag(new_ab) & other_ab > lag(other_ab)) & NP_ab >= rec_threshold
  } else {
    asymp_inf <- NP_ab > lag(NP_ab) & (new_ab > lag(new_ab) & other_ab > lag(other_ab)) & NP_ab >= naive_threshold
  }
  first_increase <- ifelse(any(asymp_inf), which(asymp_inf)[1], NA)

  # Filter data to remove any potential asymptomatic infections
  if (!is.na(first_increase)) {
    if (!df$Event.Name[first_increase] %in% c("V1", "V2", "V3", "B1", "B2", "BB1", "BB2")) {
      new_ab <- new_ab[1:(first_increase-1)]
      new_time <- new_time[1:(first_increase-1)]
      new_tau <- new_tau[new_tau <= max(new_time)]
    }
  }
  
  return(list(ID = individual, no_reinf_y = new_ab, no_reinf_t = new_time, no_reinf_tau = new_tau, size_tau = length(new_tau), size_no_reinf = length(new_time), prior_inf_data = length(prior_inf_data), prior_inf_time = time[prior_inf_data], prior_inf_ab = antibodies[prior_inf_data], size_prior_inf = length(prior_inf_data), inf_date = prior_inf_date))
}



prep_stan_data_sp_model <- function(individual, dat, protein, ig) {
  # Filer data per individual and antigen
  df <- dat %>% filter(Record.ID == individual, IgClass %in% c(ig, NA)) %>% select(Spike, NP, RBD, Date.of.blood.sample.collection, When.did.you.receive.your.diagnosis.., Dose.1, Dose.2, Dose.3.., Dose.4, Dose.5, days, type, Event.Name)
  df <- df %>% arrange(Date.of.blood.sample.collection)
  
  # Extract all vaccination and diagnosis dates during follow-up
  tau <- as.numeric(as.Date(c(df$Dose.1, df$Dose.2, df$Dose.3.., df$Dose.4, df$Dose.5, df$When.did.you.receive.your.diagnosis..)) - min(df$Date.of.blood.sample.collection, na.rm = T))
  tau <- sort(tau[!is.na(tau) & tau >= 0])
  
  # Extract any diagnosis dates for natural infections before study enrollment
  dd <- as.numeric(as.Date(df$When.did.you.receive.your.diagnosis..) - min(df$Date.of.blood.sample.collection, na.rm = T))
  dd <- abs(dd[!is.na(dd) & dd < 0])
  
  # Extract all natural infections that occurred during follow-up
  nat_inf <- as.numeric(as.Date(c(df$When.did.you.receive.your.diagnosis..)) - min(df$Date.of.blood.sample.collection, na.rm = T))
  nat_inf <- sort(nat_inf[!is.na(nat_inf) & nat_inf > 0]) 
  
  # Drop any NAs
  df <- df %>% drop_na(protein) %>% filter(!is.na(Date.of.blood.sample.collection))
  
  if (unique(df$type) == "recovered") { # if individual is in the recovered cohort
    nat_inf <- c(0, nat_inf + dd) # day 0 is from first exposure before study enrollment
    time <- df$days + dd
    tau <- c(0, tau + dd)
    prior_inf_data <- which(time < tau[2]) # want only samples before any natural infections during follow up
    if (length(tau) < 2) { # if no natural infections occurred 
      prior_inf_data <- 1:length(time)
    }
    if (length(nat_inf) > 1) { # if a natural infection occurred, filter for data only before that
      non_reinf_data <- which(time < nth(nat_inf,2))
    } else {
      non_reinf_data <- 1:length(time)
    }
  } else { # if individual is in the naive cohort
    nat_inf <- nat_inf
    time <- df$days 
    tau <- tau
    prior_inf_data <- 0
    if (length(nat_inf) > 0) { # if a natural infection occurred during follow-up, filter for data before that
      non_reinf_data <- which(time < min(nat_inf, na.rm = T))
    } else {
      non_reinf_data <- 1:length(time)
    }
  }
  
  # Extract data for samples before any natural infections during follow up (if they occurred)
  prior_inf_date <- df$When.did.you.receive.your.diagnosis..
  new_time <- time[non_reinf_data]
  new_tau <- tau[tau <= max(new_time)]
  antibodies <- df[[protein]]
  new_ab <- antibodies[non_reinf_data]
  
  # Search for any asymptomatic infections that may not have been captured
  NP_ab <- df$NP[non_reinf_data]
  if (protein == "Spike") {
    other_ab <- df$RBD[non_reinf_data]
  } else {
    other_ab <- df$Spike[non_reinf_data]
  }
  
  if (unique(df$type) == "recovered") {
    asymp_inf <- NP_ab > lag(NP_ab) & (new_ab > lag(new_ab) & other_ab > lag(other_ab)) & NP_ab >= rec_threshold
  } else {
    asymp_inf <- NP_ab > lag(NP_ab) & (new_ab > lag(new_ab) & other_ab > lag(other_ab)) & NP_ab >= naive_threshold
  }
  first_increase <- ifelse(any(asymp_inf), which(asymp_inf)[1], NA)
  
  # Filter data to remove any potential asymptomatic infections
  if (!is.na(first_increase)) {
    if (!df$Event.Name[first_increase] %in% c("V1", "V2", "V3", "B1", "B2", "BB1", "BB2")) {
      new_ab <- new_ab[1:(first_increase-1)]
      new_time <- new_time[1:(first_increase-1)]
      new_tau <- new_tau[new_tau <= max(new_time)]
    }
  }
  
  segments <- split(new_time, findInterval(new_time, new_tau, left.open = TRUE))
  ab_segments <- list()
  time_segments <- list()
  og_final_times <- list()
  
  for (i in 1:length(segments)) {
    abs_i <- new_ab[which(new_time %in% unlist(segments[i]))]
    times_i <- new_time[which(new_time %in% unlist(segments[i]))]
    final_times_i <- times_i[which.max(abs_i):length(abs_i)]
    final_abs_i <- abs_i[which.max(abs_i):length(abs_i)]
    if (length(final_abs_i) > 1) {
      og_final_times <- append(og_final_times, list(final_times_i))
      final_times_i <- final_times_i - final_times_i[1]
      ab_segments <- append(ab_segments, list(final_abs_i))
      time_segments <- append(time_segments, list(final_times_i))
    }
  }
  
  return(list(segments = segments, ab_segments = ab_segments, time_segments = time_segments, size_seg = length(ab_segments), size_all_seg = lengths(ab_segments), og_final_times = og_final_times, ID = individual))
}



calculate_antibody_boosts <- function(stan_output, stan_input) {
  df <- tibble(median_boost = colMedians(stan_output$beta), id = stan_input$ids, exposure = 1)  
  pos <- 1
  
  for (i in 1:length(stan_input$ids)) {
    for (j in 1:stan_input$size_multi[i]) {
      boost <- median(stan_output$beta[,i] * stan_output$multiplier[,pos])
      df <- rbind(df, tibble(median_boost = boost, id = stan_input$ids[i], exposure = j+1))
      pos <- pos + 1
    }
  }
  df <- df %>% arrange(match(id, stan_input$ids))
  
  return(df)
}  









