# Load oxygen files
directory <- "~/Desktop/Projects/Seagrass/Data/Oxygen"
files <- list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)

require(tidyverse)
O2 <- files %>%
  map(~ read.csv(., skip = 1, header = TRUE) %>%
        drop_na(Value) %>%
        mutate(delta_t = as.numeric(delta_t),
               delta_t_c = delta_t - mean(delta_t))) %>%
  set_names(str_remove(basename(files), "\\.csv$") %>% make.names)
str(O2$X220915_B1)

require(tidybayes)
O2_list <- O2 %>%
  map(~ select(., Value, delta_t_c) %>%
        compose_data())

O2_stan <- "
data{
  int n;
  vector<lower=0>[n] Value;
  vector[n] delta_t_c;
}

parameters{
  real<lower=0> alpha;
  real beta;
  real<lower=0> sigma;
}

model{
  // Priors (from Rose et al. 2012, doi 10.5194/os-8-545-2012 and Woo & Pattiaratchi 2008, doi 10.1016/j.dsr.2008.05.005)
  alpha ~ gamma( 227.9194^2 / 20^2, 227.9194 / 20^2 ); // reparameterised with mean and sd
  beta ~ normal( 0, 5 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = alpha + beta * delta_t_c[i];
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}
"

require(cmdstanr)
O2_mod <- cmdstan_model(stan_file = write_stan_file(code = O2_stan))

O2_samples <- O2_list %>%
  map(~ O2_mod$sample(data = .,
                      seed = 100,
                      chains = 8,
                      parallel_chains = parallel::detectCores(),
                      iter_warmup = 1e4,
                      iter_sampling = 1e4))

# check Rhat, effective sample size and chains
O2_summary <- O2_samples %>%
  map(~ .$summary())

O2_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

O2_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# decent Rhat and effective sample size

O2_draws <- O2_samples %>%
  map(~ .$draws(format = "df"))

require(bayesplot)
require(patchwork)
O2_draws %>% 
  map(~ mcmc_rank_overlay(.)) %>%
  wrap_plots() %>%
  ggsave(filename = "rank.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

O2_draws %>% 
  sample(6) %>%
  map(~ mcmc_pairs(., pars = c("alpha", "beta"))) %>%
  wrap_plots() %>%
  ggsave(filename = "pairs.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 50, height = 50, units = "cm")


# Prior-posterior comparison
O2_prior_posterior <- O2_draws %>%
  map(~ spread_draws(., alpha, beta, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1), 
          beta_prior = rnorm(length(.draw), 0, 5),
          alpha_prior = rgamma(length(.draw), 227.9194^2 / 20^2, 227.9194 / 20^2)
        )
      )

O2_prior_posterior %>%
  map(~ pivot_longer(., cols = c("beta", "beta_prior", "alpha", "alpha_prior"),
                        names_to = c("parameter", "distribution"),
                        names_sep = "_", # this will produce NAs and throw a warning message 
                        values_to = "samples") %>%
        mutate(parameter = fct(parameter),
               distribution = fct(ifelse(is.na(distribution), # here the NAs are dealt with
                                         "posterior", distribution))) %>%
        ggplot(aes(samples, fill = distribution)) +
          facet_wrap(~ parameter, scales = "free", nrow = 1) +
          geom_density(colour = NA, alpha = 0.5) +
          theme_minimal() +
          theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() %>%
  ggsave(filename = "ppcomp.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

O2_mu <- O2_prior_posterior %>% 
  map2(O2,
       ~ select(.x, alpha, beta) %>%
         expand_grid(delta_t_c = seq(.y %>% pull(delta_t_c) %>% min(),
                                     .y %>% pull(delta_t_c) %>% max(),
                                     length.out = 50)) %>%
         mutate(mu = alpha + beta * delta_t_c) %>%
         select(mu, delta_t_c)
       )

O2_mu_summary <- O2_mu %>%
  map(~ group_by(., delta_t_c) %>%
        mean_qi(mu, .width = c(.5, .8, .9))
      )
         
O2 %>%
  map2(O2_mu_summary,
       ~ ggplot() +
          geom_point(data = .x, aes(delta_t_c, Value)) +
          geom_line(data = .y, aes(delta_t_c, mu)) +
          geom_ribbon(data = .y, aes(delta_t_c, ymin = .lower, ymax = .upper,
                                     alpha = factor(.width))) +
          scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
          theme_minimal()) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "fit.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")
# fit looks good

O2_initial <- O2_mu %>%
  map(~ slice_min(., delta_t_c))

O2_initial %>%
  map(~ ggplot(., aes(mu)) +
          geom_density() +
          theme_minimal()) %>%
  wrap_plots() %>%
  ggsave(filename = "initial.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")


V <- read.csv("~/Desktop/Projects/Seagrass/Data/Volume.csv") %>%
  mutate(Vs = (Volume - mean(Volume)) / sd(Volume)) # standardise Volume

V_stan <- "
data{
  int n;
  vector[n] Vs;
}

parameters{
  real Vsmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Vsmu ~ normal( 0, 1 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Vsmu;
  }

  // Likelihood
  Vs ~ normal( mu , sigma );
}
"

V_mod <- cmdstan_model(stan_file = write_stan_file(code = V_stan))

V_samples <- V_mod$sample(data = compose_data(V),
                          seed = 100,
                          chains = 8,
                          parallel_chains = parallel::detectCores(),
                          iter_warmup = 1e4,
                          iter_sampling = 1e4)

# check Rhat, effective sample size and chains
V_summary <- V_samples$summary()

V_summary %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

V_draws <- V_samples$draws(format = "df")

V_draws %>% mcmc_rank_overlay()

V_draws %>% mcmc_pairs()

# Prior-posterior comparison
V_prior_posterior <- V_draws %>%
  spread_draws(Vsmu, sigma) %>%
  ungroup() %>%
  mutate(
    sigma_prior = rexp(length(.draw), 1), 
    Vsmu_prior = rnorm(length(.draw), 0, 1)
    )

V_prior_posterior %>%
  pivot_longer(., cols = c("Vsmu", "Vsmu_prior", "sigma", "sigma_prior"),
                     names_to = c("parameter", "distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message 
                     values_to = "samples") %>%
  mutate(parameter = fct(parameter),
         distribution = fct(ifelse(is.na(distribution), # here the NAs are dealt with
                                   "posterior", distribution))) %>%
  ggplot(aes(samples, fill = distribution)) +
    facet_wrap(~ parameter, scales = "free", nrow = 1) +
    geom_density(colour = NA, alpha = 0.5) +
    theme_minimal() +
    theme(legend.position = "top", legend.justification = 0)

require(magrittr)
V_prior_posterior %<>%
  mutate(Vmu = V %$% (Vsmu * sd(Volume) + mean(Volume))) %>%
  select(Vmu)


O2 %>%
  map(~ summarise(., Tmean = mean(Temp),
                  Tsd = sd(Temp),
                  Pmean = mean(Pressure),
                  Psd = sd(Pressure))) %>%
  bind_rows() %>%
  filter(Psd == 0 | Tsd == 0)
# in several cases pressure did not vary at all across the incubation,
# so measurement error in pressure can not be included in the model

O2 %<>%
  map(~ mutate(., Ts = (Temp - mean(Temp)) / sd(Temp)))

T_list <- O2 %>%
  map(~ select(., Ts) %>%
        compose_data())

T_stan <- "
data{
  int n;
  vector[n] Ts;
}

parameters{
  real Tsmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Tsmu ~ normal( 0, 1 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Tsmu;
  }

  // Likelihood
  Ts ~ normal( mu , sigma );
}
"

T_mod <- cmdstan_model(stan_file = write_stan_file(code = T_stan))

T_samples <- T_list %>%
  map(~ T_mod$sample(data = .,
                     seed = 100,
                     chains = 8,
                     parallel_chains = parallel::detectCores(),
                     iter_warmup = 1e4,
                     iter_sampling = 1e4))

# check Rhat, effective sample size and chains
T_summary <- T_samples %>%
  map(~ .$summary())

T_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)

# no Rhat above 1.001

T_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# decent Rhat and effective sample size

T_draws <- T_samples %>%
  map(~ .$draws(format = "df"))

T_draws %>% 
  map(~ mcmc_rank_overlay(.)) %>%
  wrap_plots() %>%
  ggsave(filename = "Trank.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

T_draws %>% 
  sample(6) %>%
  map(~ mcmc_pairs(., pars = c("Tsmu", "sigma"))) %>%
  wrap_plots() %>%
  ggsave(filename = "Tpairs.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 50, height = 50, units = "cm")


# Prior-posterior comparison
T_prior_posterior <- T_draws %>%
  map(~ spread_draws(., Tsmu, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1), 
          Tsmu_prior = rnorm(length(.draw), 0, 1)
        )
  )

T_prior_posterior %>%
  map(~ pivot_longer(., cols = c("Tsmu", "Tsmu_prior", "sigma", "sigma_prior"),
                     names_to = c("parameter", "distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message 
                     values_to = "samples") %>%
        mutate(parameter = fct(parameter),
               distribution = fct(ifelse(is.na(distribution), # here the NAs are dealt with
                                         "posterior", distribution))) %>%
        ggplot(aes(samples, fill = distribution)) +
        facet_wrap(~ parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() %>%
  ggsave(filename = "Tppcomp.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

T_prior_posterior %<>%
  map2(O2,
       ~ mutate(.x, Tmu = .y %$% (Tsmu * sd(Temp) + mean(Temp))) %>%
         select(Tmu)
       )

# combine posterior distributions into one tibble

P <- O2_prior_posterior %>%
  map2(O2_initial,
       ~ bind_cols(.x, .y)) %>%
  map2(T_prior_posterior,
       ~ bind_cols(.x, .y, V_prior_posterior)) %>%
  map2(O2 %>% 
         map(~ summarise(., Salinity = mean(Salinity),
                         Pressure = mean(Pressure))),
       ~ bind_cols(.x, .y)) %>%
  map2(O2 %>% 
         map(~ select(., Date, Time) %>%
               slice(1)),
       ~ bind_cols(.x, .y)) %>%
  imap(~ mutate(.x, ID = .y)) %>%
  bind_rows()

# clean up, build explanatory variables from ID

P %<>%
  select(.draw, ID, Date, Time, beta, mu, Tmu, Vmu, Salinity, Pressure) %>%
  rename(O2 = mu, Temperature = Tmu, Volume = Vmu) %>%
  mutate(Date = Date %>% mdy(),
         Time = Time %>% hms(),
         Species = case_when(
           str_extract(ID, "(?<=_)\\p{L}(?=\\d*)") == "B" ~ "Blank",
           str_extract(ID, "(?<=_)\\p{L}(?=\\d*)") == "A" ~ "Amphibolis antarctica",
           str_extract(ID, "(?<=_)\\p{L}(?=\\d*)") == "H" ~ "Halophila ovalis"
           ),
         Leaf = str_extract(ID, "(?<=_[A-Za-z])\\d") %>% as.numeric())

# correct photosynthesis for blank, volume and mass (µM min^-1 to µmol g^-1 h^-1)
M <- read.csv("~/Desktop/Projects/Seagrass/Data/Mass.csv")
P %<>%
  filter(Species != "Blank") %>%
  left_join(P %>% filter(Species == "Blank") %>%
              select(.draw, Date, Time, beta) %>%
              rename(blank = beta),
            by = c(".draw", "Date", "Time"),
            relationship = "many-to-one") %>%
  left_join(M %>% mutate(Date = Date %>% dmy()), 
            by = c("Date", "Species", "Leaf"),
            relationship = "many-to-one") %>%
  mutate(Pm = (beta - blank) * Volume * 60 / Mass,
         Pl = if_else(Species == "Halophila ovalis", 
                      (beta - blank) * Volume * 60 / 10,
                      (beta - blank) * Volume * 60)) %>%
  select(-c(beta, Volume, blank))

# calculate detrital age (days) from dates
P %<>%
  group_by(Species) %>%
  mutate(Day = min(Date) %--% Date %>%
           time_length("day"))

# summarise P for modelling purposes
P_summary <- P %>%
  group_by(ID, Day, Species, Leaf) %>%
  summarise(Pm_mean = mean(Pm),
            Pm_sd = sd(Pm),
            Pl_mean = mean(Pl),
            Pl_sd = sd(Pl),
            O2_mean = mean(O2),
            O2_sd = sd(O2),
            T_mean = mean(Temperature),
            T_sd = sd(Temperature),
            P_mean = mean(Pressure),
            S = mean(Salinity),
            M = mean(Mass),
            n = length(.draw)) %>%
  ungroup()

P_summary %>%
  mutate(Incubation = if_else(nchar(ID) == 12, "Second", "First")) %>%
  ggplot() +
    geom_pointrange(aes(Day, Pm_mean, 
                        ymin = Pm_mean - Pm_sd, 
                        ymax = Pm_mean + Pm_sd,
                        colour = Species,
                        shape = Incubation)) +
    geom_smooth(aes(Day, Pm_mean, colour = Species),
                se = F) +
    theme_minimal()

P_summary %>%
  mutate(Incubation = if_else(nchar(ID) == 12, "Second", "First")) %>%
  ggplot() +
    geom_pointrange(aes(Day, Pl_mean, 
                        ymin = Pl_mean - Pl_sd, 
                        ymax = Pl_mean + Pl_sd,
                        colour = Species,
                        shape = Incubation)) +
    geom_smooth(aes(Day, Pl_mean, colour = Species),
                se = F) +
    facet_wrap(~Species, scales = "free") +
    theme_minimal()

P_summary %>%
  filter(nchar(ID) == 10) %>%
  ggplot() +
  geom_pointrange(aes(Day, Pm_mean, 
                      ymin = Pm_mean - Pm_sd, 
                      ymax = Pm_mean + Pm_sd,
                      colour = Species)) +
  geom_smooth(aes(Day, Pm_mean, colour = Species),
              se = F) +
  theme_minimal()

P_summary %>%
  filter(nchar(ID) == 10) %>%
  ggplot() +
  geom_pointrange(aes(Day, Pl_mean, 
                      ymin = Pl_mean - Pl_sd, 
                      ymax = Pl_mean + Pl_sd,
                      colour = Species)) +
  geom_smooth(aes(Day, Pl_mean, colour = Species),
              se = F) +
  facet_wrap(~Species, scales = "free") +
  theme_minimal()

# clearly the second incubation skews the data at 14 d, perhaps because of the higher initial O2
# -> proceed with only the first incubation
P %<>%
  filter(nchar(ID) == 10)
P_summary %<>%
  filter(nchar(ID) == 10)

# clean up
rm(O2, O2_draws, O2_initial, O2_list, O2_mod, O2_mu, O2_mu_summary,
   O2_prior_posterior, O2_samples, O2_summary, T_draws, T_list, T_mod,
   T_prior_posterior, T_samples, T_summary, V, V_draws, V_mod, V_samples,
   V_prior_posterior, V_summary, directory, files, O2_stan, T_stan, V_stan)

# Prior simulation
Prior <- read.csv("~/Desktop/Projects/Seagrass/Data/Prior.csv")

# calculate as species-specific means and sds
Prior %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Fl = Flux * Leafmass) %>%
  group_by(Species, Variable) %>%
  summarise(Fm_mean = mean(Flux),
            Fm_sd = sd(Flux),
            Fm_median = median(Flux),
            Fl_mean = mean(Fl),
            Fl_sd = sd(Fl),
            Fl_median = median(Fl),
            n = length(Flux))

# calculate overall mean and sd 
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Fl = Flux * Leafmass) %>%
  summarise(Pm_mean = mean(Flux),
            Pm_sd = sd(Flux),
            Pm_median = median(Flux),
            Pl_mean = mean(Fl),
            Pl_sd = sd(Fl),
            Pl_median = median(Fl),
            n = length(Flux))

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Fl = Flux * M %>% 
           mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                     Mass / 10, Mass)) %>%
           pull(Leafmass) %>% mean()) %>%
  summarise(Rm_mean = mean(Flux),
            Rm_sd = sd(Flux),
            Rm_median = median(Flux),
            Rl_mean = mean(Fl),
            Rl_sd = sd(Fl),
            Rl_median = median(Fl),
            n = length(Flux))

# based on the mean-median comparison photosynthesis and 
# detrital respiration are right-skewed

# visualise distributions
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Pl = Flux * Leafmass) %>%
  rename(Pm = Flux) %>%
  pivot_longer(cols = c(Pl, Pm), names_to = "Metric", values_to = "Flux") %>%
  ggplot(aes(Flux, fill = Metric)) +
    geom_density(colour = NA, alpha = 0.5) +
    geom_vline(xintercept = 35.2173, colour = "#00bfc4") + # calculated means
    geom_vline(xintercept = 3.165822, colour = "#f8766d") +
    facet_wrap(~ Metric, scales = "free") +
    theme_minimal()

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Rl = Flux * M %>% 
           mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                     Mass / 10, Mass)) %>%
           pull(Leafmass) %>% mean()) %>%
  rename(Rm = Flux) %>%
  pivot_longer(cols = c(Rl, Rm), names_to = "Metric", values_to = "Flux") %>%
  ggplot(aes(Flux, fill = Metric)) +
    geom_density(colour = NA, alpha = 0.5) +
    geom_vline(xintercept = 5.355247, colour = "#00bfc4") + # calculated means
    geom_vline(xintercept = 1.692589, colour = "#f8766d") +
    facet_wrap(~ Metric, scales = "free") +
    theme_minimal()
# the means seem to overestimate in all cases since the distributions are right-skewed
# -> better go with an approximate mode

# prior simulation for mass-based estimates
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 35.2173, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 35.2173^2 / 5^2, rate = 35.2173 / 5^2)), aes(x)) + # arbitrary sd
    theme_minimal()
# overlaying a normal distribution with a mean of the empirical mean doesn't really work
# however, I also don't want to use a gamma distribution because that would exclude
# negative values, which are theoretically possible but highly unlikely
# -> shift and widen distribution to capture the most probable peak while maintaining positive values
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 35.2173, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 20^2 / 10^2, rate = 20 / 10^2)), aes(x)) +
    theme_minimal()
# looks better

# for respiration I am using a gamma distribution because negative values are 
# not possible once photosynthesis has ceased
Prior %>%
  filter(Variable == "Detrital respiration") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 5.355247, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 5.355247^2 / 5^2, rate = 5.355247 / 5^2)), aes(x)) + # arbitrary sd
    theme_minimal()
# already looks good but could be shifted a little more to capture the more probable peak

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  ggplot() +
  geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 5.355247, colour = "#00bfc4") +
  geom_density(data = tibble(x = rgamma(n = 1e3, shape = 5^2 / 5^2, rate = 5 / 5^2)), aes(x)) +
  # the lognormal alternative produces a very similar result:
  # geom_density(data = tibble(x = rlnorm(n = 1e3, meanlog = log( 2^2 / sqrt(2^2 + 1.5^2) ), 
  #                                       sdlog = sqrt( log(1 + 1.5^2 / 2^2) ) )), aes(x)) +
  theme_minimal()

Pm_prior <- 
  tibble(n = 1:1e3,
         alpha = rnorm(n = 1e3, mean = 30, sd = 10),
         tau = rgamma(n = 1e3, shape = 2^2 / 1.5^2, rate = 2 / 1.5^2), # reparameterised as mean and sd
         k = rgamma(n = 1e3, shape = 0.12^2 / 0.05^2, rate = 0.12 / 0.05^2), # see Fig. 2b in doi 10.1093/aob/mcad167
         mu = rnorm(n = 1e3, mean = 15, sd = 10)) %>% # mean(P_summary$Day) = 15
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pm_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
    geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
    geom_line(alpha = 0.05) +
   coord_cartesian(expand = F, clip = "off") +
   theme_minimal()

# prior simulation for leaf-based estimates
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Flux = Flux * Leafmass) %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 3.165822, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 3.165822^2 / 5^2, rate = 3.165822 / 5^2)), aes(x)) + # arbitrary sd
    theme_minimal()
# shift and narrow distribution to capture the most probable peak while maintaining positive values
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Flux = Flux * Leafmass) %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 3.165822, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 3^2 / 1^2, rate = 3 / 1^2)), aes(x)) +
    theme_minimal()
# looks better

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Flux = Flux * M %>% 
           mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                     Mass / 10, Mass)) %>%
           pull(Leafmass) %>% mean()) %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 1.692589, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 1.692589^2 / 2^2, rate = 1.692589 / 2^2)), aes(x)) + # arbitrary sd
    theme_minimal()
# already looks good but could be shifted a little more to capture the more probable peak

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Flux = Flux * M %>% 
           mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                     Mass / 10, Mass)) %>%
           pull(Leafmass) %>% mean()) %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 1.692589, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 1.7^2 / 1.7^2, rate = 1.7 / 1.7^2)), aes(x)) +
    theme_minimal()

Pl_prior <- 
  tibble(n = 1:1e3,
         alpha = rgamma(n = 1e3, shape = 3^2 / 3^2, rate = 3 / 3^2),
         tau = rgamma(n = 1e3, shape = 1.7^2 / 1.7^2, rate = 1.7 / 1.7^2), 
         k = rgamma(n = 1e3, shape = 0.12^2 / 0.1^2, rate = 0.12 / 0.1^2), 
         mu = rgamma(n = 1e3, shape = 17^2 / 10^2, rate = 17 / 10^2)) %>% 
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pl_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
    geom_hline(yintercept = P_summary %$% c(min(Pl_mean), 0, max(Pl_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()



# models for mass-based estimates
Pm_stan <- "
data{
  int n;
  int n_Species;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pm;
  
  // Species parameters
  vector[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector[n_Species] mu;
  
  // Pooled parameters
  real alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> k_mu;
  real mu_mu;
  
  real<lower=0> alpha_sigma;
  real<lower=0> tau_sigma;
  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;
  
  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Pooled priors
  alpha_mu ~ normal( 30, 10 );
  tau_mu ~ gamma( 2^2 / 1.5^2, 2 / 1.5^2 ); // reparameterised with mean and sd
  k_mu ~ gamma( 0.12^2 / 0.05^2, 0.12 / 0.05^2 );
  // lognormal alternative priors for tau and k produces similar results
  // tau_mu ~ lognormal( log( 2^2 / sqrt(2^2 + 1.5^2) ), sqrt( log(1 + 1.5^2 / 2^2) ) ); // reparameterised with mean and sd
  // k_mu ~ lognormal( log( 0.12^2 / sqrt(0.12^2 + 0.05^2) ), sqrt( log(1 + 0.05^2 / 0.12^2) ) );
  mu_mu ~ normal( 15, 10 );
  
  alpha_sigma ~ exponential( 1 );
  tau_sigma ~ exponential( 1 );
  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );
  
  // Species priors
  alpha ~ normal( alpha_mu, alpha_sigma );
  tau ~ normal( tau_mu, tau_sigma );
  k ~ normal( k_mu, k_sigma );
  mu ~ normal( mu_mu, mu_sigma );
  
  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );
  
  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] = 
    ( alpha[Species[i]] + tau[Species[i]] ) / 
    ( 1 + exp( k[Species[i]] *  ( Day[i] - mu[Species[i]] ) ) ) 
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu, Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}
"

Pm_mod <- cmdstan_model(stan_file = write_stan_file(code = Pm_stan))

Pm_samples <- Pm_mod$sample(data = P_summary %>% 
                              select(Pm_mean, Pm_sd, Day, Species) %>% 
                              compose_data(),
                            seed = 100,
                            chains = 8,
                            parallel_chains = parallel::detectCores(),
                            iter_warmup = 1e4,
                            iter_sampling = 1e4)

Pm_summary <- Pm_samples$summary()

Pm_draws <- Pm_samples$draws(format = "df")

Pm_draws %>% mcmc_rank_overlay()

Pm_draws %>% mcmc_pairs(pars = c("alpha[1]", "tau[1]", "k[1]", "mu[1]"))
Pm_draws %>% mcmc_pairs(pars = c("alpha[2]", "tau[2]", "k[2]", "mu[2]"))

# correlation between alpha and mu
# -> centre predictor variable Days

P_summary %<>%
  mutate(Day_c = Day - mean(Day))

Pm_prior_c <- 
  tibble(n = 1:1e3,
         alpha = rnorm(n = 1e3, mean = 30, sd = 10),
         tau = rgamma(n = 1e3, shape = 2^2 / 1.5^2, rate = 2 / 1.5^2),
         k = rgamma(n = 1e3, shape = 0.12^2 / 0.05^2, rate = 0.12 / 0.05^2), # see Fig. 2b in doi 10.1093/aob/mcad167
         mu = rnorm(n = 1e3, mean = 0, sd = 10)) %>% # note change here: normal centred on zero because Day is centred
  expand_grid(Day_c = P_summary %$% seq(min(Day_c), max(Day_c), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day_c - mu ) ) ) - tau)

Pm_prior_c %>%
  ggplot(aes(x = Day_c, y = P_mu, group = n)) +
    geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()


Pm_stan_c <- "
data{
  int n;
  int n_Species;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] Day_c;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pm;
  
  // Species parameters
  vector[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector[n_Species] mu;
  
  // Pooled parameters
  real alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> k_mu;
  real mu_mu;
  
  real<lower=0> alpha_sigma;
  real<lower=0> tau_sigma;
  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;
  
  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Pooled priors
  alpha_mu ~ normal( 30, 10 );
  tau_mu ~ gamma( 2^2 / 1.5^2, 2 / 1.5^2 ); // reparameterised with mean and sd
  k_mu ~ gamma( 0.12^2 / 0.05^2, 0.12 / 0.05^2 ); // reparameterised with mean and sd
  mu_mu ~ normal( 0, 10 );
  
  alpha_sigma ~ exponential( 1 );
  tau_sigma ~ exponential( 1 );
  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );
  
  // Species priors
  alpha ~ normal( alpha_mu, alpha_sigma );
  tau ~ normal( tau_mu, tau_sigma );
  k ~ normal( k_mu, k_sigma );
  mu ~ normal( mu_mu, mu_sigma );
  
  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );
  
  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] = 
    ( alpha[Species[i]] + tau[Species[i]] ) / 
    ( 1 + exp( k[Species[i]] *  ( Day_c[i] - mu[Species[i]] ) ) ) 
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu, Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}
"

Pm_mod_c <- cmdstan_model(stan_file = write_stan_file(code = Pm_stan_c))

Pm_samples_c <- Pm_mod_c$sample(data = P_summary %>% 
                                  select(Pm_mean, Pm_sd, Day_c, Species) %>% 
                                  compose_data(),
                                seed = 100,
                                chains = 8,
                                parallel_chains = parallel::detectCores(),
                                iter_warmup = 1e4,
                                iter_sampling = 1e4)

Pm_summary_c <- Pm_samples_c$summary()

Pm_draws_c <- Pm_samples_c$draws(format = "df")

Pm_draws_c %>% mcmc_rank_overlay()

Pm_draws_c %>% mcmc_pairs(pars = c("alpha[1]", "tau[1]", "k[1]", "mu[1]"))
Pm_draws_c %>% mcmc_pairs(pars = c("alpha[2]", "tau[2]", "k[2]", "mu[2]"))

# centring the predictor does not actually help in this non-linear case
# and there were more divergent transitions
# -> try non-centred parameterisation of the hierarchical model

Pm_stan_nc <- "
data{
  int n;
  int n_Species;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pm;
  
  // Species parameters as z-scores
  vector[n_Species] alpha_z;
  vector[n_Species] tau_z;
  vector[n_Species] k_z;
  vector[n_Species] mu_z;
  
  // Pooled parameters
  real alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> k_mu;
  real mu_mu;
  
  real<lower=0> alpha_sigma;
  real<lower=0> tau_sigma;
  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;
  
  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Pooled priors
  alpha_mu ~ normal( 30, 10 );
  tau_mu ~ gamma( 2^2 / 1.5^2, 2 / 1.5^2 ); // reparameterised with mean and sd
  k_mu ~ gamma( 0.12^2 / 0.05^2, 0.12 / 0.05^2 ); // reparameterised with mean and sd
  mu_mu ~ normal( 0, 10 );
  
  alpha_sigma ~ exponential( 1 );
  tau_sigma ~ exponential( 1 );
  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );
  
  // Species priors as z-scores
  alpha_z ~ normal( 0, 1 );
  tau_z ~ normal( 0, 1 );
  k_z ~ normal( 0, 1 );
  mu_z ~ normal( 0, 1 );
  
  // Meaningful species parameters calculated from z-scores
  vector[n_Species] alpha;
  vector[n_Species] tau;
  vector[n_Species] k;
  vector[n_Species] mu;
  
  alpha = alpha_z * alpha_sigma + alpha_mu;
  tau = tau_z * tau_sigma + tau_mu;
  k = k_z * k_sigma + k_mu;
  mu = mu_z * mu_sigma + mu_mu;
  
  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );
  
  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] = 
    ( alpha[Species[i]] + tau[Species[i]] ) / 
    ( 1 + exp( k[Species[i]] *  ( Day[i] - mu[Species[i]] ) ) ) 
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu, Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}

generated quantities{
  // Saved species parameters
  vector[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector[n_Species] mu;
  
  alpha = alpha_z * alpha_sigma + alpha_mu;
  tau = tau_z * tau_sigma + tau_mu;
  k = k_z * k_sigma + k_mu;
  mu = mu_z * mu_sigma + mu_mu;
}
"

Pm_mod_nc <- cmdstan_model(stan_file = write_stan_file(code = Pm_stan_nc))

Pm_samples_nc <- Pm_mod_nc$sample(data = P_summary %>% 
                                  select(Pm_mean, Pm_sd, Day, Species) %>% 
                                  compose_data(),
                                seed = 100,
                                chains = 8,
                                parallel_chains = parallel::detectCores(),
                                iter_warmup = 1e4,
                                iter_sampling = 1e4)
# lots of sampling outside of allowed non-negative probability space

Pm_summary_nc <- Pm_samples_nc$summary()

Pm_draws_nc <- Pm_samples_nc$draws(format = "df")

Pm_draws_nc %>% mcmc_rank_overlay() # fails because of NAs

Pm_draws_nc %>% mcmc_pairs(pars = c("alpha[1]", "tau[1]", "k[1]", "mu[1]")) # fails because of NAs
Pm_draws_nc %>% mcmc_pairs(pars = c("alpha[2]", "tau[2]", "k[2]", "mu[2]")) # fails because of NAs

# the non-centred parameterisation actually performs worse because 
# using z-scores forfeits control over the non-negativity of tau and k
# -> try changing the higher level priors to something other than normal,
# adjust prior uncertainty and sampling parameters

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 5.355247, colour = "#00bfc4") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 2^2 / 1.5^2, rate = 2 / 1.5^2)), aes(x)) +
    # the lognormal alternative produces a very similar result:
    geom_density(data = tibble(x = rlnorm(n = 1e3, meanlog = log( 2^2 / sqrt(2^2 + 1.5^2) ),
                                          sdlog = sqrt( log(1 + 1.5^2 / 2^2) ) )), aes(x)) +
    theme_minimal()

Pm_prior <- 
  tibble(n = 1:1e3,
         alpha = rnorm(n = 1e3, mean = 30, sd = 5), # halve uncertainty around alpha
         tau = rgamma(n = 1e3, shape = 2^2 / 1.5^2, rate = 2 / 1.5^2), # reduce uncertainty around tau
         k = rgamma(n = 1e3, shape = 0.12^2 / 0.05^2, rate = 0.12 / 0.05^2),
         mu = rnorm(n = 1e3, mean = 15, sd = 5)) %>% # halve uncertainty around mu
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pm_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
    geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()

# gamma and lognormal effectively create the same prior given the same mean and sd,
# but they can lead to different sampling behaviour, so need to be systematically compared
# -> compare gamma and lognormal distributions using model summaries
Pm_stan <- "
data{
  int n;
  int n_Species;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pm;
  
  // Species parameters
  vector[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector[n_Species] mu;
  
  // Pooled parameters
  real alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> k_mu;
  real mu_mu;
  
  real<lower=0> alpha_sigma;
  real<lower=0> tau_sigma;
  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;
  
  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Pooled priors
  alpha_mu ~ normal( 30, 5 );
  tau_mu ~ gamma( 2^2 / 1.5^2, 2 / 1.5^2 );
  k_mu ~ gamma( 0.12^2 / 0.05^2, 0.12 / 0.05^2 );
  mu_mu ~ normal( 15, 5 );
  
  alpha_sigma ~ exponential( 1 );
  tau_sigma ~ exponential( 1 );
  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );
  
  // Species priors
  alpha ~ normal( alpha_mu, alpha_sigma );
  tau ~ gamma( tau_mu^2 / tau_sigma^2 , tau_mu / tau_sigma^2 );
  k ~ gamma( k_mu^2 / k_sigma^2 , k_mu / k_sigma^2 );
  mu ~ normal( mu_mu, mu_sigma );
  
  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );
  
  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] = 
    ( alpha[Species[i]] + tau[Species[i]] ) / 
    ( 1 + exp( k[Species[i]] *  ( Day[i] - mu[Species[i]] ) ) ) 
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu, Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}
"

Pm_mod <- cmdstan_model(stan_file = write_stan_file(code = Pm_stan))

Pm_samples <- Pm_mod$sample(data = P_summary %>% 
                              select(Pm_mean, Pm_sd, Day, Species) %>% 
                              compose_data(),
                            seed = 100,
                            chains = 8,
                            parallel_chains = parallel::detectCores(),
                            iter_warmup = 1e4,
                            iter_sampling = 1e4,
                            max_treedepth = 12, # increased maximum tree depth based on warning
                            adapt_delta = 0.9) # increase target acceptance rate to reduce divergences

Pm_summary_gamma <- Pm_samples$summary()

Pm_stan <- "
data{
  int n;
  int n_Species;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pm;
  
  // Species parameters
  vector[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector[n_Species] mu;
  
  // Pooled parameters
  real alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> k_mu;
  real mu_mu;
  
  real<lower=0> alpha_sigma;
  real<lower=0> tau_sigma;
  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;
  
  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Pooled priors
  alpha_mu ~ normal( 30, 5 );
  tau_mu ~ lognormal( log( 2^2 / sqrt(2^2 + 1.5^2) ), sqrt( log(1 + 1.5^2 / 2^2) ) );
  k_mu ~ lognormal( log( 0.12^2 / sqrt(0.12^2 + 0.05^2) ), sqrt( log(1 + 0.05^2 / 0.12^2) ) );
  mu_mu ~ normal( 15, 5 );
  
  alpha_sigma ~ exponential( 1 );
  tau_sigma ~ exponential( 1 );
  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );
  
  // Species priors
  alpha ~ normal( alpha_mu, alpha_sigma );
  tau ~ lognormal( log( tau_mu^2 / sqrt(tau_mu^2 + tau_sigma^2) ), sqrt( log(1 + tau_sigma^2 / tau_mu^2) ) );
  k ~ lognormal( log( k_mu^2 / sqrt(k_mu^2 + k_sigma^2) ), sqrt( log(1 + k_sigma^2 / k_mu^2) ) );
  mu ~ normal( mu_mu, mu_sigma );
  
  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );
  
  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] = 
    ( alpha[Species[i]] + tau[Species[i]] ) / 
    ( 1 + exp( k[Species[i]] *  ( Day[i] - mu[Species[i]] ) ) ) 
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu, Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}
"

Pm_mod <- cmdstan_model(stan_file = write_stan_file(code = Pm_stan))

Pm_samples <- Pm_mod$sample(data = P_summary %>% 
                              select(Pm_mean, Pm_sd, Day, Species) %>% 
                              compose_data(),
                            seed = 100,
                            chains = 8,
                            parallel_chains = parallel::detectCores(),
                            iter_warmup = 1e4,
                            iter_sampling = 1e4,
                            max_treedepth = 12, # increased maximum tree depth based on warning
                            adapt_delta = 0.9) # increase target acceptance rate to reduce divergences

Pm_summary_lnorm <- Pm_samples$summary()

# compare model Rhat and effective sample sizes
Pm_summary_gamma %>%
  bind_rows(Pm_summary_lnorm, .id = "Distribution") %>%
  mutate(Distribution = if_else(Distribution == 1, "gamma", "lognormal")) %>%
  select(Distribution, variable, rhat, ess_bulk, ess_tail) %>%
  pivot_wider(names_from = Distribution, values_from = c(rhat, ess_bulk, ess_tail)) %>%
  ggplot(aes(rhat_gamma, rhat_lognormal)) +
    geom_point() +
    geom_abline(slope = 1) +
    xlim(1, 1.1) +
    theme_minimal()

Pm_summary_gamma %>%
  bind_rows(Pm_summary_lnorm, .id = "Distribution") %>%
  mutate(Distribution = if_else(Distribution == 1, "gamma", "lognormal")) %>%
  select(Distribution, variable, rhat, ess_bulk, ess_tail) %>%
  pivot_wider(names_from = Distribution, values_from = c(rhat, ess_bulk, ess_tail)) %>%
  ggplot(aes(ess_bulk_gamma, ess_bulk_lognormal)) +
    geom_point() +
    geom_abline(slope = 1) +
    theme_minimal()

Pm_summary_gamma %>%
  bind_rows(Pm_summary_lnorm, .id = "Distribution") %>%
  mutate(Distribution = if_else(Distribution == 1, "gamma", "lognormal")) %>%
  select(Distribution, variable, rhat, ess_bulk, ess_tail) %>%
  pivot_wider(names_from = Distribution, values_from = c(rhat, ess_bulk, ess_tail)) %>%
  ggplot(aes(ess_tail_gamma, ess_tail_lognormal)) +
    geom_point() +
    geom_abline(slope = 1) +
    theme_minimal()


# Pm_prior <- 
#   tibble(n = 1:1e3,
#          alpha = rnorm(n = 1e3, mean = 35, sd = 5), 
#          tau = rgamma(n = 1e3, shape = 5^2 / 5^2, rate = 5 / 5^2), 
#          k = rgamma(n = 1e3, shape = 0.2^2 / 0.2^2, rate = 0.2 / 0.2^2),
#          t0 = rnorm(n = 1e3, mean = 5, sd = 3)) %>%
#   expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
#   mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * Day - t0 ) ) - tau)
# 
# Pm_prior %>%
#   ggplot(aes(x = Day, y = P_mu, group = n)) +
#   geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
#   geom_line(alpha = 0.05) +
#   coord_cartesian(expand = F, clip = "off") +
#   theme_minimal()
# 
# Pm_stan <- "
# data{
#   int n;
#   int n_Species;
#   vector[n] Pm_mean;
#   vector[n] Pm_sd;
#   vector[n] Day;
#   array[n] int Species;
# }
# 
# parameters{
#   // Estimate of Pm for measurment error
#   vector[n] Pm;
#   
#   // Species parameters
#   vector[n_Species] alpha;
#   vector<lower=0>[n_Species] tau;
#   vector<lower=0>[n_Species] k;
#   
#   // Pooled parameters
#   real alpha_mu;
#   real<lower=0> tau_mu;
#   real<lower=0> k_mu;
#   real t0_mu;
#   
#   real<lower=0> alpha_sigma;
#   real<lower=0> tau_sigma;
#   real<lower=0> k_sigma;
#   
#   // Likelihood uncertainty parameter
#   real<lower=0> Pm_sigma;
# }
# 
# model{
#   // Pooled priors
#   alpha_mu ~ normal( 35, 5 );
#   tau_mu ~ gamma( 5^2 / 5^2, 5 / 5^2 );
#   k_mu ~ gamma( 0.2^2 / 0.2^2, 0.2 / 0.2^2 );
#   t0_mu ~ normal( 5, 3 );
#   
#   alpha_sigma ~ exponential( 1 );
#   tau_sigma ~ exponential( 1 );
#   k_sigma ~ exponential( 1 );
#   
#   // Species priors
#   alpha ~ normal( alpha_mu, alpha_sigma );
#   tau ~ gamma( tau_mu^2 / tau_sigma^2 , tau_mu / tau_sigma^2 );
#   k ~ gamma( k_mu^2 / k_sigma^2 , k_mu / k_sigma^2 );
#   
#   // Likelihood uncertainty prior
#   Pm_sigma ~ exponential( 1 );
#   
#   // Model
#   vector[n] Pm_mu;
#   for ( i in 1:n ) {
#     Pm_mu[i] = 
#     ( alpha[Species[i]] + tau[Species[i]] ) / 
#     ( 1 + exp( k[Species[i]] * Day[i] - t0_mu ) ) 
#     - tau[Species[i]];
#   }
# 
#   // Likelihood incorporating measurement error
#   Pm ~ normal( Pm_mu, Pm_sigma );
#   Pm_mean ~ normal( Pm , Pm_sd );
# }
# 
# generated quantities{
#   // Calculate mu from t0_mu and k
#   
#   // Species parameter
#   vector<lower=0>[n_Species] mu;
#   mu = t0_mu / k;
#   
#   // Pooled parameter
#   real<lower=0> mu_mu;
#   mu_mu = t0_mu / k_mu;
# }
# "


Pm_prior <- 
  tibble(n = 1:1e3,
         alpha = rgamma(n = 1e3, shape = 35^2 / 2^2, rate = 35 / 2^2), 
         tau = rgamma(n = 1e3, shape = 5^2 / 5^2, rate = 5 / 5^2), 
         k = rgamma(n = 1e3, shape = 0.12^2 / 0.1^2, rate = 0.12 / 0.1^2),
         mu = rgamma(n = 1e3, shape = 17^2 / 10^2, rate = 17 / 10^2)) %>% 
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pm_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
  geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
  geom_line(alpha = 0.05) +
  coord_cartesian(expand = F, clip = "off") +
  theme_minimal()

Pm_stan <- "
data{
  int n;
  int n_Species;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pm;

  // Species parameters
  vector<lower=0>[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector<lower=0>[n_Species] mu;

  // Pooled parameters
  real<lower=0> alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> k_mu;
  real<lower=0> mu_mu;

  real<lower=0> alpha_sigma;
  real<lower=0> tau_sigma;
  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;

  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Pooled priors
  alpha_mu ~ gamma( 35^2 / 2^2, 35 / 2^2 );
  tau_mu ~ gamma( 5^2 / 5^2, 5 / 5^2 );
  k_mu ~ gamma( 0.12^2 / 0.1^2, 0.12 / 0.1^2 );
  mu_mu ~ gamma( 17^2 / 10^2, 17 / 10^2 );

  alpha_sigma ~ exponential( 1 );
  tau_sigma ~ exponential( 1 );
  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );

  // Species priors
  alpha ~ gamma( alpha_mu^2 / alpha_sigma^2 , alpha_mu / alpha_sigma^2 );
  tau ~ gamma( tau_mu^2 / tau_sigma^2 , tau_mu / tau_sigma^2 );
  k ~ gamma( k_mu^2 / k_sigma^2 , k_mu / k_sigma^2 );
  mu ~ gamma( mu_mu^2 / mu_sigma^2 , mu_mu / mu_sigma^2 );

  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );

  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] =
    ( alpha[Species[i]] + tau[Species[i]] ) /
    ( 1 + exp( k[Species[i]] * ( Day[i] - mu[Species[i]] ) ) )
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu, Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}
"
Pm_mod <- cmdstan_model(stan_file = write_stan_file(code = Pm_stan))

Pm_samples <- Pm_mod$sample(data = P_summary %>%
                              select(Pm_mean, Pm_sd, Day, Species) %>%
                              compose_data(),
                            seed = 100,
                            chains = 8,
                            parallel_chains = parallel::detectCores(),
                            iter_warmup = 1e4,
                            iter_sampling = 1e4,
                            max_treedepth = 12, # increased maximum tree depth based on warning
                            adapt_delta = 0.9) # increase target acceptance rate to reduce divergences

# despite these adjustments the lowest percentage of divergent transitions attained is 8%
Pm_summary <- Pm_samples$summary()
Pm_summary %>%
  filter(rhat > 1.001) # 21 rhat above 1.001
Pm_summary %>%
  filter(rhat > 1.01) # 3 above 1.01

Pm_draws <- Pm_samples$draws(format = "df")

Pm_draws %>% mcmc_rank_overlay() # chains look good

Pm_draws %>% mcmc_pairs(pars = c("alpha[1]", "tau[1]", "k[1]", "mu[1]"))
Pm_draws %>% mcmc_pairs(pars = c("alpha[2]", "tau[2]", "k[2]", "mu[2]"))
# this is chosen as the optimal model

rm(M, Pl_prior, Pm_draws_c, Pm_draws_nc, Pm_mod_c, Pm_mod_nc,
   Pm_prior, Pm_prior_c, Pm_samples_c, Pm_samples_nc, Pm_summary_c,
   Pm_summary_nc, Pm_summary_gamma, Pm_summary_lnorm, Pm_stan_c, Pm_stan_nc)


# Prior-posterior comparison
Pm_posterior <- Pm_samples %>%
  gather_draws(alpha[Species], tau[Species], k[Species], mu[Species],
               alpha_mu, tau_mu, k_mu, mu_mu) %>%
  mutate(Species = case_when(
                    Species == 1 ~ "Amphibolis antarctica",
                    Species == 2 ~ "Halophila ovalis",
                    is.na(Species) ~ "Seagrasses"
                    ),
         Parameter = case_when(
                      .variable %in% c("alpha", "alpha_mu") ~ "alpha",
                      .variable %in% c("tau", "tau_mu") ~ "tau",
                      .variable %in% c("k", "k_mu") ~ "k",
                      .variable %in% c("mu", "mu_mu") ~ "mu"
                      ),
         Distribution = "Posterior"
         )

Pm_prior <- tibble(Species = c("Amphibolis antarctica", "Halophila ovalis", "Seagrasses") %>%
                     rep(each = 8e4),
                   .chain = 1:8 %>% rep(each = 1e4) %>% rep(3),
                   .iteration = 1:1e4 %>% rep(8*3),
                   .draw = 1:8e4 %>% rep(3),
                   alpha = rgamma(8e4, 35^2 / 2^2, 35 / 2^2) %>% rep(3),
                   tau = rgamma(8e4, 5^2 / 5^2, 5 / 5^2) %>% rep(3),
                   k = rgamma(8e4, 0.12^2 / 0.1^2, 0.12 / 0.1^2) %>% rep(3),
                   mu = rgamma(8e4, 17^2 / 10^2, 17 / 10^2) %>% rep(3),
                   Distribution = "Prior") %>%
  pivot_longer(cols = c("alpha", "tau", "k", "mu"), values_to = ".value", names_to = ".variable") %>%
  mutate(.variable = fct_relevel(.variable, "alpha", "tau", "k"),
         Parameter = .variable) %>%
  arrange(.variable)


Pm_prior_posterior <- Pm_posterior %>% bind_rows(Pm_prior)

Pm_prior_posterior %>%
  ggplot(aes(.value, fill = Distribution)) +
    facet_wrap(Species ~ Parameter, scales = "free") +
    geom_density(colour = NA, alpha = 0.5) +
    theme_minimal()

# change so that Distribution and Species are merged into one grouping variable
Pm_prior_posterior_merged <- Pm_prior_posterior %>%
  ungroup() %>%
  filter(Distribution == "Prior" & Species == "Seagrasses" |
           Distribution == "Posterior") %>%
  mutate(Group = if_else(Distribution == "Prior", "Prior", Species) %>%
           fct_relevel("Seagrasses", "Amphibolis antarctica", "Halophila ovalis", "Prior"),
         Parameter = fct_relevel(Parameter, "alpha", "tau", "k")) %>%
  select(-c(Distribution, .variable, Species))

# calculate P_mu from parameters
Pm_mu <- Pm_prior_posterior_merged %>%
  pivot_wider(names_from = Parameter, values_from = .value) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pm_mu_summary <- Pm_mu %>%
  group_by(Day, Group) %>%
  mean_qi(P_mu, .width = c(.5, .8, .9))


mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(0, .5, .2, 0), "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing = unit(1, "cm"),
                 text = element_text(family = "Futura"))

require(ggridges)
Fig_1a_top <- Pm_prior_posterior_merged %>%
  ggplot() +
    geom_density_ridges(aes(.value, Group, colour = Group, fill = Group),
                        quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.001, bandwidth = c(0.2, 0.2, 0.01, 0.2),
                        from = c(25, 0, 0, 0), to = c(50, 15, 0.4, 40)) +
    scale_colour_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                        labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                   expression(italic("Halophila ovalis")), "Prior"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                      labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                 expression(italic("Halophila ovalis")), "Prior"),
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(~ Parameter, scales = "free_x",
               labeller = labeller(Parameter = as_labeller(
                 c("alpha" = "italic('a')*' (µmol g'^-1*' h'^-1*')'",
                   "tau" = "italic('t')*' (µmol g'^-1*' h'^-1*')'",
                   "k" = "italic('k')*' (d'^-1*')'",
                   "mu" = "italic('µ')*' (d)'"),
                 label_parsed))
               ) +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())

Fig_1a_bottom <- ggplot() +
                  geom_hline(yintercept = 0) +
                  geom_violin(data = P,
                              aes(Day, Pm, fill = Species,
                                  colour = Species, group = ID),
                              alpha = 0.2, position = "identity", width = 2) +
                  # too much overplotting, so I reduced the prior to the 0.05 and 0.95 quantiles
                  geom_ribbon(data = Pm_mu_summary %>% filter(Group == "Prior", .width == 0.9),
                              aes(Day, ymin = .lower, ymax = .upper), fill = NA, colour = "#b5b8ba") +
                  geom_line(data = Pm_mu_summary %>% filter(Group != "Prior"),
                            aes(Day, P_mu, colour = fct_relevel(Group, "Halophila ovalis",
                                                                "Amphibolis antarctica"))) +
                  geom_ribbon(data = Pm_mu_summary %>% filter(Group != "Prior"),
                              aes(Day, ymin = .lower, ymax = .upper,
                                  fill = fct_relevel(Group, "Halophila ovalis",
                                                     "Amphibolis antarctica"),
                              alpha = factor(.width)), colour = NA) +
                  labs(y = expression(italic(P["max"])*" (µmol O"[2]*" g"^-1*" h"^-1*")"),
                       x = "Detrital age (d)") +
                  scale_colour_manual(values = c("#bdd268", "#4a7518", "#363538"),
                                      guide = "none") +
                  scale_fill_manual(values = c("#bdd268", "#4a7518", "#363538"),
                                    guide = "none") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
                  scale_y_continuous(breaks = seq(-10, 50, by = 10),
                                     labels = scales::label_number(style_negative = "minus")) +
                  coord_cartesian(xlim = c(-0.8, 35), ylim = c(-10, 50), expand = FALSE, clip = "off") +
                  mytheme

Fig_1a <- ( Fig_1a_top / Fig_1a_bottom ) +
          plot_layout(heights = c(0.3, 1)) +
          plot_annotation(tag_levels = list(c("A", ""))) &
          theme(plot.tag = element_text(family = "Futura", size = 20, face = "bold"))

ggsave(plot = Fig_1a, filename = "Fig_1a.pdf", device = cairo_pdf, path = "~/Desktop",
       width = 20, height = 12, units = "cm")


# leaf-based estimates

# calculate overall mean and sd 
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Mass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Fl = Flux * Leafmass) %>%
  summarise(Pm_mean = mean(Flux),
            Pm_sd = sd(Flux),
            Pm_median = median(Flux),
            Pl_mean = mean(Fl),
            Pl_sd = sd(Fl),
            Pl_median = median(Fl),
            n = length(Flux))

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Fl = Flux * M %>% 
           pull(Mass) %>% mean()) %>%
  summarise(Rm_mean = mean(Flux),
            Rm_sd = sd(Flux),
            Rm_median = median(Flux),
            Rl_mean = mean(Fl),
            Rl_sd = sd(Fl),
            Rl_median = median(Fl),
            n = length(Flux))

# prior simulation for leaf-based estimates
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              group_by(Species) %>%
              summarise(Mass = mean(Mass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Flux = Flux * Mass) %>%
  ggplot() +
  geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 15.85307, colour = "#00bfc4") +
  geom_density(data = tibble(x = rgamma(n = 1e3, shape = 16^2 / 5^2, rate = 16 / 5^2)), aes(x)) +
  theme_minimal()


Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Flux = Flux * M %>% 
           pull(Mass) %>% mean()) %>%
  ggplot() +
  geom_density(aes(Flux), fill = "#00bfc4", alpha = 0.5, colour = NA) +
  geom_vline(xintercept = 2.782084, colour = "#00bfc4") +
  geom_density(data = tibble(x = rgamma(n = 1e3, shape = 3^2 / 3^2, rate = 3 / 3^2)), aes(x)) + 
  theme_minimal()
# already looks good but could be shifted a little more to capture the more probable peak


Pl_prior <- 
  tibble(n = 1:1e3,
         alpha = rgamma(n = 1e3, shape = 16^2 / 5^2, rate = 16 / 5^2),
         tau = rgamma(n = 1e3, shape = 3^2 / 3^2, rate = 3 / 3^2), 
         k = rgamma(n = 1e3, shape = 0.12^2 / 0.1^2, rate = 0.12 / 0.1^2), 
         mu = rgamma(n = 1e3, shape = 17^2 / 10^2, rate = 17 / 10^2)) %>% 
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pl_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
  geom_hline(yintercept = P_summary %$% c(min(Pl_mean), 0, max(Pl_mean))) +
  geom_line(alpha = 0.05) +
  coord_cartesian(expand = F, clip = "off") +
  theme_minimal()

Pl_stan <- "
data{
  int n;
  int n_Species;
  vector[n] Pl_mean;
  vector[n] Pl_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pl;

  // Species parameters
  vector<lower=0>[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector<lower=0>[n_Species] mu;

  // Pooled parameters
  real<lower=0> k_mu;
  real<lower=0> mu_mu;

  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;
  
  // Likelihood uncertainty parameter
  real<lower=0> Pl_sigma;
}

model{
  // Pooled priors
  k_mu ~ gamma( 0.12^2 / 0.1^2 , 0.12 / 0.1^2 );
  mu_mu ~ gamma( 17^2 / 10^2 , 17 / 10^2 );

  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );

  // Species priors
  alpha ~ gamma( 16^2 / 5^2 , 16 / 5^2 );
  tau ~ gamma( 3^2 / 3^2 , 3 / 3^2 );
  k ~ gamma( k_mu^2 / k_sigma^2 , k_mu / k_sigma^2 );
  mu ~ gamma( mu_mu^2 / mu_sigma^2 , mu_mu / mu_sigma^2 );
  
  // Likelihood uncertainty prior
  Pl_sigma ~ exponential( 1 );

  // Model
  vector[n] Pl_mu;
  for ( i in 1:n ) {
    Pl_mu[i] =
    ( alpha[Species[i]] + tau[Species[i]] ) /
    ( 1 + exp( k[Species[i]] * ( Day[i] - mu[Species[i]] ) ) )
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pl ~ normal( Pl_mu, Pl_sigma );
  Pl_mean ~ normal( Pl , Pl_sd );
}
"

Pl_stan <- "
data{
  int n;
  int n_Species;
  vector[n] Pl_mean;
  vector[n] Pl_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pl;

  // Species parameters
  vector<lower=0>[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector<lower=0>[n_Species] mu;
  
  // Likelihood uncertainty parameter
  real<lower=0> Pl_sigma;
}

model{
  // Species priors
  alpha ~ gamma( 16^2 / 5^2 , 16 / 5^2 );
  tau ~ gamma( 3^2 / 3^2 , 3 / 3^2 );
  k ~ gamma( 0.12^2 / 0.1^2 , 0.12 / 0.1^2 );
  mu ~ gamma( 17^2 / 10^2 , 17 / 10^2 );
  
  // Likelihood uncertainty prior
  Pl_sigma ~ exponential( 1 );

  // Model
  vector[n] Pl_mu;
  for ( i in 1:n ) {
    Pl_mu[i] =
    ( alpha[Species[i]] + tau[Species[i]] ) /
    ( 1 + exp( k[Species[i]] * ( Day[i] - mu[Species[i]] ) ) )
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pl ~ normal( Pl_mu, Pl_sigma );
  Pl_mean ~ normal( Pl , Pl_sd );
}
"
Pl_mod <- cmdstan_model(stan_file = write_stan_file(code = Pl_stan))

Pl_samples <- Pl_mod$sample(data = P_summary %>%
                              select(Pl_mean, Pl_sd, Day, Species) %>%
                              mutate(Pl_mean = if_else(Species == "Halophila ovalis",
                                                       Pl_mean * 10, Pl_mean),
                                     Pl_sd = if_else(Species == "Halophila ovalis",
                                                     Pl_sd * 10, Pl_sd)) %>%
                              compose_data(),
                            seed = 100,
                            chains = 8,
                            parallel_chains = parallel::detectCores(),
                            iter_warmup = 1e4,
                            iter_sampling = 1e4,
                            max_treedepth = 12, 
                            adapt_delta = 0.9)

Pl_summary <- Pl_samples$summary()
Pl_summary %>%
  filter(rhat > 1.001) # 5 rhat above 1.001

Pl_summary %>%
  filter(rhat > 1.01) # none above 1.01

Pl_draws <- Pl_samples$draws(format = "df")

Pl_draws %>% mcmc_rank_overlay() # chains look good

Pl_draws %>% mcmc_pairs(pars = c("alpha[1]", "tau[1]", "k[1]", "mu[1]"))
Pl_draws %>% mcmc_pairs(pars = c("alpha[2]", "tau[2]", "k[2]", "mu[2]"))
# this is chosen as the optimal model

# Prior-posterior comparison
Pl_posterior <- Pl_samples %>%
  gather_draws(alpha[Species], tau[Species], k[Species], mu[Species]) %>%
  mutate(Species = case_when(
                    Species == 1 ~ "Amphibolis antarctica",
                    Species == 2 ~ "Halophila ovalis",
                    ),
         .value = if_else(Species == "Halophila ovalis" & 
                            .variable %in% c("alpha", "tau"),
                          .value / 10, .value),
         Distribution = "Posterior")

Pl_prior <- tibble(Species = c("Amphibolis antarctica", "Halophila ovalis") %>%
                     rep(each = 8e4),
                   .chain = 1:8 %>% rep(each = 1e4) %>% rep(2),
                   .iteration = 1:1e4 %>% rep(8*2),
                   .draw = 1:8e4 %>% rep(2),
                   alpha = rgamma(8e4, 16^2 / 5^2, 16 / 5^2) %>% rep(2),
                   tau = rgamma(8e4, 3^2 / 3^2, 3 / 3^2) %>% rep(2),
                   k = rgamma(8e4, 0.12^2 / 0.1^2, 0.12 / 0.1^2) %>% rep(2),
                   mu = rgamma(8e4, 17^2 / 10^2, 17 / 10^2) %>% rep(2),
                   Distribution = "Prior") %>%
  pivot_longer(cols = c("alpha", "tau", "k", "mu"), values_to = ".value", names_to = ".variable") %>%
  mutate(.variable = fct_relevel(.variable, "alpha", "tau", "k"),
         .value = if_else(Species == "Halophila ovalis" & 
                            .variable %in% c("alpha", "tau"),
                          .value / 10, .value)
         ) %>%
  arrange(.variable)


Pl_prior_posterior <- Pl_posterior %>% bind_rows(Pl_prior)

Pl_prior_posterior %>%
  ggplot(aes(.value, fill = Distribution)) +
    facet_wrap(Species ~ .variable, scales = "free") +
    geom_density(colour = NA, alpha = 0.5) +
    theme_minimal()

# calculate P_mu from parameters
Pl_mu <- Pl_posterior %>%
  ungroup() %>%
  select(-Distribution) %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pl_mu_summary <- Pl_mu %>%
  group_by(Day, Species) %>%
  mean_qi(P_mu, .width = c(.5, .8, .9))

Fig_1b_left <- ggplot() +
                  geom_hline(yintercept = 0) +
                  geom_violin(data = P %>%
                                filter(Species == "Halophila ovalis"),
                              aes(Day, Pl, group = ID), colour = "#bdd268",
                              fill = "#bdd268", alpha = 0.2,
                              position = "identity", width = 4) +
                  geom_line(data = Pl_mu_summary %>%
                              filter(Species == "Halophila ovalis"),
                            aes(Day, P_mu), colour = "#bdd268") +
                  geom_ribbon(data = Pl_mu_summary %>%
                                filter(Species == "Halophila ovalis"),
                              aes(Day, ymin = .lower, ymax = .upper,
                                  alpha = factor(.width)), fill = "#bdd268", colour = NA) +
                  labs(y = expression(italic(P["max"])*" (µmol O"[2]*" leaf"^-1*" h"^-1*")"),
                       x = "Detrital age (d)") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
                  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5),
                                     labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1))) +
                  coord_cartesian(xlim = c(-1.8, 35), ylim = c(0, 2.5), expand = FALSE, clip = "off") +
                  mytheme

Fig_1b_right <- ggplot() +
                  geom_hline(yintercept = 0) +
                  geom_violin(data = P %>%
                                filter(Species == "Amphibolis antarctica"),
                              aes(Day, Pl, group = ID), colour = "#4a7518",
                              fill = "#4a7518", alpha = 0.2,
                              position = "identity", width = 4) +
                  geom_line(data = Pl_mu_summary %>%
                              filter(Species == "Amphibolis antarctica"),
                            aes(Day, P_mu), colour = "#4a7518") +
                  geom_ribbon(data = Pl_mu_summary %>%
                                filter(Species == "Amphibolis antarctica"),
                              aes(Day, ymin = .lower, ymax = .upper,
                              alpha = factor(.width)), fill = "#4a7518", colour = NA) +
                  labs(y = expression(italic(P["max"])*" (µmol O"[2]*" leaf"^-1*" h"^-1*")"),
                       x = "Detrital age (d)") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
                  scale_y_continuous(breaks = seq(-10, 40, by = 10),
                                     labels = scales::label_number(style_negative = "minus")) +
                  coord_cartesian(xlim = c(-1.8, 35), ylim = c(-10, 40), expand = FALSE, clip = "off") +
                  mytheme +
                  theme(axis.title.y = element_blank())


Fig_1 <- ( Fig_1a_top / Fig_1a_bottom ) / ( Fig_1b_left | Fig_1b_right ) +
         plot_layout(heights = c(0.3, 1, 1)) +
         plot_annotation(tag_levels = list(c("A", "", "B"))) &
         theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"))

ggsave(plot = Fig_1, filename = "Fig_1.pdf", device = cairo_pdf, path = "~/Desktop",
       width = 20, height = 20, units = "cm")


