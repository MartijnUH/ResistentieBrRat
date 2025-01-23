library(cmdstanr); library(parallel)

input <- readRDS("stan/input.rds")
mod_type <- "multinom2d" # else choose: "multinom"
mod <- cmdstan_model(stan_file = paste0("stan/", mod_type, ".stan"))
ncores <- detectCores()-1

init_values <- list(
  bs0 = c(0.05, 0.4, -0.3, 0.2),
  alpha_u = c(0.9, 1.8, 1.8, 3.5),
  rho_u = c(0.8, 0.2, 0.3, 0.9),
  alpha_b = c(0.8, 1.9),
  rho_b = c(0.9, 0.9)
)

fit <- mod$sample(data = input, seed = 2025, 
                  chains = 1, parallel_chains = 1,
                  iter_warmup = 50, iter_sampling = 50, refresh = 1, 
                  adapt_delta = 0.90, init = function() init_values) 
fit$save_output_files(dir = "output", basename = mod_type)