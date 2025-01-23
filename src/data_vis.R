library(dplyr)
library(tidyr)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(ggplot2)

# load grid layers
vl <- st_read("data/Kaartjes/Vlaanderen.shp")
vl <- st_transform(vl, crs = 31370)
grid <- st_read("data/Kaartjes/hexgrid.shp")
grid <- st_intersection(grid, vl)

# define temporal range
yrange <- 2013:2019

# define mutant categories
mcat <- c("W", "M1", "M2", "M3")

# fit summary
stan_vars <- fit$metadata()$stan_variables
keep_vars <- c("bs0", "alpha_u", "rho_u", "alpha_b", "rho_b", 
               "rho_re", "rho_rs", "sigma_re", "sigma_rs")
fit_sum <- fit$summary(keep_vars[keep_vars %in% stan_vars])
fit_sum

# traceplots
mcmc_trace(fit$draws(keep_vars[keep_vars %in% stan_vars]))

# loo
loo <- fit$loo()
loo

# fit statistics
fit_stat <- fit %>% spread_draws(fit, fit_new)

# ideally the points are scattered around the diagonal line:
# this indicates a good model fit
ggplot(fit_stat) + 
  geom_point(aes(x = fit, y = fit_new, col = factor(.chain))) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.6) +
  labs(x = "Chi-squared discrepancy in data", 
       y = "Chi-squared discrepancy in replicated data") +
  guides(col = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 11))
ggsave("fig/fit_statistic.png", width = 12, height = 12, units = "cm", dpi = 300)

# proportions of mutants
tidy_props <- fit %>% spread_draws(prop[id, year, k]) # this may take a while...
props <- tidy_props %>% mean_qi()

props_grid <- props %>% 
  mutate(year = yrange[year], k = mcat[k]) %>%
  left_join(grid, by = "id") %>% st_as_sf() 

for (m in mcat) {
  ggplot(filter(props_grid, k == m)) + 
    geom_sf(data = vl, aes(), fill = NA, linewidth = 1.5) +
    geom_sf(aes(fill = prop), col = NA) + 
    facet_wrap(~year) +
    labs(fill = paste("Proportion of", m)) +
    scale_fill_gradient(low = "grey98", high = "firebrick") +
    theme_void()
  ggsave(paste0("fig/props_spatial_", m, ".png"), 
         width = 12, height = 16, units = "cm", dpi = 300)
}

# proportions of resistant mutants
props_res <- tidy_props %>%
  group_by(.draw, id, year) %>%
  summarize(
    W = prop[k == 1],  
    M = sum(prop[k %in% 2:4])
  ) %>%
  group_by(id, year) %>%
  mean_qi(W, M)

props_res_grid <- props_res %>% 
  mutate(year = yrange[year]) %>%
  left_join(grid, by = "id") %>% st_as_sf()

for (m in c("W", "M")) {
  ggplot(props_res_grid) + 
    geom_sf(data = vl, aes(), fill = NA, linewidth = 1.5) +
    geom_sf(aes(fill = !!sym(m)), col = NA) + 
    facet_wrap(~year) +
    labs(fill = paste("Proportion of", m)) +
    scale_fill_gradient(low = "grey98", high = "firebrick") +
    theme_void()
  ggsave(paste0("fig/props_res_spatial_", m, ".png"), 
         width = 12, height = 16, units = "cm", dpi = 300)
}

# overall yearly trend in proportions of mutants
yprops <- tidy_props %>%
  group_by(.draw, year, k) %>%
  summarize(yprop = mean(prop)) %>%
  group_by(year, k) %>%
  mean_qi() %>%
  mutate(year = yrange[year], k = mcat[k])

ggplot(yprops) + 
  geom_lineribbon(aes(x = year, y = yprop, ymin = .lower, ymax = .upper, 
                      col = k, fill = k), alpha = 0.2, linewidth = 1) +
  labs(col = "Proportion of", fill = "Proportion of") +
  facet_wrap(~k) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("fig/props_trend.png", width = 12, height = 12, units = "cm", dpi = 300)

# overall yearly trend in proportions of resitant mutants
yprops_res <- tidy_props %>%
  group_by(.draw, id, year) %>%
  summarize(
    W = prop[k == 1],  
    M = sum(prop[k %in% 2:4])
  ) %>%
  group_by(.draw, year) %>%
  summarize(W = mean(W), M = mean(M)) %>%
  group_by(year) %>%
  mean_qi() %>%
  mutate(year = yrange[year])

yprops_res <- yprops_res %>%
  rename(M.mean = M, W.mean = W) %>%
  pivot_longer(
    cols = starts_with("W") | starts_with("M"),
    names_to = c("k", ".value"),
    names_sep = "\\."
  ) %>%
  arrange(year, k)

ggplot(yprops_res) + 
  geom_lineribbon(aes(x = year, y = mean, ymin = lower, ymax = upper, 
                      col = k, fill = k), alpha = 0.2, linewidth = 1) +
  labs(col = "Proportion of", fill = "Proportion of") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("fig/props_res_trend.png", width = 12, height = 12, units = "cm", dpi = 300)

# Year-to-year percent changes in proportions of mutants
pct_changes <- tidy_props %>%
  group_by(.draw, id, k) %>%
  arrange(year) %>%
  mutate(
    pct_change = (prop - lag(prop)) / lag(prop) * 100
  ) %>%
  filter(!is.na(pct_change)) %>%
  group_by(year, id, k) %>%
  mean_qi()
  
pct_changes_grid <- pct_changes %>% 
  mutate(year = yrange[year], k = mcat[k]) %>%
  left_join(grid, by = "id") %>% st_as_sf() 

for (m in mcat) {
  m <- "M3"
  ggplot(filter(pct_changes_grid, k == m)) + 
    geom_sf(data = vl, aes(), fill = NA, linewidth = 1.5) +
    geom_sf(aes(fill = pct_change), col = NA) + 
    facet_wrap(~year) +
    labs(fill = paste("% change in proportions of", m)) +
    scale_fill_gradient2() +
    theme_void()
  ggsave(paste0("fig/pct_change_", m, ".png"), 
         width = 12, height = 16, units = "cm", dpi = 300)
}

# Year-to-year percent changes in proportions of resistant mutants
pct_changes_res <- tidy_props %>%
  group_by(.draw, id, year) %>%
  summarize(
    W = prop[k == 1],  
    M = sum(prop[k %in% 2:4])
  ) %>%
  group_by(.draw, id) %>%
  arrange(year) %>%
  mutate(
    pct_change.W = (W - lag(W)) / lag(W) * 100,
    pct_change.M = (M - lag(M)) / lag(M) * 100
  ) %>%
  filter(!is.na(pct_change.W), !is.na(pct_change.M)) %>%
  group_by(year, id) %>%
  mean_qi()

pct_changes_res_grid <- pct_changes_res %>% 
  mutate(year = yrange[year]) %>%
  left_join(grid, by = "id") %>% st_as_sf() 

for (m in c("W", "M")) {
  ggplot(pct_changes_res_grid) + 
    geom_sf(data = vl, aes(), fill = NA, linewidth = 1.5) +
    geom_sf(aes(fill = !!sym(m)), col = NA) + 
    facet_wrap(~year) +
    labs(fill = paste("% change in proportions of", m)) +
    scale_fill_gradient2() +
    theme_void()
  ggsave(paste0("fig/pct_change_res_", m, ".png"), 
         width = 12, height = 16, units = "cm", dpi = 300)
}

# Final note: the % changes can also be expressed as the actual change in proportions,
# which would be: (% change * previous' year proportion) / 100
# --> If you consider this easier to interpret, I will modify the script accordingly
