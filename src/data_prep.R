library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(sf)
library(ggplot2)
library(mapview)
library(BelgiumMaps.StatBel)
library(data.table)

data <- read_excel("data/2013 - 2019 BMK.xlsx")
data <- data %>%
  rename(X = X_lambert,
         Y = Y_lambert,
         location = Bekkennummer) %>%
  mutate(location = as.factor(location)) %>%
  rename(year = jaar) %>%
  mutate(
    iyear = .data$year - min(.data$year) + 1,
    cyear = .data$year - (min(.data$year) + (diff(range(.data$year))/2)),
    X = .data$X, Y = .data$Y,
    secondary = NA_real_,
    resistent = 1*(!(mutatie == 'WW')),
    MutatieM1 = 1*(str_detect(mutatie,regex('M1'))),
    MutatieM2 = 1*(str_detect(mutatie,regex('M2'))),
    MutatieM3 = 1*(str_detect(mutatie,regex('M3'))),
    mutation = case_when(
      MutatieM1 == 1 | (MutatieM1 == 1 & MutatieM2 == 1) | (MutatieM1 == 1 & MutatieM3 == 1) ~ "M1",
      MutatieM2 == 1 | (MutatieM2 == 1 & MutatieM3 == 1) ~ "M2",
      MutatieM3 == 1 ~ "M3",
      TRUE ~ "W"
    ),
    mutationF = factor(mutation, levels = c("W", "M1", "M2", "M3"), labels = 1:4)
  ) %>%
  select(.data$year, .data$Bekken, .data$X, .data$Y, .data$mutatie,
         .data$location, resistent, mutation, mutationF,
         iyear, cyear, X, Y, secondary)

data_sf <- st_as_sf(data, coords = c("X", "Y"), crs = 31370)
ggplot(data_sf) + geom_sf(aes(col = mutation))
mapview::mapview(data_sf, zcol = "mutation", legend = F)
ggplot(data_sf) + geom_sf(aes(col = mutation)) + facet_wrap(~year)

# Visualize grid layer used
gfiles <- unique(stringr::str_sub(list.files("data/Kaartjes/Hokken"), end = -5))
gfiles <- gfiles[gfiles != ""]
grid <- do.call("rbind", lapply(paste0("data/Kaartjes/Hokken/", gfiles ,".shp"), st_read))

ggplot(grid) + geom_sf(aes()) + 
  geom_sf(data = data_sf, aes(col = factor(location))) +
  theme(legend.position = "none")
mapview::mapview(grid) + mapview::mapview(data_sf, zcol = "location", legend = F)

# Generate hexagonal grid
vl <- st_read("data/Kaartjes/Vlaanderen.shp")
vl <- st_transform(vl, crs = 31370)
grid <- st_as_sf(st_make_grid(vl, cellsize = 5e3, square = F), crs = st_crs(vl))
ovl <- apply(st_intersects(grid, vl, sparse = F), 1, any)
grid <- mutate(grid[ovl,], id = row_number())
ggplot(vl) + geom_sf(aes()) + geom_sf(data = grid, fill = NA) +
  geom_sf(data = st_as_sf(data, coords = c("X", "Y"), crs = 31370), aes())
grid_ids <- apply(st_within(data_sf, grid, sparse = F), 1, which)
data$id <- grid_ids
st_write(grid, "data/Kaartjes/hexgrid.shp")

# Generate categorical counts per grid cell, year
griddata <- data %>% group_by(id, year, mutationF) %>% count()

# Temporal inforamtion
nyear <- length(unique(data$year))

# Spatial information
xy <- st_coordinates(data_sf)
xy_grid <- st_coordinates(st_centroid(grid))
# Standardise the spatial coordinates
standardise <- function(X, Y) {
  range_x <- diff(range(X))
  range_y <- diff(range(Y))
  x <- 2 * (X - min(X) - range_x/2) / (max(range_x, range_y))
  y <- 2 * (Y - min(Y) - range_y/2) / (max(range_x, range_y))
  return(cbind(X_std = x, Y_std = y))
}
xy_std <- standardise(xy[,1], xy[,2])
xy_grid_std <- standardise(xy_grid[,1], xy_grid[,2])
# Check for an optimal lon-lat ratio
16*max(xy_grid_std[,1])/6*max(xy_grid_std[,2])
# Define the approximate GP design matrix
xy_design <- as.matrix(expand.grid(1:16,1:6))

# Counts and proportions
count_data <- expand.grid(id = 1:nrow(grid), year = unique(data$year), mutationF = 1:4) %>%
  left_join(mutate(griddata, mutationF = as.numeric(mutationF)), by = c("id", "year", "mutationF")) %>%
  pivot_wider(names_from = "mutationF", values_from = "n")

counts <- count_data %>%
  select(-c(id, year)) %>%
  data.matrix()
counts[is.na(counts)] <- 0

prop_data <- count_data
colnames(prop_data)[3:6] <- paste0("K", 1:4)
prop_data[,3:6] <- (prop_data[,3:6]+1e-6) /
  apply(prop_data[,3:6]+1e-6, 1, sum, na.rm = TRUE)
allNA <- apply(is.na(prop_data[,3:6]), 1, all)
prop_data[!allNA,] <-
  prop_data[!allNA,] %>%
  replace_na(list(K1 = 0, K2 = 0, K3 = 0, K4 = 0))
E <- filter(prop_data, year == min(year)) %>%
  summarise(across(K1:K4, \(x) mean(x, na.rm = TRUE)))

stan_input <- list(
  N = nrow(counts),
  R = nrow(grid),
  T = length(unique(griddata$year)),
  K = length(levels(griddata$mutationF)),
  mu = as.numeric(E),
  y = counts,
  index = matrix(1:nrow(counts), nrow = nrow(grid), ncol = nyear),
  Rs = length(unique(griddata$id)),
  rs = unique(griddata$id),
  rns = (1:nrow(grid))[-unique(griddata$id)],
  year = (1:nyear) - (nyear/2),
  n_bf = nrow(xy_design),       
  xy_design = xy_design,       
  xy_std = xy_grid_std,         
  xy_lim = c(apply(xy_grid_std, 2, max)*6/4),
  scaling_factor = sqrt(prod(apply(xy_grid_std, 2, var)))
)
saveRDS(stan_input, "stan/input.rds")
