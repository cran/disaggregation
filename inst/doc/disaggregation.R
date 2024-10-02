## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  fig.width = 7,
  eval = TRUE
)

## -----------------------------------------------------------------------------
library(SpatialEpi, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(disaggregation, quietly = TRUE)
library(ggplot2)
library(sf)
library(terra)

polygons <- sf::st_as_sf(NYleukemia$spatial.polygon)

df <- cbind(polygons, NYleukemia$data)


ggplot() + geom_sf(data = df, aes(fill = cases / population)) + scale_fill_viridis_c(lim = c(0, 0.003))


## ----fig.show='hold'----------------------------------------------------------

bbox <- sf::st_bbox(df)

extent_in_km <- 111*(bbox[c(3, 4)] - bbox[c(1, 2)])
n_pixels_x <- floor(extent_in_km[[1]])
n_pixels_y <- floor(extent_in_km[[2]])

r <- terra::rast(ncols = n_pixels_x, nrows = n_pixels_y)
terra::ext(r) <- terra::ext(df)

data_generate <- function(x){
  rnorm(1, ifelse(x %% n_pixels_x != 0, x %% n_pixels_x, n_pixels_x), 3)
}

terra::values(r) <- sapply(seq(terra::ncell(r)), data_generate)
r2 <- terra::rast(ncol = n_pixels_x, nrow = n_pixels_y)
terra::ext(r2) <- terra::ext(df)
terra::values(r2) <- sapply(seq(terra::ncell(r2)), 
                     function(x) rnorm(1, ceiling(x/n_pixels_y), 3))


cov_stack <- terra::rast(list(r, r2))
cov_stack <- terra::scale(cov_stack)
names(cov_stack) <- c('layer1', 'layer2')


## ----fig.show='hold'----------------------------------------------------------
extracted <- terra::extract(r, terra::vect(df$geometry), fun = sum)
n_cells <- terra::extract(r, terra::vect(df$geometry), fun = length)
df$pop_per_cell <- df$population/n_cells$lyr.1
pop_raster <- terra::rasterize(terra::vect(df), cov_stack, field = 'pop_per_cell')


## ----fig.show='hold'----------------------------------------------------------
df <- sf::st_buffer(df, dist = 0)

## ----fig.show='hold'----------------------------------------------------------
data_for_model <- prepare_data(polygon_shapefile = df,
                               covariate_rasters = cov_stack,
                               aggregation_raster = pop_raster,
                               response_var = 'cases',
                               id_var = 'censustract.FIPS',
                               mesh_args = list(cutoff = 0.01,
                                                offset = c(0.1, 0.5),
                                                max.edge = c(0.1, 0.2),
                                                resolution = 250),
                               na_action = TRUE)

## ----fig.show='hold'----------------------------------------------------------
plot(data_for_model)

## ----fig.show='hold'----------------------------------------------------------
model_result <- disag_model(data_for_model,
                            iterations = 1000,
                            family = 'poisson',
                            link = 'log',
                            priors = list(priormean_intercept = 0,
                                          priorsd_intercept = 2,
                                          priormean_slope = 0.0,
                                          priorsd_slope = 0.4,
                                          prior_rho_min = 3,
                                          prior_rho_prob = 0.01,
                                          prior_sigma_max = 1,
                                          prior_sigma_prob = 0.01,
                                          prior_iideffect_sd_max = 0.05,
                                          prior_iideffect_sd_prob = 0.01))

## ----fig.show='hold'----------------------------------------------------------
plot(model_result)

## ----fig.show='hold'----------------------------------------------------------
preds <- predict(model_result, 
                 N = 100, 
                 CI = 0.95)

plot(preds)

