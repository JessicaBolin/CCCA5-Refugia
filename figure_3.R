# Figure 3: projected habitat suitability for sardine
# Authors: Nerea Lezama Ochoa and Barb Muhling

# Dependencies ------------------------------------------------------------

# Clean environment
rm(list = ls())

# Load required packages
library(ncdf4)
library(ggplot2)
library(viridis)
library(reshape2)
library(gridExtra)
library(maps)

# -------------------------------
# Directory setup and file lists
# -------------------------------
nc_dir <- "/Volumes/Triple_Bottom_Line/Nerea_working/Climate_resilence/sardine/output/"
nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)

# Split by time period
nc_files_hist <- nc_files[grepl("199[0-9]|200[0-9]|2010", nc_files)] # historical
nc_files_fut <- nc_files[grepl("207[0-9]|208[0-9]|209[0-9]|2100", nc_files)] # long-term future

# Variables
variables_noST <- c("pred_sdmTMB_GFDL_noST", "pred_sdmTMB_HADL_noST", "pred_sdmTMB_IPSL_noST")
variables_withST <- c("pred_sdmTMB_GFDL_ST", "pred_sdmTMB_HADL_ST", "pred_sdmTMB_IPSL_ST")

# -------------------------------
# Get grid coordinates
# -------------------------------
nc <- nc_open(nc_files_hist[1])
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

# -------------------------------
# Helper functions
# -------------------------------
extract_predictions <- function(files, varname) {
  preds <- list()
  for (file in files) {
    nc <- nc_open(file)
    if (varname %in% names(nc$var)) {
      data <- ncvar_get(nc, varname, start = c(1, 1, 100), count = c(-1, -1, 1))
      preds[[length(preds) + 1]] <- data
    }
    nc_close(nc)
  }
  return(preds)
}

compute_stats <- function(files, variables) {
  all_data <- list()
  for (v in variables) {
    preds <- extract_predictions(files, v)
    if (length(preds) > 0) {
      all_data <- c(all_data, preds)
    }
  }
  arr <- simplify2array(all_data)
  mean_vals <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  sd_vals <- apply(arr, c(1, 2), sd, na.rm = TRUE)
  list(mean = mean_vals, sd = sd_vals)
}

# -------------------------------
# Compute stats for noST and withST
# -------------------------------
stats_hist_noST <- compute_stats(nc_files_hist, variables_noST)
stats_fut_noST <- compute_stats(nc_files_fut, variables_noST)
delta_noST <- stats_fut_noST$mean - stats_hist_noST$mean

stats_hist_withST <- compute_stats(nc_files_hist, variables_withST)
stats_fut_withST <- compute_stats(nc_files_fut, variables_withST)
delta_withST <- stats_fut_withST$mean - stats_hist_withST$mean

# -------------------------------
# Prepare data frames
# -------------------------------
lon_lat_grid <- expand.grid(lon = lon, lat = lat)

to_df <- function(mean_vals, sd_vals, delta_vals) {
  list(
    df_mean_hist = data.frame(lon_lat_grid, pred = as.vector(mean_vals$mean)),
    df_mean_fut = data.frame(lon_lat_grid, pred = as.vector(sd_vals$mean)),
    df_sd_hist = data.frame(lon_lat_grid, sd = as.vector(mean_vals$sd)),
    df_sd_fut = data.frame(lon_lat_grid, sd = as.vector(sd_vals$sd)),
    df_delta = data.frame(lon_lat_grid, diff = as.vector(delta_vals))
  )
}

dfs_noST <- to_df(stats_hist_noST, stats_fut_noST, delta_noST)
dfs_withST <- to_df(stats_hist_withST, stats_fut_withST, delta_withST)

# -------------------------------
# Map style layers
# -------------------------------
usa_states <- map_data("state")

# Grey landmass
land_layer <- geom_polygon(
  data = usa_states,
  aes(x = long, y = lat, group = group),
  fill = "grey80", color = NA
)

# State outlines (U.S.)
outline_layer <- geom_path(
  data = usa_states,
  aes(x = long, y = lat, group = group),
  color = "black",
  linewidth = 0.3
)

# Emphasize California
california_outline <- geom_path(
  data = subset(usa_states, region == "california"),
  aes(x = long, y = lat, group = group),
  color = "black",
  linewidth = 0.5
)

# General map theme
coord_limits <- coord_quickmap(xlim = c(-127, -115), ylim = c(33, 43))
theme_map <- theme(
  plot.title = element_text(size = 10),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8),
  panel.background = element_rect(fill = "white"),
  panel.grid = element_blank()
)

# -------------------------------
# Plotting function
# -------------------------------
plot_maps <- function(df_mean_hist, df_mean_fut, df_sd_hist, df_sd_fut, df_delta, label) {
  
  # A helper function for consistency across all maps
  base_theme <- list(
    outline_layer,
    california_outline,
    coord_limits,
    theme_map
  )
  
  # Define the same land layer for all (use a single consistent color)
  land_layer_top <- geom_polygon(
    data = usa_states,
    aes(x = long, y = lat, group = group),
    fill = "grey80", color = NA
  )
  
  p1 <- ggplot() +
    geom_tile(data = df_mean_hist, aes(x = lon, y = lat, fill = pred)) +
    scale_fill_gradientn(colors = viridis::magma(255, direction = 1),
                         limits = c(0, 0.7), na.value = NA, oob = scales::squish) +
    land_layer_top +
    base_theme +
    labs(title = paste("Mean Prediction (1990–2010)", label),
         x = "Longitude", y = "Latitude")
  
  p2 <- ggplot() +
    geom_tile(data = df_mean_fut, aes(x = lon, y = lat, fill = pred)) +
    scale_fill_gradientn(colors = viridis::magma(255, direction = 1),
                         limits = c(0, 0.7), na.value = NA, oob = scales::squish) +
    land_layer_top +
    base_theme +
    labs(title = paste("Mean Prediction (2070–2100)", label),
         x = "Longitude", y = "Latitude")
  
  p3 <- ggplot() +
    geom_tile(data = df_sd_hist, aes(x = lon, y = lat, fill = sd)) +
    scale_fill_viridis(limits = c(0, 0.1), na.value = NA, oob = scales::squish) +
    land_layer_top +
    base_theme +
    labs(title = paste("SD Across Models (1990–2010)", label),
         x = "Longitude", y = "Latitude")
  
  p4 <- ggplot() +
    geom_tile(data = df_sd_fut, aes(x = lon, y = lat, fill = sd)) +
    scale_fill_viridis(limits = c(0, 0.1), na.value = NA, oob = scales::squish) +
    land_layer_top +
    base_theme +
    labs(title = paste("SD Across Models (2070–2100)", label),
         x = "Longitude", y = "Latitude")
  
  p5 <- ggplot() +
    geom_tile(data = df_delta, aes(x = lon, y = lat, fill = diff)) +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0,
                         limits = c(-0.4, 0.4), oob = scales::squish) +
    land_layer_top +   # Ensures same land color
    base_theme +
    labs(title = paste("Difference (Future - Historical)", label),
         x = "Longitude", y = "Latitude")
  
  list(p1, p2, p3, p4, p5)
}


# -------------------------------
# Generate and save combined plots
# -------------------------------
row_noST <- plot_maps(dfs_noST$df_mean_hist, dfs_noST$df_mean_fut,
                      dfs_noST$df_sd_hist, dfs_noST$df_sd_fut, dfs_noST$df_delta, "(noST)")
row_withST <- plot_maps(dfs_withST$df_mean_hist, dfs_withST$df_mean_fut,
                        dfs_withST$df_sd_hist, dfs_withST$df_sd_fut, dfs_withST$df_delta, "(withST)")

combined_plot <- grid.arrange(grobs = c(row_noST, row_withST), ncol = 5)

setwd("/Volumes/Triple_Bottom_Line/Nerea_working/Climate_resilence/results/new_results")
ggsave("summary_sdm_maps_rows_reportNere.tif", combined_plot, width = 18, height = 8, dpi = 300)
