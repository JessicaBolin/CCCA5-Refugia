# Figure 14: Change in persistence of refugia in the Channel Islands relative to the baseline period (1990–2010)
# Author: Jessica Bolin

# Dependencies ------------------------------------------------------------

source("scripts/_helpers.R")
ca <- st_read("shapefiles/processed/usa_contiguous.shp")
eez <- st_read("shapefiles/processed/westcoast_eez.shp")
bathy2 <- sf::st_read("shapefiles/processed/cont_shelf_50m.shp")
test <- sf::st_read("shapefiles/processed/ca_stateboundary_tigris.shp")

# Function ----------------------------------------------------------------

suffix = "far_20702099-19902009"
suffix = "near_20202049-19902009"
thresh_perc = 50
model = "ens"
area = "channel_islands"
what = "_all_defs"

if (model == "^zoom") { model2 <- "zoom" } else { model2 <- model }

r <- rast(paste0(refugia_out_wd , "/_ens/", what, "/_3_delta_rasts/", model2, 
                 "_delta_rast_", thresh_perc, suffix, ".nc"))

cropped <- r
cropped <- crop(cropped, channel_islands)
e <- ext(cropped)

# Expand by 10% on all sides
yrange <- ymax(e) - ymin(e)

ext_zoomed <- st_bbox(c(
  xmin = xmin(e),
  xmax = xmax(e),
  ymin = ymin(e),
  ymax = ymax(e) + 0.2 * yrange
), crs = st_crs(cropped))  # attach CRS from raster

vals <- c(global(cropped, "min", na.rm = TRUE)[[1]],
          global(cropped, "max", na.rm = TRUE)[[1]]) 
max_abs <- max(abs(vals))  # symmetric around zero

ttt <- colorspace::hcl_palettes(type = "diverging")

# Pick a named diverging palette from the list
my_palette <- colorspace::diverging_hcl(n = 7, palette = "Blue-Red 2")

t1 <- tm_shape(cropped, bbox = ext_zoomed) +
  tmap::tm_raster(col.scale = tmap::tm_scale_categorical(),
                  col = "white", 
                  col.legend = tmap::tm_legend_hide()) +      
  tm_graticules(ticks = T, lwd = 0.5, col = "grey50", labels.size = 0.5) +
  tmap::tm_shape(ca) +
  tmap::tm_polygons() +
  tm_shape(cropped) +
  tm_raster(col.scale = tm_scale_continuous(values = rev(my_palette),
                                            midpoint = 0,
                                            limits = c(-30, 30))) +
  tm_shape(ca) +
  tm_polygons(fill_alpha = 1) +
  tm_scalebar(position = c("left", "bottom"), text.size = 0.7) +
  tmap::tm_shape(bathy2[0]) +
  tmap::tm_lines(lty= "88") +
  tm_shape(test[0]) +
  tm_lines() 

t1

tmap_save(t1, 
          paste0("final_report_CFCCA/figures/red_abalone/out/delta_", 
                 suffix, thresh_perc, ".png"))

