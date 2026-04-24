# Figures 9-13: Projected persistence of climate refugia for ______
# Jessica Bolin

# Dependencies ------------------------------------------------------------

source("scripts/_helpers.R")
tmap_options(show.messages = F)
ca <- st_read("shapefiles/processed/usa_contiguous.shp")
eez <- st_read("shapefiles/processed/westcoast_eez.shp")
bathy2 <- sf::st_read("shapefiles/processed/cont_shelf_50m.shp")
test <- sf::st_read("shapefiles/processed/ca_stateboundary_tigris.shp")

# Change this to area/inputs of interest -----------------------------------

area = "sonoma"
what = "alldefs"
thresh = "lib"
yr_range = 1990:2010

# FUNCTION ----------------------------------------------------------------

listy <- list.files(paste0(refugia_out_wd, "/_ens/_all_defs_mp/_2_persistence_rasts/"),
                    full.names = T)
listyy <- listy[grep(yr_range %>% min, listy)]
listyyy <- listyy[grep(thresh, listyy)]
listyyyy <- listyyy[grep(paste0(area, "_"), listyyy)]
listyyyy <- listyyyy[grep(paste0(what, "_"), listyyyy)]
r <- rast(listyyyy)

pos = c("right", "top")
breaks <- c(0, 20, 40, 60, 80, 100)

t1 <- tm_shape(r) + #, bbox = ext_zoomed) +
  tmap::tm_raster(col.scale = tmap::tm_scale_categorical(),
                  col = "white", 
                  col.legend = tmap::tm_legend_hide()) +     
  tm_graticules(ticks = T, lwd = 0.5, col = "grey50", labels.size = 1) +
  tmap::tm_shape(ca) +
  tmap::tm_polygons() +
  tm_shape(r) +
  tmap::tm_raster(col.scale = 
                    tmap::tm_scale_intervals(breaks = breaks,
                                             values = viridis::magma(length(breaks)-1)), 
                  col.legend = tmap::tm_legend(title = "Persistence (%)",
                                               text.size = 1,
                                               frame = T,
                                               frame.lwd = 0.001,
                                               title.size = 1.1,
                                               bg.alpha = 0,
                                               width = 5,
                                               height = 8)) +
  tm_shape(ca) +
  tm_polygons(fill_alpha = 0.9) +
  tm_scalebar(position = c("left", "bottom"), text.size = 0.5) +
  tmap::tm_shape(bathy2[0]) +
  tmap::tm_lines(lty = "88") +
  tm_shape(test[0]) +
  tm_lines() +
  tm_shape(ca) +
  tm_polygons(fill_alpha = 0.9) 

t1

tmap_save(t1, 
          paste0("final_report_CFCCA/figures/red_abalone/out/", 
                 what, "_", area, "_", thresh, "_",  yr_range %>% min, "-", 
                 yr_range %>% max, "_.png"))


