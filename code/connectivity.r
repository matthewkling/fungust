
library(tidyverse)
library(raster)
library(windscape)
library(gdistance)


## functions #############

direction <- "downwind"


unwrap <- function(x, width=20){
  x <- x %>% crop(extent(180-width, 180, -90, 90)) %>% shift(-360) %>% merge(x)
  x <- x %>% crop(extent(-180, -180+width, -90, 90)) %>% shift(360) %>% merge(x)
  x
}


flatten <- function(x, i, varname="dist"){
  # x = distance matrix
  # i = vector of site numbers
  x <- as.data.frame(x)
  names(x) <- paste0("p", i)
  x$from <- names(x)
  x <- gather(x, to, dist, -from)
  x$from <- as.integer(sub("p", "", x$from))
  x$to <- as.integer(sub("p", "", x$to))
  names(x)[3] <- varname
  x
}

wind_hours <- function(sites, infile, water = 1, disagg = NULL){
  
  # create raster of weights for land/water
  # weights <- raster("f:/cfsr/land.tif") %>%  
  #   rotate() %>% 
  #   unwrap(180) %>%
  #   reclassify(c(NA, NA, water))
  
  # load windrose data and apply weights
  rose <- stack(infile) %>% 
    rotate() %>% 
    unwrap(width = 180)# %>%
    # "*"(weights) 
  
  if(!is.null(disagg)) rose <- rose %>% 
    disaggregate(disagg, "bilinear")
  
  # make transition graph
  wind <- rose %>% 
    add_coords() %>% 
    transition_stack(windflow, directions = 8, 
                     symm = F, direction = direction)
  
  # calculate wind distances, convert to hours
  costDistance(wind, sites) / 3600
}



## analysis ###################

# prep site locations
sites <- read_csv("df_sample_Data_coordinates.csv") %>%
  mutate(index = 1:n())
coordinates(sites) <- c("Longitude", "Latitude")
crs(sites) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
sites <- spTransform(sites, crs(raster(#"f:/cfsr/land.tif",
  "e:/wind/winds_of_change/tailwinds/data/windrose/windrose_p1_wnd10m_SON.tif")))

w10m_annual <- sites %>%
  wind_hours("e:/wind/winds_of_change/tailwinds/data/windrose/windrose_p1_wnd10m.tif")
# w10m_annual <- sites %>%
#   wind_hours("e:/wind/winds_of_change/tailwinds/data/windrose/windrose_p1_wnd10m.tif",
#              disagg = 2) %>% "/"(2)
w10m_annual %>% as.data.frame() %>% write_csv("output/wind_matrix_annual.csv")

w10m_autumn <- sites %>% 
  wind_hours("e:/wind/winds_of_change/tailwinds/data/windrose/windrose_p1_wnd10m_SON.tif")
# w10m_autumn <- sites %>% 
#   wind_hours("e:/wind/winds_of_change/tailwinds/data/windrose/windrose_p1_wnd10m_SON.tif",
#              disagg = 2) %>% "/"(2)
w10m_autumn %>% as.data.frame() %>% write_csv("output/wind_matrix_SON.csv")

dst <- geosphere::distm(sites, sites) / 1000


d <- flatten(w10m_annual, sites$index, "wind_hours_10m_annual") %>%
  left_join(flatten(w10m_autumn, sites$index, "wind_hours_10m_SON")) %>%
  left_join(flatten(dst, sites$index, "distance_km"))
d %>% write_csv("output/wind_connectivity.csv")


sites$cell_id <- raster("e:/wind/winds_of_change/tailwinds/data/windrose/windrose_p1_wnd10m.tif") %>%
  rotate() %>%
  setValues(1:ncell(.)) %>%
  raster::extract(sites)
sites <- as.data.frame(sites)
sites %>% write_csv("output/sites.csv")
