
library(tidyverse)
library(windscape)


build_rose <- function(rose, xy, fact = 1, expand = 1.5){
      rose %>%
            rast() %>% # load file from disk
            rotate() %>% # convert to standard coordinate range
            tesselate(30) %>%
            crop(extent(xy) * expand) %>% # limit domain for computational efficiency
            as("wind_rose") %>% # declare as wind rose object
            downscale(fact) # increase resolution
}

flatten <- function(x, i, varname = "dist"){
      # x = distance matrix
      # i = vector of site numbers
      x <- as.data.frame(x)
      names(x) <- paste0("p", i)
      x$from <- names(x)
      x <- gather(x, to, dist, -from)
      x$from <- as.integer(sub("p", "", x$from))
      x$to <- as.integer(sub("p", "", x$to))
      names(x)[3] <- varname
      as_tibble(x)
}



# load sampling sites
sites <- read_csv("data/df_sample_Data_coordinates.csv") %>%
      mutate(index = 1:n())
ll <- sites %>%
      dplyr::select(Longitude, Latitude) %>%
      as.matrix()

# compute wind connectivity
w <- "~/data/CFSR/multidecadal_roses/wnd10m_annual_p1.tif" %>%
      build_rose(ll, fact = 5) %>%
      picd(ll)

# export distance matrices
w$wind_dist %>% as.data.frame() %>% write_csv("output/wind_matrix.csv")
w$point_dist %>% as.data.frame() %>% write_csv("output/point_dist_matrix.csv")
w$cell_dist %>% as.data.frame() %>% write_csv("output/cell_dist_matrix.csv")

# export combined results in long format
d <- flatten(w$wind_dist, sites$index, "wind_hours") %>%
      left_join(flatten(w$point_dist, sites$index, "point_distance_km")) %>%
      left_join(flatten(w$cell_dist, sites$index, "cell_distance_km"))
d %>% write_csv("output/wind_connectivity.csv")

# export sites
sites %>% write_csv("output/sites.csv")
