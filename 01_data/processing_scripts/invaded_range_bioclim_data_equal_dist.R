
##
## Prediction grid - equal spaced by distance
##

library(dplyr)
library(sf)
library(raster)
library(readr)
library(ggplot2)
library(viridis)

# Load in invaded-range observations
genotype_locations <-
  read.csv("data/BioclimateOfOrigin_AllGenotypes.csv", header=T) %>%
  dplyr::select("source" = site_code, genotype, "longitude" = lon, "latitude" = lat)

# Join with key matrix for connecting with genotype matrix
genotype_codes <-
  read.csv("data/BRTEcg_genotypesCode.csv", header = T) %>%
  left_join(genotype_locations)  %>%
  filter(!is.na(longitude)) %>%
  arrange(SNPmatrix_column) %>%
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>%
  ### Keep lat, long as columns too
  dplyr::mutate(
    longitude = sf::st_coordinates(.)[, 1],
    latitude = sf::st_coordinates(.)[, 2]
  )

# Extract invaded range observation coords
invaded_samps <- genotype_codes %>% dplyr::select(x = longitude, y = latitude) %>% as.data.frame() #%>% st_as_sf(coords = c('x', 'y'), crs = 4326)

# Convert points to sf and transform to projected CRS EPSG:5070 (meters - CRS commonly used for continental US)
invaded_sf <- st_as_sf(invaded_samps, coords = c("x", "y"), crs = 4326)
invaded_proj <- st_transform(invaded_sf, crs = 5070)  # US Equal Area (meters)

# Set dimensions of grid here:
num_rows <- 100
num_cols <- 100

grid_centers <- st_make_grid(
  invaded_proj,
  n = c(num_cols, num_rows),     # number of cols (x), rows (y)
  what = "centers"               # return cell centroids
)
grid_centers_sf <- st_sf(geometry = grid_centers)
st_crs(grid_centers_sf) <- st_crs(invaded_proj)

# Transform grid centers back to WGS84 for raster extraction
grid_centers_wgs <- st_transform(grid_centers_sf, crs = 4326)

# Bioclimatic covariates were gathered from the CHELSA data repository, which can be found here: https://chelsa-climate.org/bioclim/
# Full variable names can be found here: https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf

# Extract variable names for raster stack
cov_names <- colnames(invaded)[6:24]

rlist <- lapply(1:19, function(i) {
  path <- paste0("chelsav2/GLOBAL/climatologies/2011-2040/GFDL-ESM4/ssp126/bio/CHELSA_bio",
                 i, "_2011-2040_gfdl-esm4_ssp126_V.2.1.tif")
  r <- raster(path)
  names(r) <- cov_names[i]
  return(r)
})
stacked <- stack(rlist)

# Extract covariate values at each grid center
grid_coords <- st_coordinates(grid_centers_wgs)
grid_vals <- raster::extract(stacked, grid_coords, df = TRUE)

# Combine coordinates and values
grid_df <- cbind(as.data.frame(grid_coords), grid_vals) %>%
  rename(lon = X, lat = Y) %>%
  dplyr::select(lon, lat, all_of(cov_names))