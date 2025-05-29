library(dplyr)
library(sf)
library(raster)
library(readr)
library(ggplot2)
library(viridis)

##
## Native Range prediction grid 
##

# Load and filter original data
all_brte_samples <- read_csv(here('01_data', 'data', 'bioclim_and_other','BRTEclim.csv'), show_col_types = FALSE)
native <- all_brte_samples %>%
  filter(range == 'native') %>%
  left_join(read_csv(here('01_data', 'data', 'bioclim_and_other','307tips.csv'), show_col_types = FALSE), by = 'PopNum') %>%
  filter(!is.na(Longitude)) %>%
  arrange(IBS.id) %>%
  mutate(snp_map_col = IBS.id + 3) # + 3 to account for the first three non-genotype columns

# Take a smaller spatial chunk 
# Set up lat and long bounds here:
min_long <- 0
max_long <- 22
min_lat <- 45
max_lat <- 55

native_samps <- native %>%
  dplyr::select(NewSiteCode, PopNum, snp_map_col, x = Longitude, y = Latitude) %>%
  as.data.frame() %>% filter(y > min_lat & y < max_lat & x > min_long & x < max_long)

# Convert points to sf and transform coordinate reference system (EPSG:3035 is good for Europe)
native_sf <- st_as_sf(native_samps, coords = c("x", "y"), crs = 4326)
native_proj <- st_transform(native_sf, crs = 3035)  

# Set dimensions of grid here:
num_rows <- 100
num_cols <- 100

# Construct grid
grid_centers <- st_make_grid(
  native_proj,
  n = c(num_cols, num_rows),     # number of cols (x), rows (y)
  what = "centers"               # return cell centroids
)
grid_centers_sf <- st_sf(geometry = grid_centers)
st_crs(grid_centers_sf) <- st_crs(native_proj)

# Convert points to sf and transform coordinate reference system (EPSG:3035 is good for Europe)
native_sf <- st_as_sf(native_samps, coords = c("x", "y"), crs = 4326)
native_proj <- st_transform(native_sf, crs = 3035)  

# Set dimensions of grid here:
num_rows <- 100
num_cols <- 100

# Construct grid
grid_centers <- st_make_grid(
  native_proj,
  n = c(num_cols, num_rows),     # number of cols (x), rows (y)
  what = "centers"               # return cell centroids
)
grid_centers_sf <- st_sf(geometry = grid_centers)
st_crs(grid_centers_sf) <- st_crs(native_proj)

# Transform grid centers back to WGS84 for raster extraction
grid_centers_wgs <- st_transform(grid_centers_sf, crs = 4326)

# Load bioclimatic rasters
cov.names <- read_csv(here('01_data', 'data', 'bioclim_and_other', 'BioclimateOfOrigin_AllGenotypes.csv'), show_col_types = FALSE) %>%
  colnames() %>% .[6:24]

rlist <- lapply(1:19, function(i) {
  path <- paste0(here('01_data', 'data', 'bioclim_and_other','chelsav2/GLOBAL/climatologies/2011-2040/GFDL-ESM4/ssp126/bio/CHELSA_bio'),
                 i, "_2011-2040_gfdl-esm4_ssp126_V.2.1.tif")
  r <- raster(path)
  names(r) <- cov.names[i]
  return(r)
})
stacked <- stack(rlist)

# Extract covariate values at each grid center
grid_coords <- st_coordinates(grid_centers_wgs)
grid_vals <- raster::extract(stacked, grid_coords, df = TRUE)

# Combine coordinates and bioclim values
grid_df <- cbind(as.data.frame(grid_coords), grid_vals) %>%
  rename(lon = X, lat = Y) %>%
  dplyr::select(lon, lat, all_of(cov.names))

# Check
# ggplot() +
#   geom_point(data = grid_df, aes(x = lon, y = lat, color = prc.wet.q)) +
#   scale_color_viridis() +
#   coord_sf(crs = 4326) +
#   theme_minimal()
# 
# ggplot() +
#   geom_sf(data = grid_centers_sf, size = 0.2) +
#   coord_sf(crs = 3035)

# Save
save_path <- here('01_data', 'data', 'bioclim_and_other', 'native_range_bioclim')
write_csv(grid_df, file = save_path)