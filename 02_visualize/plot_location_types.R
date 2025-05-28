library(maps)
library(tigris)
library(USA.state.boundaries)

#setwd('bromecast-data/')

# Load in state/country borders for plotting
data("state_boundaries_wgs84")
north_america <- map_data("world", region = c("USA", "Canada", "Mexico"))

### 1. Read in reference genotype locations ------------------------------------
ref_genos <- read.csv('gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv') %>%
  mutate(type = 'Reference Genotype Location') %>%
  st_as_sf(., coords = c("lon", "lat"), crs = 4326) %>%
  dplyr::mutate(
    longitude = sf::st_coordinates(.)[, 1],
    latitude = sf::st_coordinates(.)[, 2]
  )

# length(unique(ref_genos$genotype))
# length(unique(ref_genos$site_code))
# idx <- which(duplicated(ref_genos$site_code))
# ref_genos[idx, ] # same site code but different genotype number 


## 2. Read in satellite info for each year and join
df1 <- read.csv("https://github.com/pbadler/bromecast-data/raw/main/satellites/rawdata/SiteInfo_2020-2021.csv")
df1 <- df1 %>% dplyr::select(site = 'Site.code', lon = 'Longitude..decimal.degrees.', lat = 'Latitude..decimal.degrees.')

df2 <- read.csv("https://github.com/pbadler/bromecast-data/raw/main/satellites/rawdata/SiteInfo_2021-2022.csv")
df2 <- df2 %>% dplyr::select(site = 'Site.code', lon = 'Longitude..decimal.degrees.', lat = 'Latitude..decimal.degrees.')

df3 <- read.csv("https://github.com/pbadler/bromecast-data/raw/main/satellites/rawdata/SiteInfo_2022-2023.csv")
df3 <- df3 %>% dplyr::select(site = 'Site.code', lon = 'Longitude..decimal.degrees.', lat = 'Latitude..decimal.degrees.')

df4 <- read.csv("https://github.com/pbadler/bromecast-data/raw/main/satellites/rawdata/SiteInfo_2023-2024.csv")
df4 <- df4 %>% dplyr::select(site = 'Site.code', lon = 'Longitude..decimal.degrees.', lat = 'Latitude..decimal.degrees.')

prelim_sat_sites <- rbind(df1, df2, df3, df4)  

bromecast_sites <- read.csv(here('01_data', 'data', 'bioclim_and_other', 'BromecastSites.csv')) %>%
  dplyr::select(site = 'Site.code', lon = 'Longitude..decimal.degrees.', lat = 'Latitude..decimal.degrees.')

anti_join(prelim_sat_sites, bromecast_sites, by = 'site')
anti_join(bromecast_sites, prelim_sat_sites, by = 'site')

joined_sat_sites <- rbind(prelim_sat_sites, bromecast_sites) %>% 
  distinct(site) 

### 2. Load in satellite site locations ----------------------------------------
sat_sites <- read.csv('satellites/rawdata/BromecastSites.csv') %>% 
  rename(lon = Longitude..decimal.degrees., lat = Latitude..decimal.degrees.) %>%
  st_as_sf(., coords = c("lon", "lat"), crs = 4326) %>%
  dplyr::mutate(
    longitude = sf::st_coordinates(.)[, 1],
    latitude = sf::st_coordinates(.)[, 2]
  ) %>% 
  mutate(type = 'Satellite Site') %>% 
  dplyr::select(Site.code, latitude, longitude, type)

### 3. Read in common garden locations -----------------------------------------
comm_gardens <- read.csv('gardens/rawdata/garden_info.csv') %>%
  mutate(type = 'Common Garden') %>%
  st_as_sf(.,
           coords = c("garden_lon", "garden_lat"),
           crs = 4326) %>%
  dplyr::mutate(
    longitude = sf::st_coordinates(.)[, 1],
    latitude = sf::st_coordinates(.)[, 2]
  )


################################################################################
## JOINT PLOT  -----------------------------------------------------------------
################################################################################

p <- 
  ggplot() +
  geom_polygon(data = north_america, aes(x = long, y = lat, group = group), fill = 'white', color = 'black') +
  geom_sf(data = state_boundaries_wgs84, fill = NA, color = 'black', linewidth = 0.4) +
  geom_point(data = comm_gardens, aes(x = longitude, y = latitude, fill = type, shape = garden_name), 
             color = 'blue', size = 5, alpha = 0.6) + # common garden sites
  geom_point(data = ref_genos, aes(x = longitude, y = latitude, fill = type),  # reference genotype sites
             shape = 21, color = 'black', size = 3, alpha = 0.7) + 
  geom_point(data = sat_sites, aes(x = longitude, y = latitude, fill = type), # satellite sites
             shape = 21, color = 'black', size = 3, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22, 23, 24)) + # Use unique shapes for `garden_name`
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  coord_sf(xlim = c(-130, -95), ylim = c(30, 52), crs = 4326) +
  labs(x = '', y = '', fill = 'Type', shape = 'Common garden name') +
  theme_bw()
p

ggsave('location_type.jpg', dpi = 600)


# ggplot() +
#   geom_sf(data = state_boundaries_wgs84, fill = NA, color = 'black', linewidth = 0.4) +
#   geom_point(data = comm_gardens, aes(x = longitude, y = latitude, fill = type, shape = garden_name), 
#              color = 'black', size = 5, alpha = 0.6) + # Map "Common Garden" explicitly
#   geom_point(data = ref_genos, aes(x = longitude, y = latitude, fill = type), 
#              shape = 21, color = 'black', size = 2, alpha = 0.7) + 
#   geom_point(data = sat_sites, aes(x = longitude, y = latitude, fill = type), 
#              shape = 21, color = 'black', size = 2, alpha = 0.7) +
#   scale_shape_manual(values = c(21, 22, 23, 24)) + # Unique shapes for `garden_name`
#   scale_fill_manual(
#     values = c(
#       "Common Garden" = "#00AFBB",  # Blue for Common Garden
#       "Reference Genotype Locations" = "#E7B800",         # Yellow for other types
#       "Satellite Site" = "#FC4E07"          # Orange for other types
#     )
#   ) +
#   coord_sf(xlim = c(-130, -75), ylim = c(30, 52), crs = 4326) +
#   labs(x = '', y = '', fill = 'Type', shape = 'Common Garden Name') +
#   theme_bw()
# 
# 
# ggplot() +
#   geom_polygon(data = north_america, aes(x = long, y = lat, group = group), fill = 'white', color = 'black') +
#   geom_sf(data = state_boundaries_wgs84, fill = NA, color = 'black', linewidth = 0.4) +
#   geom_point(data = comm_gardens, aes(x = longitude, y = latitude, fill = type, shape = garden_name), 
#              color = 'black', size = 5, alpha = 0.6) + # Map "Common Garden" explicitly
#   geom_point(data = ref_genos, aes(x = longitude, y = latitude, fill = type), 
#              shape = 21, color = 'black', size = 2, alpha = 0.7) + 
#   geom_point(data = sat_sites, aes(x = longitude, y = latitude, fill = type), 
#              shape = 21, color = 'black', size = 2, alpha = 0.7) +
#   scale_shape_manual(values = c(21, 22, 23, 24)) + # Unique shapes for `garden_name`
#   scale_fill_manual(
#     values = c(
#       "Common Garden" = "#00AFBB",  # Blue for Common Garden
#       "Reference Genotype Locations" = "#E7B800", # Yellow for Reference Genotypes
#       "Satellite Site" = "#FC4E07"  # Orange for Satellite Sites
#     )
#   ) +
#   coord_sf(xlim = c(-130, -75), ylim = c(30, 52), crs = 4326) +
#   labs(x = 'Longitude', y = 'Latitude', fill = 'Type', shape = 'Common Garden Name') +
#   theme_bw()
# 
