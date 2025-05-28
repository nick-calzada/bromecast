library(readr)
library(rnaturalearth)
library(USA.state.boundaries)
library(tidyverse)

data("state_boundaries_wgs84")
world <- ne_countries(scale = "medium", returnclass = "sf")
north_america <- world[which(world$continent == "North America"),]

# Load in locations and codes
genotype_locations <- 
  read_csv(here('01_data', 'data', 'bioclim_and_other', 'BioclimateOfOrigin_AllGenotypes.csv')) %>%
  dplyr::select("source" = site_code, genotype, "longitude" = lon, "latitude" = lat)

# Join with key matrix for connecting with genotype matrix 
genotype_codes <- 
  read_csv(here('01_data', 'data', 'bioclim_and_other', 'BRTEcg_genotypesCode.csv')) %>%
  left_join(genotype_locations)  %>%
  filter(!is.na(longitude)) %>%
  arrange(SNPmatrix_column) %>%
  st_as_sf(.,
           coords = c("longitude", "latitude"), 
           crs = 4326) %>%
  ### Keep lat, long as columns too 
  dplyr::mutate(longitude = sf::st_coordinates(.)[,1],
                latitude = sf::st_coordinates(.)[,2])

# Load in genetic dissimilarities...
load(here('01_data', 'data', 'snps', 'invaded_Z_and_D.RData'))

# Convert genetic distance matrix into a dataframe
d_long <- as.data.frame(as.table(invaded_D)) %>%
  rename(start = Var1, stop = Var2, gd = Freq) 

# Add in 'start' coordinates
d_long <- left_join(d_long, genotype_codes, by = c('start' = 'source')) 
d_long <- d_long %>% rename(start_lon = longitude, start_lat = latitude, PopNum_start = PopNum)

# Add in 'stop' coordinates 
d_long <- left_join(d_long, genotype_codes, by = c('stop' = 'source'))
d_long <- d_long %>% rename(end_lon = longitude, end_lat = latitude, PopNum_end = PopNum)
d_long <- d_long %>% dplyr::select(start, stop, start_lon, start_lat, end_lon, end_lat, gd, PopNum_start, PopNum_end) %>% filter(start != stop) # Do not want pairs of self-self

# Gather unique pairings 
u_pairs <- d_long %>% 
  rowwise() %>%
  mutate(
    pair_start = min(start, stop),
    pair_stop = max(start, stop)
  ) %>%
  ungroup() %>% 
  distinct(pair_start, pair_stop, .keep_all = TRUE)

# Calculate distance between start and stop coordinates 
u_pairs$distance_km <- distHaversine(
  cbind(u_pairs$start_lon, u_pairs$start_lat), 
  cbind(u_pairs$end_lon, u_pairs$end_lat)
) / 1000

# Separate into groups based on gd 
under_100 <- u_pairs %>% filter(gd < 100)
in_between <- u_pairs %>% filter(gd >= 100 & gd <= 800)
over_800 <- u_pairs %>% filter(gd > 800)


# Plot
a <- 1
lw1 <- 0.15
lw2 <- 0.3

p <- ggplot() +
  geom_sf(data = north_america,
          fill = NA,
          color = 'black') +
  geom_sf(
    data = state_boundaries_wgs84,
    fill = NA,
    color = 'black',
    linewidth = 0.3,
    alpha = 0.2
  ) +
  geom_segment(
    data = over_800,
    aes(
      x = start_lon,
      y = start_lat,
      xend = end_lon,
      yend = end_lat,
      color = '> 800'
    ),
    linewidth = lw1,
    alpha = a
  ) +
  
  geom_segment(
    data = under_100,
    aes(
      x = start_lon,
      y = start_lat,
      xend = end_lon,
      yend = end_lat,
      color = '< 100'
    ),
    linewidth = lw2,
    alpha = a
  ) +
  geom_point(
    data = u_pairs,
    aes(x = start_lon, y = start_lat),
    color = 'black'
  ) +
  scale_color_manual(values = c('> 800' = 'lightblue', '< 100' = 'deeppink')) +
  labs(x = '', y = '', color = 'Genetic Distance') +
  coord_sf(
    xlim = c(-125, -100),
    ylim = c(34, 51),
    crs = 4326,
    expand = FALSE
  ) +
  theme_bw()

p
