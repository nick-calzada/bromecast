library(readr)
library(rnaturalearth)
library(tidyverse)

# Load in country geometries for plotting
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
Asia <- world[which(world$continent == "Asia"),]
Africa <- world[which(world$continent == 'Africa'),]

# Load in native range meta data
native <- read_csv(here('01_data', 'data', 'bioclim_and_other', 'BRTEclim.csv')) %>% filter(range == 'native')

# Load in previoulsy constructed D matrix for native range
load(here('01_data', 'data', 'snps', 'native_Z_and_D.RData'))

# IBS.id column contains snp matrix columns
cols <- read_csv(here('01_data', 'data', 'bioclim_and_other', '307tips.csv'))
native <- left_join(native, cols, by = 'PopNum') %>% filter(!is.na(Longitude))

d_long <- as.data.frame(as.table(native_D)) %>%
  rename(start = Var1, stop = Var2, gd = Freq) 

# Add in 'start' coordinates
d_long <- left_join(d_long, native, by = c('start' = 'NewSiteCode')) 
d_long <- d_long %>% rename(start_lon = Longitude, start_lat = Latitude, PopNum_start = PopNum)

# Add in 'stop' coordinates 
d_long <- left_join(d_long, native, by = c('stop' = 'NewSiteCode'))
d_long <- d_long %>% rename(end_lon = Longitude, end_lat = Latitude, PopNum_end = PopNum)
d_long <- d_long %>% dplyr::select(start, stop, start_lon, start_lat, end_lon, end_lat, gd, PopNum_start, PopNum_end) %>% filter(start != stop) # do not want pairs of self-self

all_unique_pairs <- d_long %>% 
  rowwise() %>%
  mutate(
    pair_start = min(start, stop),
    pair_stop = max(start, stop)
  ) %>%
  ungroup() %>% 
  distinct(pair_start, pair_stop, .keep_all = TRUE)

# Split into genetic distance groups
over_800 <- all_unique_pairs %>% filter(gd > 800)
under_100 <- all_unique_pairs %>% filter(gd < 100)

# Plot
a <- 1 # alpha
lw1 <- 0.2 # lw of pairs with a big gd
lw2 <- 0.3 # lw of pairs with a low gd 


p <- ggplot() +
  geom_sf(
    data = Europe,
    fill = 'white',
    color = 'black',
    linewidth = 0.3
  ) +
  geom_sf(
    data = Asia,
    fill = 'white',
    color = 'black',
    linewidth = 0.3
  ) +
  geom_sf(
    data = Africa,
    fill = 'white',
    color = 'black',
    linewidth = 0.3
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
    # <- Mapped as a string
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
    # <- Mapped as a string
    linewidth = lw2,
    alpha = a
  ) +
  
  geom_point(data = all_unique_pairs,
             aes(x = start_lon, y = start_lat),
             color = 'black') +
  
  scale_color_manual(values = c('> 800' = 'lightblue', '< 100' = 'deeppink')) +
  labs(x = '', y = '', color = 'Genetic Distance') +
  # theme(
  #   panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  #   panel.background = element_blank()
  # ) +
  theme_bw() +
  coord_sf(xlim = c(-15, 80),
           ylim = c(28, 67),
           expand = FALSE)

p
