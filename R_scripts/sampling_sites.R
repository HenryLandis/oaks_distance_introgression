library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)

# Country basemaps.
mexico <- ne_countries(country = 'mexico',
                       scale = 'medium',
                       returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

# Shapefiles for species ranges.
spp <- c('querariz',
         'quergamb',
         'quergris',
         'querhava',
         'quermohr',
         'quermueh',
         'querturb')
rlist <- sapply(spp, function(x) {
  st_read(paste0(
    "C:/Users/Henry/Desktop/WSU courses/Thesis/Analyses/Sampling sites/",x,".shp"),
    crs = 4326) %>%
    st_cast("MULTIPOLYGON")
}, simplify = FALSE)

# Coordinates for sampling sites, the centroids of all
# individuals at the site.
sites <- data.frame(longitude = c(-114.3291, # GLV
                                  -107.8563, # SJN
                                  -105.3454, # LIN
                                  -111.3891, # PAY
                                  -105.2955, # SFE
                                  -108.2239, # GIL
                                  -110.6225), # GLO
                    latitude = c(35.11573, 
                                 37.46103, 
                                 33.54853,
                                 34.05656, 
                                 35.29932,
                                 32.89932, 
                                 33.60294))
myshapes <- c(19, 19, 1, 1, 19, 1, 1)
myfills <- c("red", "red", "blue", "blue", "red",
             "blue", "blue")
names <- c("GOLDEN VALLEY
Q. turbinella", #oatman
           "SAN JUAN
Q. gambelii", #durango
           "LINCOLN
Q. gambelii
Q. grisea", #lincoln
           "PAYSON
Q. gambelii
Q. grisea
Q. turbinella", #mt. ord
           "SANTA FE
Q. grisea", #villanueva
           "GILA
Q. gambelii
Q. grisea", #gila
           "GLOBE
Q. grisea
Q. turbinella" #globe
) #gmnp

# Species ranges only.
ggplot() +
  geom_sf(data = mexico, fill = 'lightgray') +
  geom_sf(data = states, fill = 'lightgray') +
  geom_sf(data = rlist[["quergamb"]], alpha = 0.4, 
          fill = "#d95f02", linetype = "blank") +
  geom_sf(data = rlist[["quergris"]], alpha = 0.4, 
          fill = "#7570b3", linetype = "blank") +
  geom_sf(data = rlist[["querturb"]], alpha = 0.4, 
          fill = "#1b9e77", linetype = "blank") +
  coord_sf(crs = 4326, 
           xlim = c(-123.5, -100), 
           ylim = c(28, 41)) +
  annotate('rect', xmin = -120, xmax = -115, ymin = 38, ymax = 40.5,
           fill = 'white', color = 'black') +
  annotate("text",
           x = c(-117.5, -117.5, -117.5), 
           y = c(40, 39.25, 38.5),
           label = c("Quercus gambelii",
                     "Quercus grisea",
                     "Quercus turbinella"),
           color = c("#d95f02", "#7570b3", "#1b9e77")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank())

# With the sampling sites.
ggplot() +
  geom_sf(data = mexico, fill = 'lightgray') +
  geom_sf(data = states, fill = 'lightgray') +
  geom_sf(data = rlist[["quergamb"]], alpha = 0.4, 
          fill = "#d95f02", linetype = "blank") +
  geom_sf(data = rlist[["quergris"]], alpha = 0.4, 
          fill = "#7570b3", linetype = "blank") +
  geom_sf(data = rlist[["querturb"]], alpha = 0.4, 
          fill = "#1b9e77", linetype = "blank") +
  geom_point(data = sites, aes(x = longitude, y = latitude),
             size = 3, stroke = 1, color = "black",
             shape = myshapes, fill = myfills) +
  geom_label_repel(data = sites, aes(
    x = longitude, y = latitude, label = names), 
    size = 3, max.overlaps = Inf,
    nudge_x = c(-3, 4, -2.5, -6, 4, -4, -1),
    nudge_y = c(-1, 1, -3.5, -8, -0.5, -3, 4)) +
  coord_sf(crs = 4326, 
            xlim = c(-123.5, -100), 
            ylim = c(28, 41)) +
  annotate('rect', xmin = -120, xmax = -115, ymin = 38, ymax = 40.5,
           fill = 'white', color = 'black') +
  annotate("text",
           x = c(-117.5, -117.5, -117.5), 
           y = c(40, 39.25, 38.5),
           label = c("Quercus gambelii",
                     "Quercus grisea",
                     "Quercus turbinella"),
           color = c("#d95f02", "#7570b3", "#1b9e77")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank())

# Plot only the sympatric sites and simplify labels 
# (for slide comparing to admixture plots).
sympatric_sites <- data.frame(longitude = c(-105.3454, # LIN
                                            -111.3891, # PAY
                                            -108.2239, # GIL
                                            -110.6225), # GLO
                    latitude = c(33.54853,
                                 34.05656,
                                 32.89932, 
                                 33.60294))
sympatric_shapes <- c(1, 1, 1, 1)
sympatric_names <- c('Lincoln', 'Payson', 'Gila', 'Globe')
ggplot() +
  geom_sf(data = mexico, fill = 'lightgray') +
  geom_sf(data = states, fill = 'lightgray') +
  geom_sf(data = rlist[["quergamb"]], alpha = 0.4, 
          fill = "#d95f02", linetype = "blank") +
  geom_sf(data = rlist[["quergris"]], alpha = 0.4, 
          fill = "#7570b3", linetype = "blank") +
  geom_sf(data = rlist[["querturb"]], alpha = 0.4, 
          fill = "#1b9e77", linetype = "blank") +
  geom_point(data = sympatric_sites, 
             aes(x = longitude, y = latitude),
             size = 2, stroke = 1, color = "black",
             shape = sympatric_shapes) +
  geom_label_repel(data = sympatric_sites, aes(
    x = longitude, y = latitude, label = sympatric_names), 
    size = 3, max.overlaps = Inf,
    nudge_x = c(0, 0, 0, -0.25),
    nudge_y = c(-0.75, -0.758, -0.75, -0.75)) +
  coord_sf(crs = 4326, 
           xlim = c(-123.5, -100), 
           ylim = c(28, 41)) +
  annotate('rect', xmin = -120, xmax = -115, ymin = 38, ymax = 40.5,
           fill = 'white', color = 'black') +
  annotate("text",
           x = c(-117.5, -117.5, -117.5), 
           y = c(40, 39.25, 38.5),
           label = c("Quercus gambelii",
                     "Quercus grisea",
                     "Quercus turbinella"),
           color = c("#d95f02", "#7570b3", "#1b9e77")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank())
