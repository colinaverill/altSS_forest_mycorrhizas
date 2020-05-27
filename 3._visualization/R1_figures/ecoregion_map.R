rm(list=ls())
library(rgdal)
library(rnaturalearth)
library(rmapshaper)
library(broom)
library(RColorBrewer)
library(tidyverse)
source('paths.r')

# Load shapefile.----
shp <- readOGR(dsn = EPA_L2_ecoregions_raster.path, stringsAsFactors = F)

#Specify ecoregions of interest.----
subset <- shp[shp$NA_L2NAME %in% c("ATLANTIC HIGHLANDS",
                                   "MIXED WOOD PLAINS",
                                   "MIXED WOOD SHIELD",
                                   "OZARK/OUACHITA-APPALACHIAN FORESTS",
                                   "SOUTHEASTERN USA PLAINS"),]

# Define color palette.----
col <- c("#B26A4A","#B2A877","#62B25E","#40B29A","#ef7e6b")


# Convert CRS to WGS84.-----
subset <- spTransform(subset, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Simplify the shapefile with 'ms_simplify', keeping 5% (speed).----
subset_simplified <- ms_simplify(subset, keep = 0.05, keep_shapes = T)

USA <- ne_countries(scale = 10, country = 'United States of America')

# Crop by US borders.----
subset_simplified_cropped <- raster::crop(subset_simplified, USA)

subset_simplified_df <- tidy(subset_simplified_cropped, region = "NA_L2NAME") %>% 
  data.frame() %>% 
  rename(NA_L2NAME = id) %>% 
  mutate(NA_L2NAME = str_to_title(NA_L2NAME))

# Define world outline
world <- ne_countries(scale = "small", returnclass = "sf")
USA_plotdata <- ne_countries(country = 'United States of America', scale = 10, returnclass = "sf")


#Drop the plot.----
png('test.png', width = 8, height = 4, units = 'in', res = 300)
ggplot() +
  geom_sf(data = USA_plotdata, size = 0.2) +
  theme_minimal() +
  coord_sf(xlim = c(-100, -65), ylim = c(23, 50), expand = FALSE) +         #subset plotting region to eastern US.
  geom_polygon(data = subset_simplified_df, 
               aes(x = long, y = lat, group = group, fill = NA_L2NAME),
               size = 1) +
  theme(axis.title = element_blank()) +                                     #drop X-Y axis lbels.
  theme(panel.grid = element_blank(),axis.text = element_blank()) +         #drop gridlines and axis unit labels
  scale_fill_manual(values = col, name = "Ecoregion") 
dev.off()
