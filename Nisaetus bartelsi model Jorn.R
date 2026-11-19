### Model for future occurrences Javan Hawk-eagle
# Libraries
library(terra)
library(geodata)
library(predicts)
library(dplyr)
library(ncdf4)
library(jsonlite)
library(raster)
library(elevatr)
library(viridis)
library(rayshader)
library(sp)
library(rnaturalearth)
library(RColorBrewer)
library(rworldmap)
library(spThin)
library(sf)
library(ggplot2)
library(dbscan)
library(rgbif)
library(landscapemetrics)
# set working directory
setwd("C:/Master Biology/WEC")
# plot 1 ------------------------------------------------------------------
# Load species data
eagle_data <- sp_occurrence("Nisaetus","bartelsi", geo=TRUE, removeZeros=TRUE, args = c("Year=2010,2026"))
eagle_data <- eagle_data %>%
  dplyr::select(species, lat, lon) %>%
  distinct() %>%
  filter(!is.na(lat),
         !is.na(lon),
         between(lat, -8.8,-5.9),
         between(lon, 105, 114.6))
eagle_sf <- st_as_sf(
  eagle_data,
  coords=c("lon", "lat"),
  crs=4326)
eagle_data <- eagle_data %>% distinct()
# Download polygon for country
indonesia <- geodata::gadm(
  country = "IDN",
  level = 1,
  path = tempdir() 
  ) |> st_as_sf()
java_provinces <- indonesia %>%
  filter(NAME_1 %in% c(
    "Banten",
    "Jawa Barat",
    "DKI Jakarta",
    "Jawa Tengah",
    "DI Yogyakarta",
    "Jawa Timur"
  ))
java_provinces1 <- java_provinces %>%
  mutate(
    Population = case_when(
      NAME_1 %in% c("Banten", "Jawa Barat",
                    "Jakarta Raya", "DKI Jakarta") ~ "West",
      NAME_1 %in% c("Jawa Tengah",
                    "DI Yogyakarta", "Yogyakarta") ~ "Central",
      NAME_1 == "Jawa Timur" ~ "East",
      TRUE ~ NA_character_
    ))
java_provinces1 %>%
  st_drop_geometry() %>%
  dplyr::select(NAME_1, Population)
eagle_assigned <- st_join(
  eagle_sf,
  java_provinces1 %>% dplyr::select(NAME_1, Population),
  join = st_within,
  left = TRUE
)
eagle_assigned %>%
  st_drop_geometry() %>%
  count(Population)
ggplot() +
    # Province polygons
  geom_sf(
    data = java_provinces1,
    aes(fill = Population),
    color = "black",
    linewidth = 0.3
  ) +
    # Eagle occurrence records
  geom_sf(
    data = eagle_assigned,
    aes(color = Population),
    size = 2,
    alpha = 0.8
  ) +
    # Colours
  scale_fill_manual(
    values = c(
      "West" = "#E69F00",
      "Central" = "#56B4E9",
      "East" = "#009E73"
    )
  ) +
    scale_color_manual(
    values = c(
      "West" = "blue",
      "Central" = "pink",
      "East" = "yellow"
    ),
    na.value = "red"
  ) +
    # Map extent
  coord_sf(
    xlim = c(105, 114.6),
    ylim = c(-8.8, -5.9),
    expand = FALSE
  ) +
    labs(
    title = expression(italic("Nisaetus bartelsi")),
    subtitle = "Geographical distribution of three subpopulations on Java",
    fill = "Subpopulation",
    color = "Subpopulation",
    x = "Longitude",
    y = "Latitude"
  ) +
    theme_minimal() +
    theme(
    panel.grid = element_line(color = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "italic", size = 14)
  )

# 2e plot -----------------------------------------------------------------
eagle_data <- eagle_data %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  distinct(lon, lat, .keep_all = TRUE)
eagle_sf <- st_as_sf(
  eagle_data,
  coords = c("lon", "lat"),
  crs = 4326
)
nrow(eagle_sf)
eagle_utm <- st_transform(eagle_sf, 32749)
coords <- st_coordinates(eagle_utm)
set.seed(123)
eagle_clusters <- dbscan(
  coords,
  eps = 25000,
  minPts = 1
)
eagle_utm$Population <- as.factor(eagle_clusters$cluster)
table(eagle_utm$Population)
ggplot() +
geom_sf(
    data = java_provinces,
    color = "black",
    linewidth = 0.3
  ) +
  geom_sf(
    data = eagle_utm,
    aes(color = Population),
    size = 2,
    alpha = 0.8
  ) +
  
  scale_color_viridis_d(
    option = "turbo",
    name = "Population"
  ) +
  
  labs(
    title = expression(italic("Nisaetus bartelsi")),
    subtitle = "DBSCAN spatial clustering (5 km distance)",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  coord_sf(crs = 4326) +
  
  theme_minimal() +
  
  theme(
    legend.position = "right",
    plot.title = element_text(face = "italic")
  )

populationsjava <- data.frame (id = 1:7, name = c("Ujung Kulon",
                                                  "Halimun-Salak & Gede-Pangrango",
                                                  "Ciremai & West-Central Mts",
                                                  "Dieng & Slamet",
                                                  "Merapi-Merbabu & Lawu",
                                                  "Bromo-Tengger-Semeru",
                                                  "Ijen & Baluran-Meru Betiri"
),
lon = c(105.35, 106.75, 108.40, 109.55, 110.80, 112.95, 114.20),
lat = c(-6.75,  -6.72,  -6.90,  -7.22,  -7.55,  -8.00,  -8.12),
# Rough size classes: 10, 50, 100, 250 (Total across Java ~500-1000 individuals)
size_class = c(50, 250, 100, 100, 100, 250, 100),
min_est = c(20, 180, 60, 50, 50, 150, 60),
max_est = c(80, 350, 150, 160, 150, 350, 160),
matrix_note = c(
  "Isolated coastal forest",
  "High forest cover; connected mountain block",
  "Fragmented by lowland agriculture",
  "Stepping-stone mountain peaks",
  "Heavily isolated by agricultural matrix",
  "Large montane refuge",
  "Eastern montane & lowland forest block"
)
)

# 3. Define Migration / Gene Flow Links between neighboring populations
# Note: Juvenile dispersal occurs across forest gaps, but heavily cultivated lowlands restrict movement
migration <- data.frame(
  from = c(1, 2, 3, 4, 5, 6),
  to   = c(2, 3, 4, 5, 6, 7),
  # Rough migration rate classes (migrants per generation)
  mig_class = c(
    "0.5 - 1 / gen",
    "1 - 4 / gen",
    "0.5 - 1 / gen",
    "0.1 - 0.5 / gen",
    "0.1 - 0.5 / gen",
    "0.5 - 1 / gen"
  ),
  # Uncertainty scenarios (Min - Max migrants/generation)
  mig_scenario = c(
    "0.2 - 1.5",
    "1.0 - 5.0",
    "0.3 - 1.5",
    "0.05 - 0.8",
    "0.05 - 0.8",
    "0.3 - 1.5"
  )
)

# Join coordinates to migration table for plotting segments
mig_lines <- migration %>%
  left_join(populationsjava %>% dplyr::select(id, lon, lat), by = c("from" = "id")) %>%
  rename(x = lon, y = lat) %>%
  left_join(populationsjava %>% dplyr::select(id, lon, lat), by = c("to" = "id")) %>%
  rename(xend = lon, yend = lat) %>%
  mutate(
    mid_x = (x + xend) / 2,
    mid_y = (y + yend) / 2,
    mig_class = factor(mig_class, levels = c(
      "0.1 - 0.5 / gen",
      "0.5 - 1 / gen",
      "1 - 4 / gen",
      "> 4 / gen"
    ))
  )

# 4. Load Geographic Layer (Indonesia / Java outline)
library(rnaturalearthdata)
indonesia <- ne_countries(scale = "medium", country = "Indonesia", returnclass = "sf")

# Optional: If you have the 'geodata' package installed, you can load elevation 
# to visually show the mountain "sky islands" vs. the unsuitable lowland matrix:
# library(geodata)
# elev <- elevation_30s(country = "IDN", path = tempdir())
# elev_java <- terra::crop(elev, terra::ext(105, 115, -9, -5.5))
# elev_df <- as.data.frame(elev_java, xy = TRUE)

# 5. Plot the Population & Gene Flow Map
ggplot() +
  # Base map of Java
  geom_sf(data = indonesia, fill = "#f2efe9", color = "#7f7f7f", linewidth = 0.4) +
  
  # Migration arrows / gene flow lines
  geom_curve(
    data = mig_lines,
    aes(x = x, y = y, xend = xend, yend = yend, color = mig_class, linewidth = mig_class),
    curvature = -0.15,
    alpha = 0.85
  ) +
  
  # Population points scaled by rough size class
  geom_point(
    data = populationsjava,
    aes(x = lon, y = lat, size = size_class),
    color = "#1b4332",
    fill = "#40916c",
    shape = 21,
    stroke = 1.5,
    alpha = 0.9
  ) +
  
  # Labels with Population Name + Size Class + Uncertainty [Min-Max]
  geom_text(
    data = populationsjava,
    aes(
      x = lon, y = lat,
      label = paste0(name, "\nClass: ~", size_class, " [", min_est, "-", max_est, "]")
    ),
    vjust = -1.2,
    size = 3,
    fontface = "bold",
    color = "#1d3557"
  ) +
  
  # Customize scales and legends
  scale_size_continuous(
    name = "Population Size Class\n(Individuals)",
    breaks = c(50, 100, 250),
    range = c(4, 12)
  ) +
  scale_color_manual(
    name = "Migration Rate\n(Migrants / Gen)",
    values = c(
      "0.1 - 0.5 / gen" = "#d95f02",
      "0.5 - 1 / gen"   = "#7570b3",
      "1 - 4 / gen"     = "#1b9e77",
      "> 4 / gen"       = "#006d2c"
    )
  ) +
  scale_linewidth_manual(
    name = "Migration Rate\n(Migrants / Gen)",
    values = c(
      "0.1 - 0.5 / gen" = 0.8,
      "0.5 - 1 / gen"   = 1.4,
      "1 - 4 / gen"     = 2.2,
      "> 4 / gen"       = 3.0
    )
  ) +
  
  # Zoom into Java
  coord_sf(xlim = c(105.0, 114.8), ylim = c(-8.9, -5.8), expand = FALSE) +
  
  # Titles and Model Metadata (Matrix Suitability & Linear Geography Note)
  labs(
    title = "Javan Hawk-Eagle (Nisaetus bartelsi) Subpopulations & Gene Flow Model",
    subtitle = "Rough size classes with [Min–Max] uncertainty scenarios across major montane forest blocks",
    caption = paste(
      "Matrix Suitability: LOW in lowland agricultural/urban areas (major dispersal barrier); HIGH only along forested mountain ridges.\n",
      "Spatial Structure Note: Java's elongated shape creates a 1D stepping-stone distribution, causing Fst to increase rapidly with distance."
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "#dceef8", color = NA),
    panel.grid.major = element_line(color = "#ffffff", linewidth = 0.3),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.caption = element_text(hjust = 0, size = 8.5, face = "italic")
  )

# Plot 3 ------------------------------------------------------------------
library(readr)
# Create a local folder to store downloaded rasters
dir.create("java_rasters", showWarnings = FALSE)

# 1. Get Java polygon boundary (Province level = 1, filter to Java provinces)
idn_prov <- gadm(country = "IDN", level = 1, path = "java_rasters")
java_provs <- c("Banten", "Jakarta Raya", "Jawa Barat", "Jawa Tengah", "Yogyakarta", "Jawa Timur")
java_poly <- idn_prov[idn_prov$NAME_1 %in% java_provs, ]

# 2. Download Elevation (~1 km) & compute Slope + Aspect
elev <- elevation_30s(country = "IDN", path = "java_rasters") |> crop(java_poly, mask = TRUE)
terrain_vars <- terrain(elev, v = c("slope", "aspect"), unit = "degrees")

# 3. Download WorldClim Bioclimatic layers (0.5 arc-min = ~1 km resolution)
clim <- worldclim_country(country = "IDN", var = "bio", res = 0.5, path = "java_rasters") |>
  crop(java_poly, mask = TRUE)
# Keep key variables: BIO1 (Annual Mean Temp), BIO12 (Annual Precip), BIO15 (Precip Seasonality)
clim_sub <- clim[[c(1, 12, 15)]]
names(clim_sub) <- c("mean_temp", "annual_precip", "precip_seasonality")

# 4. Download Fractional Tree Cover & Human Footprint (~1 km)
trees <- landcover(var = "trees", path = "java_rasters") |> crop(java_poly, mask = TRUE)
footprint_idx <- footprint(year = 2009, path = "java_rasters") |> crop(java_poly, mask = TRUE)

# 5. Align all layers to the exact same grid/resolution as elevation
trees_res <- resample(trees, elev, method = "bilinear")
foot_res  <- resample(footprint_idx, elev, method = "bilinear")
clim_res  <- resample(clim_sub, elev, method = "bilinear")

# 6. Combine into a single ready-to-model raster stack
env_stack <- c(elev, terrain_vars, clim_res, trees_res, foot_res)
names(env_stack) <- c("elevation", "slope", "aspect", "mean_temp", 
                      "annual_precip", "precip_seasonality", "tree_cover", "human_footprint")

# Save stack to disk so you don't have to re-download
writeRaster(env_stack, "java_rasters/java_env_stack.tif", overwrite = TRUE)
target_crs <- "EPSG:32749"
env_utm <- project(env_stack, target_crs)
pts_utm <- project(eagle_data, target_crs)
bg_pts <- spatSample(env_ut)
newmap <- geodata::getMap(resolution = "low")
plot(newmap)
extent(c(-8.8,-5.9,105,114.6))
plot(newmap, xlim = extent[c(1,2)], ylim = extent[c(3,4)])
pop_data <- read.table("nisaetus_bartelsi_populations.txt", header=TRUE)
points(pop_data$x, pop_data$y)
pop_data$pointcex <- ifelse(pop_data$size<100,1,
                            ifelse(pop_data$size<250,2))

plot(newmap,
     xlim = extent[c(1,2)],
     ylim = extent[c(3,4)],
     border="pink")
points(pop_data$x, pop_data$y, cex=pop_data$pointcex, col = "blue", pch=19)
text(pop_data$x, pop_data$y, pop_data$size, cex=1)
