######################################################################################### 
#                                                                                       #
#                 Calculating microclimate connectivity in the Amazon                   #
#                                                                                       #
######################################################################################### 
#
# This script calculates microclimate connectivity and is a direct adaptation from the R code
# by McGuire et al. (2016) (www.pnas.org/cgi/doi/10.1073/pnas.1602817113). The source
# can be accessed in https://github.com/JennyMcGuire/ClimateConnectivity
#
# The script is divided into three parts:
# (1) Preparing data
# (2) Creating temperature patches, i.e. the basic units for climate connectivity analysis
# (3) Calculating climate connectivity
# 
# The exact same workflow applies for macroclimate connectivity. Input raster data must be 
# changed accordingly and reprojected to match the microclimate rasters.
#
# -----------------------------------------------------------------------------
# Example data
# -----------------------------------------------------------------------------
#
# A small example dataset is included with this repository to demonstrate the workflow and
# allow users to run the analysis without access to the full Amazon-wide datasets used in the study.
#
# The example dataset consists of a spatial subset of the study area and is provided for testing 
# and reproducibility purposes only. Results generated from the example data will differ from those 
# reported in the manuscript because the published analyses were conducted using the complete datasets.
#
# No example dataset for macroclimate is provided.
#
# -----------------------------------------------------------------------------
# Full study data
# -----------------------------------------------------------------------------
#
# Due to file size constraints, the full datasets are not distributed with this repository.
# The full analyses used:
#
# Inputs:
# - microclimate present temperature raster (2013-2022) 
#   Historic mean in https://doi.org/10.23729/fd-4585bc48-aba4-3e58-8d45-b0c9dc2a6492
#
# - microclimate future temperature raster (SSP245, 2093-2100)
#   SSP245 mean in https://doi.org/10.23729/fd-4585bc48-aba4-3e58-8d45-b0c9dc2a6492
#
# - microclimate future temperature raster (SSP585, 2093-2100)
#  SSP585 mean in https://doi.org/10.23729/fd-4585bc48-aba4-3e58-8d45-b0c9dc2a6492
#
# - cropped elevation in  (NASADEM) # NASADEM https://doi.org/10.5069/G93T9FD9
#
# - cropped river polygons  # GRWL https://doi.org/10.1126/science.aat0636
#
# -----------------------------------------------------------------------------
# Outputs:
# - temperature patch polygons
# - microclimate connectivity raster for SSP245
# - microclimate connectivity raster for SSP585
# - microclimate connectivity CSV
# - microclimate connectivity shapefile
#
# All rasters and vector layers must:
# - Be projected to EPSG:5880
# - Cover the same spatial extent
# - Have matching resolution and alignment
#
# All lines that save outputs are initially commented.

######################################################################################### 
#                                                                                       #
#                                     1. Preparations                                   #
#                                                                                       #
######################################################################################### 

# Load libraries
library(terra)
library(sf)
library(igraph)
library(dplyr)

# Define data folders
data.dir <- "data"
result.folder <- "path"

# Create the results folder if it does not exist yet
#dir.create(result.folder, showWarnings = FALSE)

# Analysis settings
max_elevation <- 600
river_width_threshold <- 1000
patch_merge_distance <- 2000
minimum_patch_area <- 10 # NOTE! For example data, change to 0 

## Full-data preparation workflow (not executable without source datasets)

# The code assumes files are cropped to match the Amazon biome boundaries

#micro_t <- rast(file.path(data.dir, "microclimate_file")) 
#elev <- rast(file.path(data.dir, "elevation_file")) 
#rivers <- st_read(file.path(data.dir, "rivers_file"))
#micro_245 <- rast(file.path(data.dir, "ssp245_file"))
#micro_585 <- rast(file.path(data.dir, "ssp585_file")) 

# Reproject data
#micro_t <- project(micro_t, "EPSG:5880", method = "bilinear")
#elev <- project(elev, "EPSG:5880", method = "bilinear")
#rivers <- project(rivers, "EPSG:5880")
#micro_245 <- project(micro_245, micro_t, method = "bilinear")
#micro_585 <- project(micro_585, micro_t, method = "bilinear")

## Remove high-elevation pixels > 600 meters in elevation and mask them out from temperature rasters

#elev[elev > max_elevation] <- NA

# Next, mask high-elevation pixels from the temperature map
#micro_t_filtered <- mask(micro_t, elev)
#micro_245_filtered <- mask(micro_245, elev)
#micro_585_filtered <- mask(micro_585, elev)

# Save rasters for future
#writeRaster(micro_t_filtered, file.path(data.dir, "micro_t_filtered.tif"))
#writeRaster(micro_245_filtered, file.path(data.dir, "micro_245_filtered.tif"))
#writeRaster(micro_585_filtered, file.path(data.dir, "micro_585_filtered.tif"))

############################################################################################# 
#                                                                                           #
# 2. Create the basic units for climate connectivity calculation, the "temperature patches" #
#                                                                                           #
#############################################################################################

## Example datasets

micro_t_filtered <- rast(
  file.path(data.dir, "micro_t_test_data.tif")
)

micro_245_filtered <- rast(
  file.path(data.dir, "micro_245_test_data.tif")
)

micro_585_filtered <- rast(
  file.path(data.dir, "micro_585_test_data.tif")
)

rivers <- st_read(
  file.path(data.dir, "river_test_data.gpkg")
)

## Analysis workflow

# "temp" is the working file that can refer to the original microclimate, macroclimate or example data
# of current temperature conditions.

# Create working file
temp <- micro_t_filtered # Change according to your temperature data

## Reclassify temperature pixels with 0.5 degree intervals and create polygons

# Get minimum and maximum temperatures
rmin <- floor(minmax(temp)[1,1])
rmax <- ceiling(minmax(temp)[2,1])

# Create 0.5 interval breakpoints (always at .0 and .5)
breaks <- seq(rmin, rmax, by = 0.5)

# Build reclassification matrix
reclass_matrix <- cbind(breaks[-length(breaks)], breaks[-1], 1:(length(breaks)-1))

# Reclassify temperature raster
temp_classes <- classify(temp, reclass_matrix, include.lowest = TRUE)

# Convert raster to polygons
# You need to first aggregate and then disaggregate, otherwise you lose information
temp_polygons <- as.polygons(temp_classes, aggregate = TRUE)

# Disaggregate
temp_polygons_disagg <- disagg(temp_polygons)

# Add unique polygon ID column
temp_polygons_disagg[["ID"]] <- 1:nrow(temp_polygons_disagg)

# Rename value column
names(temp_polygons_disagg) <- c("Class","ID")

# Account for large rivers 

# Remove narrow river segments (when working with original data)
#rivers <- rivers[rivers$width_m >= river_width_threshold]

# Dissolve into one line polygon (when working with original data)
#rivers <- dissolve(rivers)

classes <- unique(temp_polygons_disagg$Class)

# Convert from st to sf
sf_polygons <- st_as_sf(temp_polygons_disagg)
sf_use_s2(FALSE)

sf_polygons$patch_id <- NA_integer_
next_patch_id <- 1

# Merge polygons belonging to the same class based on maximum distance
# Prevent merging across rivers wider than the defined thershold

# NOTE!
# Patch creation can be computationally intensive for large study areas.
# Runtime increases substantially with the number of temperature polygons and classes.

for (class in classes) {
  
  message("Processing class: ", class)
  
  # Subset class
  p <- sf_polygons[sf_polygons$Class == class, ]
  n <- nrow(p)
  
  # Find all neighbors within distance
  nb_dist <- st_is_within_distance(p, p, dist = patch_merge_distance)
  
  # Initialize allowed neighbors list
  nb_allowed <- vector("list", n)
  
  # Loop over each polygon and find if a river prevents merging
  
  message("Now looping through each polygon and creating lines between origin and candidates")
  
  for (i in seq_len(n)) {
    candidates <- nb_dist[[i]]
    
    keep <- logical(length(candidates))
    
    for (k in seq_along(candidates)) {
      j <- candidates[k]
      
      line <- st_nearest_points(p[i, ], p[j, ])
      keep[k] <- !st_intersects(line, rivers, sparse = FALSE)[1, 1]
      
    }
    
    nb_allowed[[i]] <- candidates[keep]
  }
  
  message("Building graph")
  g <- graph_from_adj_list(nb_allowed, mode = "all")
  cl <- components(g)
  
  # Assign global patch IDs
  message("Assigning patch IDs")
  p$patch_id <- cl$membership + next_patch_id - 1
  
  # Update the full object
  message("Updating original polygons")
  sf_polygons$patch_id[sf_polygons$Class == class] <- p$patch_id
  
  # Increment the counter
  next_patch_id <- max(sf_polygons$patch_id, na.rm = TRUE) + 1
}

# Dissolve polygons based on patch ID
polygons_dissolved <- sf_polygons |>
  dplyr::group_by(patch_id) |>
  dplyr::summarise(Class = mean(Class))

# Add area column
polygons_dissolved$area_m2 <- round(st_area(polygons_dissolved),0)

# Convert area to km2
polygons_dissolved <- polygons_dissolved %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6)

# Remove polygons with area smaller than the minimum area
polygons_filtered <- polygons_dissolved %>%
  filter(area_km2 >= minimum_patch_area)

# Save
#st_write(polygons_filtered, file.path(data.dir, "polygons_filtered.shp"))

############################################################################################# 
#                                                                                           #
#                           3. Calculate climate connectivity                               #
#                                                                                           #
#############################################################################################

## Extract current and future temperatures for polygons (i.e. patches)

# First, combine rasters
rasters <- c(micro_t_filtered, micro_245_filtered, micro_585_filtered)

# Extract mean temperature values of current, ssp245 and ssp585 temperatures per polygon
mean_values <- terra::extract(rasters, vect(polygons_filtered), fun = mean, na.rm = TRUE)

#Rename the current, ssp245 and ssp585 columns
names(mean_values)[2:4] <- c("current", "ssp245", "ssp585")

#Add mean temperatures as new columns in polygons
polygons_filtered <- cbind(polygons_filtered, round(mean_values[,-1],1))

# Create new ID column for final patches
polygons_filtered$new_patch_id <- seq_len(nrow(polygons_filtered))

# Find neighbors for each polygon
neighbors_list <- st_relate(polygons_filtered, polygons_filtered, pattern = "F***1****")

# Make into dataframe
neighbor_table <- do.call(rbind, lapply(seq_along(neighbors_list), function(i) {
  if (length(neighbors_list[[i]]) == 0) return(NULL)
  data.frame(
    src = polygons_filtered$new_patch_id[i],
    nbr = polygons_filtered$new_patch_id[neighbors_list[[i]]]
  )
}))

###### The following code is adapted from McGuire et al. (2016) (www.pnas.org/cgi/doi/10.1073/pnas.1602817113) 
###### available at https://github.com/JennyMcGuire/ClimateConnectivity

# First, we will determine the Origin (warmer) & Dest (cooler) Cores for all adjacent core pairs
neighbors<-cbind(neighbor_table,Origin=0, Dest=0)

# Joining the temperature of each set of cores (Core1 & Core2)
colnames(neighbors)[1]<-"new_patch_id"
neighbors <- neighbors %>%
  inner_join(polygons_filtered %>% select(new_patch_id, current), by = "new_patch_id")
colnames(neighbors)[1]<-"new_patch_id_src"

colnames(neighbors)[2]<-"new_patch_id"
neighbors <- neighbors %>%
  inner_join(polygons_filtered %>% select(new_patch_id, current), by = "new_patch_id")
colnames(neighbors)[2]<-"new_patch_id_nbr"

# Remove geometry columns
neighbors <- neighbors[,-c(6,8)]

# Rename columns
colnames(neighbors)[colnames(neighbors)=="current.x"]<-"Mean1"
colnames(neighbors)[colnames(neighbors)=="new_patch_id_src"]<-"Cores1"
colnames(neighbors)[colnames(neighbors)=="current.y"]<-"Mean2"
colnames(neighbors)[colnames(neighbors)=="new_patch_id_nbr"]<-"Cores2"

n<-nrow(neighbors)

# Lists the hotter core # as the "Origin" & cooler core # as the "Dest"
for (i in 1:n)
{
  if (neighbors$Mean1[i]>neighbors$Mean2[i])
  {neighbors$Origin[i]<-neighbors$Cores1[i]
  neighbors$Dest[i]<-neighbors$Cores2[i]}
  
  else {neighbors$Origin[i]<-neighbors$Cores2[i]
  neighbors$Dest[i]<-neighbors$Cores1[i]}
}

colnames(neighbors)[5]<-"MeanO"
colnames(neighbors)[6]<-"MeanD"

# Identify neighbors
# First, transform polygons to dataframe to avoid issues with geometry
polygons_df <- st_drop_geometry(polygons_filtered)

connections <- setNames(
  lapply(polygons_df$new_patch_id, function(pid) {
    neighbors$Dest[neighbors$Origin == pid]
  }),
  polygons_df$new_patch_id
)

# Determine the ﬁnal, coolest destination patch
# Iterate over unique patch temperatures from cooler to hotter
uniquetemps <- sort(unique(polygons_df$current)) 
uniquetemps

# Set up output columns 
polygons_df$dest <- NA 
polygons_df$dest_temp <- NA 
polygons_df$inter_patch <- NA

# Copy the original dataframe
running <- polygons_df[, c("new_patch_id", "current")]

for (i in 1:length(uniquetemps)) {
  
  this_temp <- uniquetemps[i]
  inds <- which(running$current == this_temp)
  
  for (j in 1:length(inds)) {
    
    # get neighbors by patch ID
    this_pid  <- running$new_patch_id[inds[j]]
    dest_pids <- connections[[as.character(this_pid)]]
    
    if (length(dest_pids) > 0) { 
      
      # convert destination patch IDs to row indices
      dest_inds <- match(dest_pids, polygons_df$new_patch_id)
      
      # minimum destination temperature
      t <- min(running$current[dest_inds], na.rm = TRUE)
      min_ind <- dest_inds[which.min(running$current[dest_inds])]
      
      polygons_df$dest_temp[inds[j]] <- t 
      polygons_df$dest[inds[j]]      <- running$new_patch_id[min_ind]
      
      running$current[inds[j]]       <- t
      running$new_patch_id[inds[j]]  <- running$new_patch_id[min_ind]
      
      inter_patch <- polygons_df$inter_patch[min_ind]
      if (is.na(inter_patch)) {
        polygons_df$inter_patch[inds[j]] <- polygons_df$new_patch_id[min_ind]
      } else {
        polygons_df$inter_patch[inds[j]] <-
          paste(
            polygons_df$new_patch_id[min_ind],
            inter_patch[!is.na(inter_patch)],
            sep = ";"
          )
      }
      
    } else {
      # no neighbors -> stays where it started
      polygons_df$dest[inds[j]]      <- this_pid
      polygons_df$dest_temp[inds[j]] <- this_temp
    }
  }
}

# Assign destination future temperature
# SSP245
polygons_df$dest_ftemp_ssp245 <-
  vapply(1:nrow(polygons_df), function(x){
    dest <- polygons_df$dest[x]
    dest_ftemp <- polygons_df$ssp245[polygons_df$new_patch_id == dest] 
    return(dest_ftemp)
  }, FUN.VALUE = numeric(1))

# SSP585
polygons_df$dest_ftemp_ssp585 <-
  vapply(1:nrow(polygons_df), function(x){
    dest <- polygons_df$dest[x]
    dest_ftemp <- polygons_df$ssp585[polygons_df$new_patch_id == dest] 
    return(dest_ftemp)
  }, FUN.VALUE = numeric(1))

## Calculate climate connectivity by subtracting the destination future mean temperature 
# from the origin patch's current mean temperature
# Climate connectivity is calculated separately for SSP245 and SSP585
polygons_df$clim_conn_245 <- round(polygons_df$current - polygons_df$dest_ftemp_ssp245,1)
polygons_df$clim_conn_585 <- round(polygons_df$current - polygons_df$dest_ftemp_ssp585,1)

# Merge connectivity information to temperature patch polygons
polygons_connectivity <- polygons_filtered %>%
  left_join(polygons_df %>%
              select(new_patch_id, starts_with("clim_conn_")),
            by = "new_patch_id")

# Convert to rasters
conn_cols <- grep("^clim_conn_", names(polygons_df), value = TRUE)
conn_rasters <- lapply(conn_cols, function(col) {
  
  message("Processing raster: ", col)
  rasterize(polygons_connectivity, temp, field = col, fun = mean)
  
})

names(conn_rasters) <- conn_cols

connectivity_245 <- conn_rasters[[1]]
connectivity_585 <- conn_rasters[[2]]

# Plot
plot(connectivity_245)
plot(connectivity_585)

# Save as TIF
#writeRaster(connectivity_245, file.path(result.folder, "connectivity_245.tif"))
#writeRaster(connectivity_585, file.path(result.folder, "connectivity_585.tif"))

# Save CSV
#write.csv(polygons_df, file.path(result.folder, "connectivity_polygons.csv"))

# Save shapefile
#st_write(polygons_connectivity, file.path(result.folder, "connectivity_polygons.shp"))
