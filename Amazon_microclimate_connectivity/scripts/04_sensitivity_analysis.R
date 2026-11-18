######################################################################################### 
#                                                                                       #
#                                   Sensitivity analysis                                #
#                                                                                       #
######################################################################################### 
#
# This script evaluates the sensitivity of microclimate connectivity to different
# warming magnitudes between present conditions and the SSP585 scenario.
#
# The workflow:
# (1) Prepare data
# (2) Simulate intermediate warming scenarios
# (3) Calculate climate connectivity for each warming scenario
# (4) Basic statistics for each warming scenario
# (5) Plot
#
# Inputs:
# - Present-day microclimate temperature raster
# - Future microclimate temperature raster (SSP585)
# - Present-day macroclimate temperature raster
# - Future macroclimate temperature raster (SSP585)
# - Temperature patch polygons generated in Script 01
#
# Outputs:
# - Scenario table CSV
# - Summary statistics table CSV
# - Area proportion table CSV
# - Manuscript Figure 2
#
# All rasters and vector layers must:
# - Be projected to EPSG:5880
# - Cover the same spatial extent
# - Have matching resolution and alignment
#
# All lines that save outputs are initially commented.
#
# -----------------------------------------------------------------------------
# Reproducibility note
# -----------------------------------------------------------------------------
#
# Script "Calculating microclimate connectivity in the Amazon" (01_calculate_climate_connectivity_amazon.R)
# is distributed with a small example dataset that allows users to test the workflow. 
# Because the example dataset represents only a small subset of the Amazon biome, the resulting
# connectivity rasters and polygons are not suitable for reproducing the sensitivity analyses
# presented in this script.
#
# Therefore, this script uses the Amazon-wide connectivity rasters generated from
# the full study datasets. These rasters are provided with the repository (or are
# available through the repository archive) and should be used as inputs for all
# analyses below.

######################################################################################### 
#                                                                                       #
#                                     1. Preparations                                   #
#                                                                                       #
######################################################################################### 

# Load libraries
library(terra)
library(sf)
library(dplyr)

# Define data folders
data.dir <- "data"
result.folder <- "result.folder"

# Load data
micro_t <- rast(file.path(data.dir, "micro_t_filtered.tif"))
micro_585 <- rast(file.path(data.dir, "micro_585_filtered.tif"))
macro_now <- rast(file.path(data.dir, "macroclimate_file")) 
macro_585 <- rast(file.path(data.dir, "macroclimate_585_file")) 
polygons_filtered <- st_read(file.path(result.folder, "polygons_filtered.shp"))

######################################################################################### 
#                                                                                       #
#                             2. Simulate warming scenarios                             #
#                                                                                       #
######################################################################################### 

# Calculate difference between now and future
dif_585 <- micro_585 - micro_t
dif_585_macro <- macro_585 - macro_now

# sequence of warming fractions
alpha_seq <- seq(0.1, 1, by = 0.1)

# Generate rasters from present to temperatures corresponding to SSP585 with 10% increases
warming_scenarios <- lapply(alpha_seq, function(a) {
  micro_t + dif_585 * a
})

# Same for macroclimate
warming_scenarios_macro <- lapply(alpha_seq, function(a) {
  macro_now + dif_585_macro * a
})

# Name scenarios "warm_1", "warm_2", etc
names(warming_scenarios) <- paste0("warm_", 1:length(alpha_seq))
names(warming_scenarios_macro) <- paste0("warm_", 1:length(alpha_seq))

# Convert lists to rasters
warming_scenarios <- rast(warming_scenarios)
warming_scenarios_macro <- rast(warming_scenarios_macro)

# Get mean microclimate temperature increase for each warming scenario
mean_warming_scenarios <- warming_scenarios - micro_t
mean_diffs_df <- as.data.frame(global(mean_warming_scenarios, fun = "mean", na.rm = TRUE))
# Give scenario names matching to the raster names
mean_diffs_df$scenario <- names(warming_scenarios)
# Secure correct order
mean_diffs_df <- mean_diffs_df[
  order(as.numeric(gsub("warm_", "", mean_diffs_df$scenario))),
]

# Get mean macroclimate temperature increase for each warming scenario
mean_warming_scenarios_macro <- warming_scenarios_macro - macro_now
mean_diffs_df_macro <- as.data.frame(global(mean_warming_scenarios_macro, fun = "mean", na.rm = TRUE))
# Give scenario names matching to the raster names
mean_diffs_df_macro$scenario <- names(warming_scenarios_macro)
names(mean_diffs_df_macro) <- c("mean_macro", "scenario")
# Secure correct order
mean_diffs_df_macro <- mean_diffs_df_macro[
  order(as.numeric(gsub("warm_", "", mean_diffs_df_macro$scenario))),
]

######################################################################################### 
#                                                                                       #
#               3. Calculate microclimate connectivity for each scenario                #
#                                                                                       #
######################################################################################### 

# Extract mean temperature values of present-day and warming scenario mean temperatures
# per polygon
mean_values <- terra::extract(warming_scenarios, vect(polygons_filtered), fun = mean, na.rm = TRUE)

# Round all warm_* columns to 1 decimal
mean_values[, -1] <- round(mean_values[, -1], 1)

# Calculate connectivity of all steps

# Add mean temperatures as new columns in polygons
polygons_filtered <- cbind(polygons_filtered, mean_values[,-1])

## Calculate climate connectivity similarly as in "01_calculate_climate_connectivity_amazon.R"

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

# First we will determine the Origin (warmer) & Dest (cooler) Cores for all adjacent core pairs
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

# Identify warming scenario columns
warm_cols <- grep("^warm_", names(polygons_df), value = TRUE)

# For each warming scenario, compute destination future temperature
for (w in warm_cols) {
  
  # name the new dest-temp column
  dest_col <- paste0("dest_", w)
  
  # fill column
  polygons_df[[dest_col]] <- vapply(
    1:nrow(polygons_df),
    function(i){
      dest <- polygons_df$dest[i]
      dest_ftemp <- polygons_df[[w]][polygons_df$new_patch_id == dest] 
      return(dest_ftemp)
    },
    FUN.VALUE = numeric(1)
  )
}

# Now compute climate connectivity for all warming scenarios
for (w in warm_cols) {
  
  dest_col <- paste0("dest_", w)
  conn_col <- paste0("clim_conn_", sub("warm_", "", w))
  
  polygons_df[[conn_col]] <- round(polygons_df$current - polygons_df[[dest_col]], 1)
}

# Save result as CSV
#write.csv(polygons_df, file.path(result.folder, "micro_connectivity_sensitivity_test.csv"))

# Merge connectivity information to forest polygons
polygons_connectivity <- polygons_filtered %>%
  left_join(polygons_df %>% 
              select(new_patch_id, starts_with("clim_conn_")),
            by = "new_patch_id"
  )

## Convert to rasters

conn_cols <- grep("^clim_conn_", names(polygons_connectivity), value = TRUE) # original

# Create reference
template <- micro_t

#Rasterize
conn_rasters <- lapply(conn_cols, function(col) {
  print(col)
  r <- rasterize(polygons_connectivity, template, field = col, fun = mean)
  names(r) <- col
  #outfile <- file.path(result.folder, paste0(col.tif"))
  #terra::writeRaster(r, outfile, overwrite = TRUE)
  #message("Saved: ", outfile)
  return(r)
})

names(conn_rasters) <- conn_cols

# Check that column names match
if ("nw_ptc_" %in% names(polygons_connectivity)) {
  polygons_connectivity <- dplyr::rename(polygons_connectivity, new_patch_id = nw_ptc_)
}

######################################################################################### 
#                                                                                       #
#                     4. Basic statistics for each warming scenario                     #
#                                                                                       #
######################################################################################### 

# Raster statistics
summary_table_sensitivity <- lapply(names(conn_rasters), function(nm) {
  
  r <- conn_rasters[[nm]]
  vals <- values(r, na.rm = TRUE)
  
  data.frame(
    scenario = nm,
    mean = mean(vals),
    sd   = sd(vals),
    min  = min(vals),
    max  = max(vals)
  )
})

summary_table_sensitivity <- do.call(rbind, summary_table_sensitivity)

# Add mean temperature increase
summary_table_sensitivity$temp_inc <- round(mean_diffs_df$mean,2)
summary_table_sensitivity$temp_inc_macro <- round(mean_diffs_df_macro$mean_macro,2)

# Add proportional increase column
summary_table_sensitivity$increase <- alpha_seq

# Save
#write.csv(summary_table_sensitivity, file.path(result.folder,  "summary_statistics_sensitivity_test.csv"))

## Proportions of positive/negative/zero (connected/unconnected area), area-based
# Function for calculating cell size in km2:
cell_area_km2 <- function(r) {
  terra::cellSize(r, unit = "km")
}

area_prop_table_sensitivity <- lapply(names(conn_rasters), function(nm) {
  
  r <- conn_rasters[[nm]]
  a <- cellSize(r, unit = "km")
  
  pos  <- a * (r > 0)
  neg  <- a * (r < 0)
  zero <- a * (r == 0)
  
  area_pos  <- global(pos,  "sum", na.rm = TRUE)[1,1]
  area_neg  <- global(neg,  "sum", na.rm = TRUE)[1,1]
  area_zero <- global(zero, "sum", na.rm = TRUE)[1,1]
  
  area_total <- area_pos + area_neg + area_zero
  
  data.frame(
    scenario = nm,
    area_pos  = area_pos,
    area_neg  = area_neg,
    area_zero = area_zero,
    prop_pos  = area_pos  / area_total,
    prop_neg  = area_neg  / area_total,
    prop_zero = area_zero / area_total
  )
}) 

area_prop_table_sensitivity <- do.call(rbind, area_prop_table_sensitivity)

# Merge mean temperature increase values
area_prop_table_sensitivity$temp_inc <- round(mean_diffs_df$mean,2)
area_prop_table_sensitivity$temp_inc_macro <- round(mean_diffs_df_macro$mean_macro,2)

#Add proportional increase column
area_prop_table_sensitivity$increase <- alpha_seq

# Save
#write.csv(area_prop_table_sensitivity, file.path(result.folder,  "area_prop_table_sensitivity_corrected.csv"))

######################################################################################### 
#                                                                                       #
#                                        5. Plot                                        #
#                                                                                       #
######################################################################################### 

## Manuscript Figure 2

## Visualization parameters and preparations

#Create plot data
plot_data <- area_prop_table_sensitivity

# Prepare matrix for stacking
mat <- as.matrix(t(plot_data[, c("prop_pos", "prop_zero", "prop_neg")] * 100))
rownames(mat) <- c("Positive", "Zero", "Negative")

# Colors
bar_colors <- c("royalblue", "lightgrey", "coral")
# Labels
bar_labels <- c("Positive", "Zero", "Negative")

# Uncomment if saving
#png(
  #filename = "connectivity_scenarios.png",
  #path = result.folder,
  #width = 190,
  #height = 150,
  #units = "mm",
  #res = 300, 
  #type = "cairo",
  #bg = "white"
#)

## Plot

par(
  mar = c(8, 4, 3, 5),
  mgp = c(2.1, 0.15, 0),
  las = 1
)

# Compute bar positions
bp <- barplot(
  mat,
  plot = FALSE
)

# Empty plotting region
plot(
  NA, xlim = range(bp) + c(-0.5, 0.5),
  ylim = c(0, 100),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  bty  = "n"
)

# Background gridlines
abline(
  h = seq(0, 100, by = 10),
  col = "#d6d6d6",
  lwd = 1
)

# Bars
barplot(
  mat,
  add = TRUE,
  beside = FALSE,
  col = bar_colors,
  border = NA,
  ylim = c(0, 100),
  xaxt = "n",
  yaxt = "n"
)

x_center <- mean(par("usr")[1:2])
y_pos <- -28

legend(
  x = x_center + 2,
  y = y_pos,
  legend = bar_labels,
  fill = bar_colors,
  border = NA,
  bty = "n",
  horiz = TRUE,
  xpd = TRUE,
  xjust = 0.5,
  yjust = 0.5
)

text(
  x = x_center - 2,
  y = y_pos,
  labels = "Climate connectivity:",
  adj = c(1, 0.5),
  xpd = TRUE
)

# Bottom axis
axis(
  1,
  at = bp,
  labels = plot_data$increase,
  col = NA,
  col.axis = "#4d4d4d",
  tick = FALSE
)

# Left axis
axis(
  2,
  col = NA,
  col.axis = "#4d4d4d",
  tick = FALSE
)

# Axis labels
mtext(
  "Proportion of temperature increase from present to SSP585 in 2100",
  side = 1,
  line = 1.4
)

mtext(
  "Proportion of area (%)",
  side = 2,
  line = 1.4,
  las = 0
)

# Temperature legend
temp_y_pos <- y_pos - 8

legend(
  x = x_center + 2,
  y = temp_y_pos,
  legend = c("Macroclimate", "Microclimate"),
  lty = c(1, 2),
  lwd = 1.5,
  pch = c(16, 1),
  bty = "n",
  horiz = TRUE,
  xpd = TRUE,
  xjust = 0.5,
  yjust = 0.5
)

text(
  x = x_center - 2,
  y = temp_y_pos,
  labels = "Temperature increase:",
  adj = c(1, 0.5),
  xpd = TRUE
)

# Temperature lines
temp_min <- floor(
  min(c(
    plot_data$temp_inc,
    plot_data$temp_inc_macro
  ), na.rm = TRUE)
)

temp_max <- ceiling(
  max(c(
    plot_data$temp_inc,
    plot_data$temp_inc_macro
  ), na.rm = TRUE)
)

par(new = TRUE)

plot(
  bp,
  plot_data$temp_inc_macro,
  type = "o",
  pch = 16,
  lwd = 1.5,
  lty = 1,
  col = "black",
  axes = FALSE,
  xlab = "",
  ylab = "",
  xlim = range(bp) + c(-0.5, 0.5),
  ylim = c(temp_min, temp_max)
)

lines(
  bp,
  plot_data$temp_inc,
  type = "o",
  pch = 1,
  lwd = 1.5,
  lty = 2,
  col = "black"
)

# Right-side temperature axis
axis(
  4,
  at = pretty(c(temp_min, temp_max)),
  col = NA,
  col.axis = "#4d4d4d",
  tick = FALSE
)

mtext(
  "Mean temperature increase by 2100 (°C)",
  side = 4,
  line = 2,
  las = 3
)

# Top axis title
mtext(
  "Sensitivity analysis of microclimate connectivity",
  side = 3, line = 1, col = "black", cex = 1.2
)

#dev.off()


