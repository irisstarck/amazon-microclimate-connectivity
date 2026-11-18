######################################################################################### 
#                                                                                       #
#                                 Patch role analysis                                   #
#                                                                                       #
######################################################################################### 
#
# This script assigns a role to each temperature patch based on the diversity of assumed species 
# movements experienced by the focal patch. The roles are as follows:
# 
# - Origin patches are patches for which temperature movements are systematically directed
# away from the focal patch (Class 1).
# - Neutral patches mean that the origin patch was not connected to any cooler patches and therefore 
# there is no possible movement to cooler conditions from this patch (Class 2). 
# - Intermediate patches are corridors on the way to the destination patch (Class 3). 
# - Destination patches mean that temperature movements are directed towards this patch.
# The destination role overrides any other role (Class 4).
#
# The script is divided into three parts:
# (1) Preparing data
# (2) Assign patch role
# (3) Patch roles in Protected Areas
# (4) Plot
#
# All data must be reprojected to EPSG:5880 with valid geometries.
#
# Inputs:
# - microclimate temperature patch polygons (output from script-01, "connectivity_polygons.shp")
# - microclimate connectivity CSV (output from script-01, "connectivity_polygons.csv")
# - protected area polygons: https://www.raisg.org/es/mapas/
# - microclimate raster for rasterization reference (the script uses microclimate connectivity raster)
# 
# Outputs:
# - microclimate temperature patch polygons with added role column
# - role raster
# - role proportion table in CSV
# - manuscript Figure 3, without country or PA borders
#
# All lines with saving ouputs are initially commented.
#
# -----------------------------------------------------------------------------
# Reproducibility note
# -----------------------------------------------------------------------------
#
# Script "Calculating microclimate connectivity in the Amazon" (01_calculate_climate_connectivity_amazon.R)
# is distributed with a small example dataset that allows users to test the workflow. 
# Because the example dataset represents only a small subset of the Amazon biome, the resulting
# connectivity rasters and polygons are not suitable for reproducing the patch role analyses
# presented in this script.
#
# Therefore, this script uses the Amazon-wide connectivity rasters and polygons generated from
# the full study datasets. These rasters are provided in 
# https://doi.org/10.23729/fd-2d2803a3-8804-3fea-a31f-c2f5e74e10c2 and should be 
# used as inputs for all analyses below.
#
######################################################################################### 
#                                                                                       #
#                                     1. Preparations                                   #
#                                                                                       #
######################################################################################### 

# Load libraries
library(data.table)
library(sf)
library(terra)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Define data folders
data.dir <- "data"
result.folder <- "result.folder"

# Read input files
polygons_connectivity <- st_read(file.path(result.folder, "connectivity_polygons.shp"))
polygons_df <- fread(file.path(result.folder, "connectivity_polygons.csv"))
pa <- st_read(file.path(data.dir, "file"))
connectivity_245 <-  rast(file.path(result.folder, "connectivity_245.tif")) # can also be original microclimate map

# Ensure same CRS
pa <- project(pa, crs(connectivity_245))
pa <- makeValid(pa)

######################################################################################### 
#                                                                                       #
#                                  2. Assign patch role                                 #
#                                                                                       #
######################################################################################### 

## Remove intermediate patch IDs that equal destination patch ID

# Create function
clean_inter_patch <- function(df) {
  df %>%
    mutate(
      inter_patch = ifelse(
        is.na(inter_patch),
        NA_character_,
        sapply(seq_along(inter_patch), function(i) {
          parts <- unlist(strsplit(inter_patch[i], ";"))
          parts <- parts[parts != as.character(dest[i])]
          paste(parts, collapse = ";")
        })
      )
    )
}

# Apply
df_cleaned <- clean_inter_patch(polygons_df)

## Find the role of each patch (origin/neutral/intermediate/destination)

# Create function for assigning role
assign_flow_role <- function(df) {
  
  # Split inter_patch string into separate rows
  long_paths <- df %>%
    mutate(inter_patch = ifelse(is.na(inter_patch), "", inter_patch)) %>%
    separate_rows(inter_patch, sep = ";") %>%
    mutate(
      inter_patch = suppressWarnings(as.integer(inter_patch))
    )
  
  # Unique lists of patch IDs by role
  origins      <- unique(df$new_patch_id)
  destinations <- unique(long_paths$dest[long_paths$dest != long_paths$new_patch_id])
  intermediates <- unique(long_paths$inter_patch[!is.na(long_paths$inter_patch)])
  neutrals <- unique(long_paths$new_patch_id[long_paths$new_patch_id == long_paths$dest])
  
  # Role assignment
  # First, all will be assigned as origins
  roles <- tibble(
    new_patch_id = unique(df$new_patch_id),
    role = "origin"
  )
  
  roles$role[roles$new_patch_id %in% intermediates &
               roles$role != "destination"] <- "intermediate"
  
  roles$role[roles$new_patch_id %in% neutrals] <- "neutral"
  
  roles$role[roles$new_patch_id %in% destinations] <- "destination"
  
  # Count how many times each patch "acts" as destinations or intermediate
  
  # Destination count:
  dest_count <- df %>%
    count(dest, name="n_dest_acts") %>%
    rename(new_patch_id = dest)
  
  # Intermediate count:
  inter_count <- long_paths %>%
    filter(!is.na(inter_patch)) %>%
    count(inter_patch, name="n_intermediate_acts") %>%
    rename(new_patch_id = inter_patch)
  
  # merge counts (patches not appearing in a role get NA -> 0)
  roles <- roles %>%
    left_join(dest_count, by="new_patch_id") %>%
    left_join(inter_count, by="new_patch_id") %>%
    mutate(
      n_dest_acts = replace_na(n_dest_acts, 0),
      n_intermediate_acts = replace_na(n_intermediate_acts, 0)
    )
  print(sum(is.na(roles$role)))
  return(roles)
}

# Apply on temperature patch csv
patch_roles <- assign_flow_role(df_cleaned)

## Join role information and patch shapefile

# Rename columns to match
if ("nw_ptc_" %in% names(polygons_connectivity)) {
  polygons_connectivity <- dplyr::rename(polygons_connectivity, new_patch_id = nw_ptc_)
}

# Add role columns to shapefile
patch_with_roles <- polygons_connectivity %>%
  left_join(patch_roles, by = "new_patch_id")

# Add numerical role class for rasterization
patch_with_roles$role_class <- NA

patch_with_roles <- patch_with_roles %>%
  mutate(role_class = case_when(
    role == "origin"       ~ 1,
    role == "neutral"      ~ 2,
    role == "intermediate" ~ 3,
    role == "destination"  ~ 4
    
  ))

# Rasterize
r.roles <- rasterize(vect(patch_with_roles), connectivity_245, field = "role_class", fun = max)

# Save raster
#writeRaster(r.roles, file.path(result.folder, "r.roles.tif"))

######################################################################################### 
#                                                                                       #
#                             3. Patch roles in Protected Areas                         #
#                                                                                       #
######################################################################################### 

# Mask role raster with protected areas
roles_pa <- mask(r.roles, pa)

# Mask areas outside protected areas
roles_opa <- mask(r.roles, pa, inverse = TRUE)

## Calculate areas inside and outside of PAs

# Function for calculating cell size in km2:
cell_area_km2 <- function(r) {
  terra::cellSize(r, unit = "km")
}

a <- cell_area_km2(r.roles)

classes <- 1:4

role_areas <- do.call(rbind, lapply(classes, function(k) {
  data.frame(
    class       = k,
    area_inside  = global(a * (roles_pa  == k), "sum", na.rm = TRUE)[1, 1],
    area_outside = global(a * (roles_opa == k), "sum", na.rm = TRUE)[1, 1]
  )
}))

role_areas$prop_inside <- role_areas$area_inside /
  (role_areas$area_inside + role_areas$area_outside)

role_areas$prop_outside <- role_areas$area_outside /
  (role_areas$area_inside + role_areas$area_outside)

# Add role information as text
role_areas$role <- "origin"
role_areas[role_areas$class == 2, "role"] <- "neutral"
role_areas[role_areas$class == 3, "role"] <- "intermediate"
role_areas[role_areas$class == 4, "role"] <- "destination"

#Save CSV
#write.csv(role_areas, file.path(result.folder, "role_areas.csv"))

######################################################################################### 
#                                                                                       #
#                                       4. Plot                                         #
#                                                                                       #
######################################################################################### 

## Manuscript Figure 3, without country or PA borders

## Visualization parameters and preparations

# Legend labels
legend_labels <- c(
  "origin"       = "Origin",
  "destination"  = "Destination",
  "intermediate" = "Intermediate",
  "neutral"      = "Neutral"
)

# Create common fill
common_fill <- scale_fill_manual(
  values = c(
    "origin"       = "#E69F00",
    "neutral"      = "white",
    "intermediate" = "#009E73",
    "destination"  = "#56B4E9"
  ), 
  labels = legend_labels,
  guide = guide_legend(
    override.aes = list(
      color = c("black", "black", "black", "black"),  # adds border to all legend keys
      linewidth  = 0.5
    )
  ), drop = FALSE
)

# Role levels
role_levels <- c("origin", "neutral", "intermediate", "destination")

# Define levels for input data
patch_with_roles <- patch_with_roles %>%
  mutate(
    role = factor(role, levels = role_levels))

## Map of patch roles

# Remove the smallest PAs to make the map nicer
# Convert PAs to sf
pa <- st_as_sf(pa)
# Filter
pa_filtered <- pa %>%
  mutate(area_m2 = st_area(geometry)) %>%        # compute area
  mutate(area_km2 = as.numeric(area_m2) / 1e6) %>%   # convert to km²
  filter(area_km2 >= 100)

# Plot
role_plot <- ggplot(patch_with_roles) + #patch_with_roles
  geom_sf(aes(fill = role), color = NA) +
  geom_sf(data = pa_filtered, fill = NA, color = "black", linewidth = 0.1)+
  labs(
    fill = "",
    tag = "a")+
  common_fill +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "lightgray", color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    plot.tag.position = c(0.15, 0.95)
  ) 

##  Barplot of area share

# Create long format data
role_long <- role_areas %>% 
  select(role, prop_inside, prop_outside) %>%
  pivot_longer(
    cols = c(prop_inside, prop_outside),
    names_to = "location",
    values_to = "prop"
  ) %>%
  mutate(
    role = factor(role, levels = role_levels),
    location = factor(recode(location,
                             prop_inside  = "PA",
                             prop_outside = "OPA"),
                      levels = c("OPA", "PA"))   # PA stacks first (left after coord_flip)
  )

# Barplot
role_barplot <- ggplot(role_long,
                             aes(x = fct_rev(role), y = prop, fill = role, alpha = location)) +
  geom_col(color = "black", width = 0.5) +
  geom_text(
    aes(x = fct_rev(role), y = prop, group = location, label = location),
    position = position_stack(vjust = 0.5),
    size = 3.5, color = "black", inherit.aes = FALSE
  ) +
  coord_flip() +
  guides(fill = "none") +  
  scale_y_continuous(limits = c(0, 1), labels = scales::percent, expand = c(0, 0)) +
  scale_alpha_manual(values = c("PA" = 1, "OPA" = 0.35), guide = "none") +
  common_fill +
  labs(x = "", y = "Proportion of area", tag = "b") +
  theme_minimal(base_size = 14) +
  theme(
    plot.tag.position = c(0.02, 0.95),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_blank(),
    plot.margin = margin(t = 5, r = 15, b = 5, l = 5)
  )

# Combine
role_plot_combined <- role_plot + role_barplot +
  plot_layout(widths = c(1, 0.5), guides = "collect") & theme(legend.position = "bottom")

# Save figure
#ggsave(
  #filename = "Figure3.png",
  #plot = role_plot_combined,
  #path = result.folder,
  #width = 20,
  #height = 12,
  #units = "cm",
  #dpi = 300,
  #bg = "white"
#)
