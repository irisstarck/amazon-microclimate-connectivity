######################################################################################### 
#                                                                                       #
#                 Microclimate connectivity statistics in the Amazon                    #
#                                                                                       #
######################################################################################### 
#
# This script includes the statistical analyses for microclimate connectivity.
# 
# The exact same workflow applies for macroclimate connectivity. Input raster data must be 
# changed accordingly and reprojected to match the microclimate rasters.
#
# The script is divided into four parts:
# (1) Preparing data
# (2) Calculating basic raster statistics for the whole Amazon, inside (PA) and outside (OPA) Protected Areas
# (3) Estimating effect of PA sizeon climate connectivity
# (4) Calculating statistical difference and effect size between PAs and OPAs
# including calculating spatial autocorrelation with Moran's I
#
# All data must be in EPSG:5880
#
# Inputs:
# - microclimate connectivity raster for SSP245
#   MICC_ssp245 in https://doi.org/10.23729/fd-2d2803a3-8804-3fea-a31f-c2f5e74e10c2
#
# - microclimate connectivity raster for SSP585
#   MICC_ssp585 in https://doi.org/10.23729/fd-2d2803a3-8804-3fea-a31f-c2f5e74e10c2
#
# - protected area polygons: https://www.raisg.org/es/mapas/
#
# Outputs:
# - summary statistics tables for the whole Amazon
# - summary statistics tables for PAs
# - summary statistics tables for OPAs
# - tables of positive/negative/zero climate connectivity for the whole Amazon
# - tables of positive/negative/zero climate connectivity for PAs
# - tables of positive/negative/zero climate connectivity for OPAs
# - table of Moran's I within spatial blocks
# - table of t-test and paired Wilcox test results using spatial blocks
# - table of Cohen's d
#
# -----------------------------------------------------------------------------
# Reproducibility note
# -----------------------------------------------------------------------------
#
# Script "Calculating microclimate connectivity in the Amazon" (01_calculate_climate_connectivity_amazon.R)
# is distributed with a small example dataset that allows users to test the workflow. 
# Because the example dataset represents only a small subset of the Amazon biome, the resulting
# connectivity rasters are not suitable for reproducing the statistical analyses
# presented in this script.
#
# Therefore, this script uses the Amazon-wide connectivity rasters generated from
# the full study datasets. These rasters are provided in 
# https://doi.org/10.23729/fd-2d2803a3-8804-3fea-a31f-c2f5e74e10c2 and should be 
# used as inputs for all analyses below.

######################################################################################### 
#                                                                                       #
#                                     1. Preparations                                   #
#                                                                                       #
######################################################################################### 

# Load libraries
library(terra)
library(tidyr)
library(mgcv)
library(dplyr)
library(effectsize)
library(ape)  

# Define data folders
data.dir <- "data"
result.folder <- "result.folder"

# Load data:
connectivity_245 <- rast(file.path(result.folder, "connectivity_245.tif"))
connectivity_585 <- rast(file.path(result.folder, "connectivity_585.tif"))

pa <- vect(file.path(data.dir, "file"))

# Ensure same CRS
pa <- project(pa, crs(connectivity_245))
pa <- makeValid(pa)

# Make into list
connectivity_r_list <- list(
  connectivity_245 = connectivity_245,
  connectivity_585 = connectivity_585
)

# Mask connectivity rasters with protected area polygons
connectivity_245_pa <- mask(connectivity_245, pa)
connectivity_585_pa <- mask(connectivity_585, pa)

#Make into list
r_pa_list <- list(
  connectivity_245_pa = connectivity_245_pa,
  connectivity_585_pa = connectivity_585_pa
)

# Mask areas outside protected areas
connectivity_245_opa <- mask(connectivity_245, pa, inverse = TRUE)
connectivity_585_opa <- mask(connectivity_585, pa, inverse = TRUE)

#Make into list
r_opa_list <- list(
  connectivity_245_opa = connectivity_245_opa,
  connectivity_585_opa = connectivity_585_opa
)

######################################################################################### 
#                                                                                       #
#                                    2. Basic statistics                                #
#                                                                                       #
######################################################################################### 

## Statistics for the whole Amazon:

# Raster statistics
summary_table <- lapply(names(connectivity_r_list), function(nm) {
  
  r <- connectivity_r_list[[nm]]
  vals <- values(r, na.rm = TRUE)
  
  data.frame(
    dataset  = sub("_(245|585)", "", nm),
    scenario = ifelse(grepl("245", nm), "SSP245", "SSP585"),
    mean = mean(vals),
    sd   = sd(vals),
    min  = min(vals),
    max  = max(vals)
  )
})

summary_table <- do.call(rbind, summary_table)

#Write csv
#write.csv(summary_table, file.path(result.folder, "summary_statistics.csv"))

# Proportions of positive/negative/zero (connected/unconnected area), area-based
# Function for calculating cell size in km2:
cell_area_km2 <- function(r) {
  terra::cellSize(r, unit = "km")
}

props_weighted <- lapply(names(connectivity_r_list), function(nm) {
  
  r <- connectivity_r_list[[nm]]
  
  # cell area
  a <- cell_area_km2(r)
  
  # area by sign
  pos  <- a * (r > 0)
  neg  <- a * (r < 0)
  zero <- a * (r == 0)
  
  # sum areas
  area_pos  <- terra::global(pos,  "sum", na.rm = TRUE)[1,1]
  area_neg  <- terra::global(neg,  "sum", na.rm = TRUE)[1,1]
  area_zero <- terra::global(zero, "sum", na.rm = TRUE)[1,1]
  
  area_total <- area_pos + area_neg + area_zero
  
  data.frame(
    dataset = sub("_(245|585)", "", nm),
    scenario = ifelse(grepl("245", nm), "SSP245", "SSP585"),
    area_pos  = area_pos,
    area_neg  = area_neg,
    area_zero = area_zero,
    prop_pos  = area_pos  / area_total,
    prop_neg  = area_neg  / area_total,
    prop_zero = area_zero / area_total
  )
}) 

props_weighted <- do.call(rbind, props_weighted)

#Write csv
#write.csv(props_weighted, file.path(result.folder, "props_weighted.csv"))

## Statistics for Protected Areas:

# Raster statistics
pa_summary_table <- lapply(names(r_pa_list), function(nm) {
  
  r <- r_pa_list[[nm]]
  vals <- values(r, na.rm = TRUE)
  
  data.frame(
    dataset  = sub("_(245|585)", "", nm),
    scenario = ifelse(grepl("245", nm), "SSP245", "SSP585"),
    mean = mean(vals),
    sd   = sd(vals),
    min  = min(vals),
    max  = max(vals)
  )
})

pa_summary_table <- do.call(rbind, pa_summary_table)

#Write csv
#write.csv(pa_summary_table, file.path(result.folder, "summary_statistics_pa.csv"))

# Proportions of positive/negative/zero (connected/unconnected area)
props_weighted_pa <- lapply(names(r_pa_list), function(nm) {
  
  r <- r_pa_list[[nm]]
  
  # cell area
  a <- cell_area_km2(r)
  
  # area by sign
  pos  <- a * (r > 0)
  neg  <- a * (r < 0)
  zero <- a * (r == 0)
  
  # sum areas
  area_pos  <- terra::global(pos,  "sum", na.rm = TRUE)[1,1]
  area_neg  <- terra::global(neg,  "sum", na.rm = TRUE)[1,1]
  area_zero <- terra::global(zero, "sum", na.rm = TRUE)[1,1]
  
  area_total <- area_pos + area_neg + area_zero
  
  data.frame(
    dataset = sub("_(245|585)", "", nm),
    scenario = ifelse(grepl("245", nm), "SSP245", "SSP585"),
    area_pos  = area_pos,
    area_neg  = area_neg,
    area_zero = area_zero,
    prop_pos  = area_pos  / area_total,
    prop_neg  = area_neg  / area_total,
    prop_zero = area_zero / area_total
  )
}) 

props_weighted_pa <- do.call(rbind, props_weighted_pa)

#Write csv
#write.csv(props_weighted_pa, file.path(result.folder, "props_weighted_pa.csv"))

## Statistics Outside of Protected Areas:

#Summary statistics
opa_summary_table <- lapply(names(r_opa_list), function(nm) {
  
  r <- r_opa_list[[nm]]
  vals <- values(r, na.rm = TRUE)
  
  data.frame(
    dataset  = sub("_(245|585)", "", nm),
    scenario = ifelse(grepl("245", nm), "SSP245", "SSP585"),
    mean = mean(vals),
    sd   = sd(vals),
    min  = min(vals),
    max  = max(vals)
  )
})

opa_summary_table <- do.call(rbind, opa_summary_table)

#Write csv
#write.csv(opa_summary_table, file.path(result.folder, "summary_statistics_opa.csv"))

#Proportions of connected and unconnected area
props_weighted_opa <- lapply(names(r_opa_list), function(nm) {
  
  r <- r_opa_list[[nm]]
  
  a <- cellSize(r, unit = "km")
  
  pos  <- a * (r > 0)
  neg  <- a * (r < 0)
  zero <- a * (r == 0)
  
  area_pos  <- global(pos,  "sum", na.rm = TRUE)[1,1]
  area_neg  <- global(neg,  "sum", na.rm = TRUE)[1,1]
  area_zero <- global(zero, "sum", na.rm = TRUE)[1,1]
  
  area_total <- area_pos + area_neg + area_zero
  
  data.frame(
    dataset = sub("_(245|585)", "", nm),
    scenario = ifelse(grepl("245", nm), "SSP245", "SSP585"),
    area_pos  = area_pos,
    area_neg  = area_neg,
    area_zero = area_zero,
    prop_pos  = area_pos  / area_total,
    prop_neg  = area_neg  / area_total,
    prop_zero = area_zero / area_total
  )
}) 

props_weighted_opa <- do.call(rbind, props_weighted_opa)

#Write csv
#write.csv(props_weighted_opa, file.path(result.folder, "props_weighted_opa.csv"))

######################################################################################### 
#                                                                                       #
#                 3. Effect of Protected Area size on climate connectivity              #
#                                                                                       #
######################################################################################### 

# Calculate area of each protected area (PA)
pa$area_km2 <- terra::expanse(pa, unit = "km")

## Extract connectivity value for each PA
# SSP245
pa_245 <- terra::extract(
  connectivity_245,
  pa,
  fun = mean,
  na.rm = TRUE
)

# SSP585
pa_585 <- terra::extract(
  connectivity_585,
  pa,
  fun = mean,
  na.rm = TRUE
)

# Add connectivity information to the vector
pa$connectivity_245 <- pa_245[, 2]
pa$connectivity_585 <- pa_585[, 2]

## Test the effect of PA size on climate connectivity using a Generalized Additive Model (GAM) and Spearman's correlation

# Check data
hist(pa$connectivity_245, breaks = 50)

plot(log10(pa$area_km2), pa$connectivity_245,
     pch = 16, col = rgb(0, 0, 0, 0.2))

plot(log10(pa$area_km2), pa$connectivity_585,
     pch = 16, col = rgb(0, 0, 0, 0.2))

# SSP245
m1 <- gam(
  connectivity_245 ~ s(log10(area_km2)),
  data = as.data.frame(pa),
  method = "REML"
)

# SSP585
m2 <- gam(
  connectivity_585 ~ s(log10(area_km2)),
  data = as.data.frame(pa),
  method = "REML"
)

## Model summaries
summary(m1)
summary(m2)

# GAM check
gam.check(m1)
gam.check(m2)

# Sensitivity test 1: increase k to 20

# SSP245
m1_check <- gam(
  connectivity_245 ~ s(log10(area_km2), k = 20),
  data = as.data.frame(pa),
  method = "REML"
)

# SSP585
m2_check <- gam(
  connectivity_585 ~ s(log10(area_km2), k = 20),
  data = as.data.frame(pa),
  method = "REML"
)

gam.check(m1_check)
gam.check(m2_check)

summary(m1_check)
summary(m2_check)

# Sensitivity test 1: increase k to 40
# SSP245
m1_check2 <- gam(
  connectivity_245 ~ s(log10(area_km2), k = 40),
  data = as.data.frame(pa),
  method = "REML"
)

# SSP585
m2_check2 <- gam(
  connectivity_585 ~ s(log10(area_km2), k = 40),
  data = as.data.frame(pa),
  method = "REML"
)

gam.check(m1_check2)
gam.check(m2_check2)

summary(m1_check2)
summary(m2_check2)

## Spearman's correlation
# SSP245
cor.test(
  pa$area_km2,
  pa$connectivity_245,
  method = "spearman",
  exact = FALSE
)

# SSP585
cor.test(
  pa$area_km2,
  pa$connectivity_585,
  method = "spearman",
  exact = FALSE
)
######################################################################################### 
#                                                                                       #
#                       4. Statistical difference between OPA and PA                    #
#                                                                                       #
#########################################################################################  

## Sampling function

# sample number
n <- 5000

# Create function
sample_vals <- function(r, n, seed = 1) {
  set.seed(seed)
  s <- spatSample(r, size = n, method = "random", na.rm = TRUE, xy = TRUE)
  names(s)[3] <- "conn"
  s
}

# Create scenario variable
scen <- c(SSP245 = 1, SSP585 = 2)

## Test spatial autocorrelation using spatial blocks

# Create function block_diffs():
# k = scenario (SSP245 or SSP585)
# block_size = side length of the square blocks, in metres
# n = number of random cells drawn inside PAs and outside PAs (separately)
# min_n = minimum number of sampled cells per group for a block to be used

# Uses the same sampling and blocking as block_test(), so the same blocks are tested.
# Returns one row per block: mean connectivity inside (m_PA) and outside (m_OPA) PAs,
# their difference (diff = PA minus OPA) and the block center coordinates in km.

block_diffs <- function(k, block_size, n = 50000, min_n = 20) {
  
  # Draw random cells with data inside PAs and outside PAs, and label them
  i <- sample_vals(r_pa_list[[k]],  n); i$inside <- "PA"
  o <- sample_vals(r_opa_list[[k]], n); o$inside <- "OPA"
  b <- rbind(i, o)
  # Assign each cell to a square block: dividing the coordinates by the block size and
  # rounding down gives a grid position; cells with the same (bx, by) share a block
  b$bx <- floor(b$x / block_size)
  b$by <- floor(b$y / block_size)
  
  b |>
    # Mean connectivity and number of sampled cells per block and group
    group_by(bx, by, inside) |>
    summarise(m = mean(conn), n = n(), .groups = "drop") |>
    # One row per block
    pivot_wider(names_from = inside, values_from = c(m, n)) |>
    # Keep only blocks with enough cells in both groups
    filter(n_PA >= min_n, n_OPA >= min_n) |>
    # Within-block difference (negative = lower connectivity inside PAs)
    mutate(diff = m_PA - m_OPA,
           xk = (bx + 0.5) * block_size / 1000, # block centres in km
           yk = (by + 0.5) * block_size / 1000)
}

# Moran's I of the within-block differences
moran_block <- function(k, block_size) {
  blk <- block_diffs(k, block_size)
  
  # Inverse-distance weight matrix between all pairs of blocks
  w <- 1 / as.matrix(dist(blk[, c("xk", "yk")])) # inverse-distance weights
  diag(w) <- 0
  
  #Moran's I
  mi <- Moran.I(blk$diff, w)
  data.frame(block_km = block_size / 1000, n_blocks = nrow(blk),
             moran_I = mi$observed, expected = mi$expected, p = mi$p.value)
}

# Run for SSP245 and SSP585
moran_results <- do.call(rbind, lapply(names(scen), function(nm) {
  do.call(rbind, lapply(c(250000, 500000, 750000), function(bs) {
    cbind(scenario = nm, moran_block(scen[[nm]], bs))
  }))
}))

moran_results

# Save result
#write.csv(moran_results, file.path(result.folder, "moran_results.csv"))

## Significance test using paired comparison within spatial blocks

# Create function block_test():
# k = scenario (SSP245 or SSP585)
# block_size = side length of the square blocks, in metres
# n = number of random cells drawn inside PAs and outside PAs (separately)
# min_n = minimum number of sampled cells per group for a block to be used

block_test <- function(k, block_size, n = 50000, min_n = 20) {
  
  # Draw random cells with data inside PAs and outside PAs
  i <- sample_vals(r_pa_list[[k]],  n); i$inside <- "PA"
  o <- sample_vals(r_opa_list[[k]], n); o$inside <- "OPA"
  b <- rbind(i, o)
  
  # Assign each cell to a square block from its coordinates
  b$block <- paste(floor(b$x / block_size), floor(b$y / block_size))
  
  blk <- b |>
    # Mean connectivity and number of sampled cells per block and group
    group_by(block, inside) |>
    summarise(m = mean(conn), n = n(), .groups = "drop") |>
    # One row per block
    pivot_wider(names_from = inside, values_from = c(m, n)) |>
    # Keep only blocks with enough cells in both groups
    filter(n_PA >= min_n, n_OPA >= min_n) |>
    # Within-block difference (negative = lower connectivity inside PAs)
    mutate(diff = m_PA - m_OPA)
  
  # Paired T-test, one difference per block
  tt <- t.test(blk$diff)
  data.frame(block_km  = block_size / 1000,
             # number of blocks used
             n_blocks  = nrow(blk),
             # mean PA-minus-OPA difference
             mean_diff = unname(tt$estimate),
             # 95% CI of the mean difference
             ci_low    = tt$conf.int[1],
             ci_high   = tt$conf.int[2],
             p_ttest   = tt$p.value,
             # Wilcoxon signed-rank test as a non-parametric check
             p_wilcox  = wilcox.test(blk$diff, exact = FALSE)$p.value)
}

# Run for both climate scenarios and for three block sizes (250, 500 and 750 km)
block_results <- do.call(rbind, lapply(names(scen), function(nm) {
  do.call(rbind, lapply(c(250000, 500000, 750000), function(bs) {
    cbind(scenario = nm, block_test(scen[[nm]], bs))
  }))
}))

# See result
block_results

# Save result
#write.csv(block_results, file.path(result.folder, "block_results.csv"))

## Test effect size with Cohen's d

# Random sample of cliamte connectivity values in PAs and OPAs
s_in  <- lapply(r_pa_list,  sample_vals, n = n)
s_out <- lapply(r_opa_list, sample_vals, n = n)

# Combine into one data frame per scenario
samples <- Map(function(i, o) {
  i$group <- 1
  o$group <- 0
  rbind(i, o)
}, s_in, s_out)

# Make group a factor
samples <- lapply(samples, function(s) {
  s$group <- factor(s$group, levels = c(1, 0), labels = c("PA", "OPA"))
  s
})

# Rename
names(samples) <- c("SSP245","SSP585")

# Test effect size with Cohen's d
cohens_results <- do.call(rbind, lapply(names(samples), function(nm) {
  d <- cohens_d(conn ~ group, data = samples[[nm]])
  cbind(scenario = nm, as.data.frame(d))
}))

# Save result
#write.csv(cohens_results, file.path(result.folder, "cohens_d_results.csv"))