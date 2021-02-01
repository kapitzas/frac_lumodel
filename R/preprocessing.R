rm(list = ls())
require(devtools)
require(WorldClimTiles)
require(rgdal)
require(gdaltools)
require(raster)
require(stringr)
require(sf)

data_path <- file.path(getwd(), "data", "data_ama_1k")
temp_path <- file.path(data_path, "temp")
if(!dir.exists(temp_path)){
  dir.create(temp_path)
}
raw_path <-  file.path(path.expand("~"), "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data")

ts <- c(1, 5, 10, 15, 20, 25, 27)
yrs <- c(1992:2018)[ts]

#-----------------------#
#### 1. Land use data####
#-----------------------#

# 1. a load clc data and Amazon boundary file

# boundary 

boundary <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/Amazon boundaries/Lim_Biogeografico.shp")

# lu 1992-2015
lu <- stack(file.path(raw_path, "Global", "CLC maps", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif"))

# crop clc data to boundary extent  
input <- file.path(raw_path, "Global", "CLC maps", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif")
out <- file.path(temp_path, "lu.ama.tif")
reproj_ras(input, out, crs = crs(lu), res = res(lu), method = "near", ext = extent(boundary))

# do the same for lu 2015-2018 (this comes as a seperate repo)
lu <- stack(list.files(file.path(raw_path, "Global", "CLC maps"), pattern = "v2.1", full.names = TRUE))

for (i in 1:nlayers(lu)){
  lu_path <- lu[[i]]@file@name
  out <- file.path(temp_path, paste0("lu_", names(lu[[i]]), ".tif"))
  reproj_ras(lu_path, out, crs = crs(lu), res = res(lu), method = "near", ext = extent(boundary))
}

# load cropped data as stack and bring into correct order
lu <- stack(c(list.files(temp_path, pattern = "1.tif", full.names = TRUE), list.files(temp_path, pattern = "lu.ama.tif", full.names = TRUE)))
names(lu) <- paste0("Y", c(2016:2018, 1992:2015))

lu <- lu[[order(names(lu))]]

# write cropped rasters to disk (to prevent memory issues)
writeRaster(lu[[ts]], bylayer = TRUE, filename = paste0(temp_path, "/lu_cropped_", names(lu[[ts]])), format = "GTiff", overwrite = TRUE)

# load a mask raster for later
mask <- raster("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data/data_ama/temp/ts1992_l01_ama.tif")

# reallocate clc land use classes to new classes
crop <- c(10, 11, 12, 20, 30)
crop_mosaic <- 40
forest <- c(50, 60, 61, 62, 70, 80, 90, 100, 160, 170) #>15% cover
grass <- c(110, 130)
wetland <- c(180)
urban <- c(190)
shrub <- c(120, 121,122)
other <- c(140, 150, 151, 152, 153, 200,201,202,220)
water <- 210

# reload cropped lu data (this validates they exist at the location)
lu <- stack(list.files(temp_path, pattern = "lu_cropped", full.names = TRUE))

# processing lu data takes a long time - the logfile keeps track of where issues occurred if this fails. We need to do this sequentially because of some memory issues.
logfile <- file.path(data_path, paste0("amazon_log.txt"))
writeLines(c(""), logfile)

# loop through time steps (we do validation on 5 year intervals to save computational time)
for(j in 1:nlayers(lu)){
  cat(paste0("Time step ", yrs[j], "\n"), file = logfile, append = TRUE)
  cat(paste0("Reclassifying...",   "\n"), file = logfile, append = TRUE)
  
  # reclassify super-high res categorical map
  r <- lu[[j]]
  r[r%in%crop] <- 1
  r[r%in%crop_mosaic] <- 2
  r[r%in%forest] <- 3
  r[r%in%grass] <- 4
  r[r%in%shrub] <- 5
  r[r%in%wetland] <- 6
  r[r%in%urban] <- 7
  r[r%in%other] <- 8
  r[r%in%water] <- 9
  
  # layerize the multi-cat map (to 9 0/1 layers) and write to save RAM
  cat(paste0("Layerizing...",   "\n"), file = logfile, append = TRUE)
  layers <- stack(layerize(r))
  
  names(layers) <- paste0("X", str_pad(1:nlayers(layers), 2, pad = "0"))
  cat(paste0("Aggregating...",   "\n"), file = logfile, append = TRUE)
  writeRaster(layers, bylayer = TRUE, filename = paste0(temp_path, "/temp_", names(layers)), format = "GTiff", overwrite = TRUE)
  removeTmpFiles(h = 0)
  
  # calculate fractional rasters by aggregating layers to target resolution (approx 10km) and calculating averages
  files <- list.files(temp_path, pattern = "temp_", full.names = TRUE)
  for(l in 1:length(files)){
    out <- file.path(temp_path, paste0(paste0("ts", str_pad(yrs[j], 2, pad = "0")), "_", paste0("l", str_pad(l, 2, pad = "0")), "_ama.tif"))
    reproj_ras(files[l], out, crs = crs(mask), res = 0.00833333, method = "average", ext = extent(mask))
  }
  unlink(files)
}

unlink(file.path(temp_path, "lu.ama.tif"))

#---------------------------#
#### 2. FINAL MASK LAYER ####
#---------------------------#

mask <- raster(file.path(temp_path, "ts1992_l01_ama.tif"))
mask[!is.na(mask[])] <- 1
mask <-raster::mask(mask, boundary)
plot(mask)
saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

#-----------------------------#
#### 3. BIOClIMATIC LAYERS ####
#-----------------------------#

# load data
biofiles <- list.files(file.path(raw_path, "Global", "wc21_30s_bio"), full.names = TRUE)

# load mask
mask <- readRDS(file.path(data_path, "mask_ama.rds"))

# reproject bioclim layers to mask
bio_names <- paste0("bio", c(1, 10:19, 2:9))

for(i in 1:length(biofiles)){
  infile <- biofiles[[i]]
  outfile <- file.path(temp_path, paste0(bio_names[i], "_ama.tif"))
  reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask)[c(1,2,3,4)], res = res(mask), method = "bilinear")
}

#--------------------------#
#### 4. Protected Areas ####
#--------------------------#
infile <- file.path(raw_path, "Global", "Global Protected Areas", "WDPA_Mar2018-shapefile-polygons.shp")
outfile <- file.path(temp_path, "PA_cropped.shp")
crop_shp(infile, outfile, ext = extent(mask))
pas <- sf::st_read(outfile)

years <- yrs[ts]
for (i in 1:length(years)){
  PA <- pas[which(pas$IUCN_CAT%in%c("Ia", "Ib", "II") & pas$STATUS_YR <= years[i]),]
  r <- rasterize(PA, mask)
  out <- mask
  out[which(!is.na(r[]))] <- 0
  out <- raster::mask(out, mask)
  
  writeRaster(out, file.path(temp_path, paste0("PA", years[i], "_ama.tif")), format = "GTiff", overwrite = TRUE)
}

#------------------#
#### 4. Terrain ####
#------------------#

mask <- readRDS(file.path(data_path, "mask_ama.rds"))
infile <-file.path(raw_path, "Global", "wc2.1_30s_elev.tif")
r <- raster(infile)

# ELevation
unlink(file.path(temp_path, "srtm_ama.tif"))
outfile <- file.path(temp_path,  "srtm_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

# slope and roughness
srtm <- raster(outfile)
slope <- terrain(srtm, "slope")
roughness <- terrain(srtm, "roughness")
writeRaster(slope, filename = file.path(temp_path, paste0("slope_ama.tif")), driver = "GTiff", overwrite = TRUE)
writeRaster(roughness, filename = file.path(temp_path, paste0("roughness_ama.tif")), driver = "GTiff", overwrite = TRUE)

#--------------------------#
#### 5. Distance layers ####
#--------------------------#

# Distance to roads
# We devide the roads feature file into individual distance rasters, one for each road type (they are liekly to have different influence on land use suitability)
infile <- paste0(file.path(raw_path, "Global", "groads-v1-americas-shp/gROADS-v1-americas.dbf"))
roads <- st_read(infile)
rid <- matrix(c(0:7, "hwy", "pri", "sec", "tert", "loc", "trail", "priv", "unspec"), ncol = 2, nrow = 8)
rid <- data.frame(rid)
road_classes <- sort(unique(roads$FCLASS))
road_classes <- road_classes[-c(5, 3)]
for (i in road_classes){
  print(paste0("Writing subset ", i))
  out <- roads[which(roads$FCLASS == i),]
  outfile <- file.path(temp_path, paste0("roads_", rid[i+1,2], "_raster.shp"))
  test <- st_crop(out, mask)

  if(nrow(test) == 0){
    message("road type not in study area")
    next}
  st_write(out, outfile)
  
  print(paste0("Rasterizing subset ", i))
  infile <- outfile
  outfile <- file.path(temp_path, paste0("roads_", rid[i+1,2], "_raster.tif"))
  rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])
  unlink(infile)
  
  print(paste0("Calculating subset ", i))
  infile <- outfile
  outfile <- file.path(temp_path, paste0("diro_", rid[i+1,2], "_ama.tif"))
  proximity_ras(infile, outfile)
  unlink(infile)
}

# Distance to built-up areas
infile <- file.path(raw_path, "Global", "Global Built up areas", "bltupa.shp")
unlink(file.path(temp_path, "builtup_raster.tif"))
outfile <- file.path(temp_path, "builtup_raster.tif")
rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])

infile <- file.path(temp_path, "builtup_raster.tif")
unlink(file.path(temp_path, "dibu_ama.tif"))
outfile <- file.path(temp_path, "dibu_ama.tif")
proximity_ras(infile, outfile)

# Distance to lakes
infile <- file.path(raw_path, "Global", "GSHHS data", "GSHHS_lakes_L2-L4.shp")
unlink(file.path(temp_path, "dila_raster.tif"))
outfile <- file.path(temp_path, "dila_raster.tif")
rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])

infile <- file.path(temp_path, "dila_raster.tif")
unlink(file.path(temp_path, "dila_ama.tif"))
outfile <- file.path(temp_path, "dila_ama.tif")
proximity_ras(infile, outfile)

#Distance to rivers
infile <- file.path(raw_path, "Global", "GSHHS data", "WDBII_rivers_global_L2-L9.shp")
unlink(file.path(temp_path, "diri_raster.tif"))
outfile <- file.path(temp_path, "diri_raster.tif")
rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])

infile <- file.path(temp_path, "diri_raster.tif")
unlink(file.path(temp_path, "diri_ama.tif"))
outfile <- file.path(temp_path, "diri_ama.tif")
proximity_ras(infile, outfile)

#----------------#
#### 6. Soils ####
#----------------#
#Full description: https://www.isric.org/explore/soilgrids/faq-soilgrids

# Organic Carbon Density
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "BLDFIE_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "ocdens_eur.tif"))
outfile <- file.path(temp_path,  "ocdens_eur.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

# Available Soil Water Capacity
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "WWP_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "wwp_ama.tif"))
outfile <- file.path(temp_path,  "wwp_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

# pH Index measured in Water Solution
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "PHIHOX_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "phihox_ama.tif"))
outfile <- file.path(temp_path,  "phihox_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

# Bulk density fine earth
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "BLDFIE_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "bldfie_ama.tif"))
outfile <- file.path(temp_path,  "bldfie_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#---------------------------#
#### 7. Final processing ####
#---------------------------#

# Check that rasters have same extent
fl <- list.files(temp_path, "ama.tif$", full.names = TRUE)
r <- list()
for(i in 1:length(fl)){
  print(i)
  r[[i]] <- extent(raster(fl[[i]]))
  if(length(unique(r))!=1){
    print(i)
    break
  }
}


# Synch NA and output new mask
mask <- readRDS(file.path(data_path, "mask_ama.rds"))
for(i in 1:length(fl)){
  r <- raster(fl[[i]])
  mask <- raster::mask(mask, r)
  print(i)
}

saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

# create NA-synched land use data matrix
mask <- readRDS(file.path(data_path, "mask_ama.rds"))
inds <- which(!is.na(mask[]))
fl_lu <- fl[grepl("ts", fl)]
lu_mat <- matrix(data = NA, nrow = length(inds), ncol = length(fl_lu))
colnames <- character()

for(i in 1:length(fl_lu)){
  print(i)
  r <- raster(fl_lu[[i]])
  colnames[i] <- gsub('.{4}$', '', names(r))
  r <- getValues(r)
  r <- r[inds]
  lu_mat[,i] <- r
}

colnames(lu_mat) <- colnames
saveRDS(lu_mat, file = file.path(data_path, "lu.rds"), compress = TRUE)

# Create NA-synched covariate matrix
inds <- which(!is.na(mask[]))
fl_cov <- fl[-which(grepl("ts", fl))]
cov_mat <- matrix(data = NA, nrow = length(inds), ncol = length(fl_cov))
colnames <- character()
for(i in 1:length(fl_cov)){
  print(i)
  r <- raster(fl_cov[[i]])
  colnames[i] <- gsub('.{4}$', '', names(r))
  r <- getValues(r)
  r <- r[inds]
  cov_mat[,i] <- r
}

colnames(cov_mat) <- colnames
saveRDS(cov_mat, file = file.path(data_path, "cov.rds"), compress = TRUE)