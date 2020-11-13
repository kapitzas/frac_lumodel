gc()
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
raw_path <-  file.path(path.expand("~"), "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data")

#lu2 <- list.files(raw_path, pattern = "nc", full.names = TRUE)

#for (i in 1:length(lu2)){
#system(paste0("gdalwarp -of Gtiff -co COMPRESS=LZW -co TILED=YES -ot Byte -te -180.0000000 -90.0000000 #180.0000000 90.0000000 -tr 0.002777777777778 0.002777777777778 -t_srs EPSG:4326 NETCDF:", paste0("'", #lu2[[i]], "'"), ":lccs_class ",  paste0("'", gsub(".nc", ".tif", lu2[[i]]), "'")))
#}

##########################
#### 1.) Land use data####
##########################

lu <- stack(file.path(raw_path, "Global", "CLC maps", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif"))
boundary <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/Amazon boundaries/Lim_Biogeografico.shp")

input <- file.path(raw_path, "Global", "CLC maps", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif")

out <- file.path(temp_path, "lu.ama.tif")

if(!dir.exists(temp_path)){
  dir.create(temp_path)
}
  
reproj_ras(input, out, crs = crs(lu), res = res(lu), method = "near", ext = extent(boundary))

#saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

lu <- stack(list.files(file.path(raw_path, "Global", "CLC maps"), pattern = "v2.1", full.names = TRUE))

for (i in 1:nlayers(lu)){
  lu_path <- lu[[i]]@file@name
  out <- file.path(temp_path, paste0("lu_", names(lu[[i]]), ".tif"))
  reproj_ras(lu_path, out, crs = crs(lu), res = res(lu), method = "near", ext = extent(boundary))
}

ts <- c(1, 5, 10, 15, 20, 25, 27)
yrs <- c(1992:2018)

lu <- stack(c(list.files(temp_path, pattern = "1.tif", full.names = TRUE), list.files(temp_path, pattern = "lu.ama.tif", full.names = TRUE)))
names(lu) <- paste0("Y", c(2016:2018, 1992:2015))

lu <- lu[[order(names(lu))]]

mask <- raster("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data/data_ama/temp/ts1992_l01_ama.tif")
writeRaster(lu[[ts]], bylayer = TRUE, filename = paste0(temp_path, "/lu_cropped_", names(lu[[ts]])), format = "GTiff", overwrite = TRUE)

lu <- stack(list.files(temp_path, pattern = "lu_cropped", full.names = TRUE))

logfile <- file.path(data_path, paste0("amazon_log.txt"))
writeLines(c(""), logfile)

uniq <- list()
for(i in 1:nlayers(lu)){
  uniq[[i]] <- unique(lu[[i]][])
}
unique(unlist(uniq))

crop <- c(10, 11, 12, 20, 30)
crop_mosaic <- 40
forest <- c(50, 60, 61, 62, 70, 80, 90, 100, 160, 170) #>15% cover
grass <- c(110, 130)
wetland <- c(180)
urban <- c(190)
shrub <- c(120, 121,122)
other <- c(140, 150, 151, 152, 153, 200,201,202,220)
water <- 210

ts <- c(1, 5, 10, 15, 20, 25, 27)
yrs <- c(1992:2018)[ts]
for(j in 1:nlayers(lu)){
  cat(paste0("Time step ", yrs[j], "\n"), file = logfile, append = TRUE)
  cat(paste0("Reclassifying...",   "\n"), file = logfile, append = TRUE)
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
  
  cat(paste0("Layerizing...",   "\n"), file = logfile, append = TRUE)
  layers <- stack(layerize(r))
  
  names(layers) <- paste0("X", str_pad(1:nlayers(layers), 2, pad = "0"))
  cat(paste0("Aggregating...",   "\n"), file = logfile, append = TRUE)
  writeRaster(layers, bylayer = TRUE, filename = paste0(temp_path, "/temp_", names(layers)), format = "GTiff", overwrite = TRUE)
  
  removeTmpFiles(h = 0)
  files <- list.files(temp_path, pattern = "temp_", full.names = TRUE)
  
  for(l in 1:length(files)){
    out <- file.path(temp_path, paste0(paste0("ts", str_pad(yrs[j], 2, pad = "0")), "_", paste0("l", str_pad(l, 2, pad = "0")), "_ama.tif"))
    reproj_ras(files[l], out, crs = crs(mask), res = 0.00833333, method = "average", ext = extent(mask))
  }
  unlink(files)
}

unlink(file.path(temp_path, "lu.ama.tif"))

####################
#### 2. MASK #######
####################

mask <- raster(file.path(temp_path, "ts1992_l01_ama.tif"))
mask[!is.na(mask[])] <- 1
mask <-raster::mask(mask, boundary)
plot(mask)
saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

####################
#### 3. BIO ########
####################

biofiles <- list.files(file.path(raw_path, "Global", "wc21_30s_bio"), full.names = TRUE)
mask <- readRDS(file.path(data_path, "mask_ama.rds"))
bio_names <- paste0("bio", c(1, 10:19, 2:9))

for(i in 1:length(biofiles)){
  infile <- biofiles[[i]]
  outfile <- file.path(temp_path, paste0(bio_names[i], "_ama.tif"))
  reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask)[c(1,2,3,4)], res = res(mask), method = "bilinear")
}

##########################
####4. Protected Areas####
##########################
infile <- file.path(raw_path, "Global", "Global Protected Areas", "WDPA_Mar2018-shapefile-polygons.shp")
outfile <- file.path(temp_path, "PA_cropped.shp")
crop_shp(infile, outfile, ext = extent(mask))
test <- sf::st_read(outfile)
#plot(test[which(test$IUCN_CAT%in%c("Ia", "Ib", "II")),], max.plot = 1, add = TRUE)

years <- yrs[ts]
for (i in 1:length(years)){
  PA <- test[which(test$IUCN_CAT%in%c("Ia", "Ib", "II") & test$STATUS_YR <= years[i]),]
  r <- rasterize(PA, mask)
  out <- mask
  out[which(!is.na(r[]))] <- 0
  out <- raster::mask(out, mask)
  
  writeRaster(out, file.path(temp_path, paste0("PA", years[i], "_ama.tif")), format = "GTiff", overwrite = TRUE)
}

####################
####4. Elevation####
####################

mask <- readRDS(file.path(data_path, "mask_ama.rds"))
infile <-file.path(raw_path, "Global", "wc2.1_30s_elev.tif")
r <- raster(infile)

unlink(file.path(temp_path, "srtm_ama.tif"))
outfile <- file.path(temp_path,  "srtm_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#slope and roughness
srtm <- raster(outfile)
slope <- terrain(srtm, "slope")
roughness <- terrain(srtm, "roughness")
writeRaster(slope, filename = file.path(temp_path, paste0("slope_ama.tif")), driver = "GTiff", overwrite = TRUE)
writeRaster(roughness, filename = file.path(temp_path, paste0("roughness_ama.tif")), driver = "GTiff", overwrite = TRUE)

###########################
####5. Distance rasters####
###########################

#Roads
# infile <- paste0(file.path(raw_path, "Global", "groads-v1-americas-shp/gROADS-v1-americas.shp"))
# unlink(file.path(temp_path, "roads_raster.tif"))
# outfile <- file.path(temp_path, "roads_raster.tif")
# rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])
# 
# infile <- file.path(temp_path, "roads_raster.tif")
# unlink(file.path(temp_path, "diro_ama.tif"))
# outfile <- file.path(temp_path, "diro_ama.tif")
# proximity_ras(infile, outfile)

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
  #ggplot() +  geom_sf(data = boundary, fill = "darkred") + geom_sf(data = test)

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

#Built-up areas
infile <- file.path(raw_path, "Global", "Global Built up areas", "bltupa.shp")
unlink(file.path(temp_path, "builtup_raster.tif"))
outfile <- file.path(temp_path, "builtup_raster.tif")
rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])

infile <- file.path(temp_path, "builtup_raster.tif")
unlink(file.path(temp_path, "dibu_ama.tif"))
outfile <- file.path(temp_path, "dibu_ama.tif")
proximity_ras(infile, outfile)

#Lakes
infile <- file.path(raw_path, "Global", "GSHHS data", "GSHHS_lakes_L2-L4.shp")
unlink(file.path(temp_path, "dila_raster.tif"))
outfile <- file.path(temp_path, "dila_raster.tif")
rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])

infile <- file.path(temp_path, "dila_raster.tif")
unlink(file.path(temp_path, "dila_ama.tif"))
outfile <- file.path(temp_path, "dila_ama.tif")
proximity_ras(infile, outfile)

#Rivers
infile <- file.path(raw_path, "Global", "GSHHS data", "WDBII_rivers_global_L2-L9.shp")
unlink(file.path(temp_path, "diri_raster.tif"))
outfile <- file.path(temp_path, "diri_raster.tif")
rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])

infile <- file.path(temp_path, "diri_raster.tif")
unlink(file.path(temp_path, "diri_ama.tif"))
outfile <- file.path(temp_path, "diri_ama.tif")
proximity_ras(infile, outfile)

################
####6. Soils####
################
#Full description: https://www.isric.org/explore/soilgrids/faq-soilgrids

#Organic Carbon Density
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "BLDFIE_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "ocdens_eur.tif"))
outfile <- file.path(temp_path,  "ocdens_eur.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#Available Soil Water Capacity
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "WWP_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "wwp_ama.tif"))
outfile <- file.path(temp_path,  "wwp_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#pH Index measured in Water Solution
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "PHIHOX_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "phihox_ama.tif"))
outfile <- file.path(temp_path,  "phihox_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#Bulk density fine earth
infile <- file.path(raw_path, "Global", "soil_data", "ISRIC", "BLDFIE_M_sl3_1km_ll.tif")
unlink(file.path(temp_path, "bldfie_ama.tif"))
outfile <- file.path(temp_path,  "bldfie_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#############################
####7. Population Density####
#############################

#infile <- file.path(raw_path, "Global", "Pop_density", "gluds00ag.bil")
#unlink(file.path(temp_path, "popdens_ama.tif"))
#outfile <- file.path(temp_path,  "popdens_ama.tif")
#reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#####################
####8. Processing####
#####################

#Check that rasters have same extent
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


#Synch NA and output mask
mask <- readRDS(file.path(data_path, "mask_ama.rds"))
for(i in 1:length(fl)){
  r <- raster(fl[[i]])
  mask <- raster::mask(mask, r)
  print(i)
}

plot(mask)

saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

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
#saveRDS(lu_mat, file = file.path(data_path, "lu.rds"), compress = TRUE)

#Getting rid of ice/snow class (class 12) and renaming water to class class 12
rm(lu_mat)
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

# 
# ice_cols <- seq(9, 8*8, 9)
# ice_rows <- which(rowSums(lu_mat[,ice_cols])!=0)
# lu_mat <- lu_mat[-ice_rows,-ice_cols]
# mask[inds[ice_rows]] <- NA
# inds <- inds[-ice_rows]
# saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

#rename land use class 13 to land use class 12
#colnames(lu_mat) <- colnames[-seq(13, 13*8, 13)]

#plot(raster("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data/data_ama/temp/diro_loc_ama.tif"))

