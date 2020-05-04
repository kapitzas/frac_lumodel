gc()
require(devtools)
require(WorldClimTiles)
require(gdaltools)
require(rgdal)
require(raster)
require(stringr)
require(SpaDES)

data_path <- file.path(getwd(), "data", "data_ama")
temp_path <- file.path(getwd(), "data" , "data_ama", "temp")
raw_path <-  file.path(path.expand("~"), "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data")


# lu2 <- list.files(raw_path, pattern = "nc", full.names = TRUE)
# 
# for (i in 1:length(lu2)){
# system(paste0("gdalwarp -of Gtiff -co COMPRESS=LZW -co TILED=YES -ot Byte -te -180.0000000 -90.0000000 180.0000000 90.0000000 -tr 0.002777777777778 0.002777777777778 -t_srs EPSG:4326 NETCDF:", paste0("'", lu2[[i]], "'"), ":lccs_class ",  paste0("'", gsub(".nc", ".tif", lu2[[i]]), "'")))
# }

lu <- stack(file.path(raw_path, "Global", "CLC maps", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif"))
boundary <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/Amazon boundaries/Lim_Biogeografico.shp")

input <- file.path(raw_path, "Global", "CLC maps", "ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif")

out <- file.path(temp_path, "lu.ama.tif")
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

lu <- stack(list.files(temp_path, pattern = "tif", full.names = TRUE))
names(lu) <- paste0("Y", c(2016:2018, 1992:2015))

mask <- lu[[1]]
writeRaster(lu, bylayer = TRUE, filename = paste0(temp_path, "/lu_cropped_", names(lu)), format = "GTiff", overwrite = TRUE)

lu <- stack(list.files(temp_path, pattern = "lu_cropped", full.names = TRUE))

logfile <- file.path(data_path, paste0("amazon_log.txt"))
writeLines(c(""), logfile)

uniq <- list()
for(i in 1:nlayers(lu)){
  uniq[[i]] <- unique(lu[[i]][])
}
table(unlist(uniq))
all(unlist(lapply(uniq, FUN = length)) == 26)

crop <- c(10, 11, 12, 20)
crop_natveg_mosaic <- c(30, 40)
tree_co <- c(50, 60, 61, 62, 70, 80) #>15% cover
tree_mixed <- c(90) #mixed leaf
treeshrub_herb_mosaic <- c(100,110)
shrub <- c(120,122)
grass <- c(130)
sparse <- c(150, 153)
tree_water <- c(160,170)
shrub_water <- c(180)
urban <- c(190)
bare <- c(200)
water <- 210
snow_ice <- c(220)


for(j in ts){
  cat(paste0("Time step ", yrs[j], "\n"), file = logfile, append = TRUE)
  
  cat(paste0("Reclassifying...",   "\n"), file = logfile, append = TRUE)
  r <- lu[[j]]
  r[r%in%crop] <- 1
  r[r%in%crop_natveg_mosaic] <- 2
  r[r%in%tree_co] <- 3
  r[r%in%tree_mixed] <- 4 #mixed leaf
  r[r%in%treeshrub_herb_mosaic] <- 5
  r[r%in%shrub] <- 6
  r[r%in%grass] <- 7
  r[r%in%sparse] <- 8
  r[r%in%tree_water] <- 9
  r[r%in%shrub_water] <- 10
  r[r%in%urban] <- 11
  r[r%in%bare] <- 12
  r[r%in%snow_ice] <- 13
  r[r%in%water] <- NA
  
  cat(paste0("Layerizing...",   "\n"), file = logfile, append = TRUE)
  layers <- stack(layerize(r))
  names(layers) <- paste0("X", str_pad(1:nlayers(layers), 2, pad = "0"))
  cat(paste0("Aggregating...",   "\n"), file = logfile, append = TRUE)
  
  writeRaster(layers, bylayer = TRUE, filename = paste0(temp_path, "/temp_", names(layers)), format = "GTiff", overwrite = TRUE)
  
  removeTmpFiles(h = 0)
  files <- list.files(temp_path, pattern = "temp_", full.names = TRUE)
  
  for(l in 1:length(files)){
    out <- file.path(temp_path, paste0(paste0("ts", str_pad(yrs[j], 2, pad = "0")), "_", paste0("l", str_pad(l, 2, pad = "0")), "_ama.tif"))
    reproj_ras(files[l], out, crs = crs(mask), res = 0.0833333, method = "average", ext = extent(mask))
  }
  unlink(files)
}

unlink(file.path(temp_path, "lu.ama.tif"))
####################
###PREPROCESSING####
####################
mask <- raster("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data_ama/temp/ts1992_l01_ama.tif")
mask[!is.na(mask[])] <- 1
mask <-mask(mask, boundary)
plot(mask)
saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

##########################
#### 1.) Land use data####
##########################
#unlink(list.files(temp_path, full.names = TRUE), recursive = TRUE)

mask <- readRDS(file.path(data_path, "mask_ama.rds"))
tiles <- tile_name(mask, "worldclim")
bio <- tile_get(tiles = tiles, var = "bio", path = temp_path)
bio <- tile_merge(bio)
bio <- crop(bio, mask)

bio_names <- names(bio)
infile <- bio@file@name
i <- 1
for(i in 1:nlayers(bio)){
  infile <- file.path(temp_path, "temp_bio.tif")
  writeRaster(bio[[i]], infile, format = "GTiff", bylayer = TRUE)
  outfile <- file.path(temp_path, paste0(bio_names[i], "_ama.tif"))
  reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask)[c(1,2,3,4)], res = res(mask), method = "bilinear")
  unlink(infile)
}


####################
####4. Elevation####
####################

mask <- readRDS(file.path(data_path, "mask_ama.rds"))
infile <-file.path(raw_path, "Global", "topo30", "topo30.grd")
r <- raster(infile)
r <- rotate(r)

writeRaster(r, gsub("grd", "tif", infile), driver = "GTiff")
infile <- gsub("grd", "tif", infile)

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
infile <- paste0("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/groads-v1-americas-shp/gROADS-v1-americas.shp")
unlink(file.path(temp_path, "roads_raster.tif"))
outfile <- file.path(temp_path, "roads_raster.tif")
rasterize_shp(infile, outfile, res = res(mask)[1], ext = extent(mask)[c(1,2,3,4)])

infile <- file.path(temp_path, "roads_raster.tif")
unlink(file.path(temp_path, "diro_ama.tif"))
outfile <- file.path(temp_path, "diro_ama.tif")
proximity_ras(infile, outfile)

t <- mask
t[which(is.na(mask[]))] <- dat[,1]

plot(mask)
nrow(lu)
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

infile <- file.path(raw_path, "Global", "Pop_density", "gluds00ag.bil")
unlink(file.path(temp_path, "popdens_ama.tif"))
outfile <- file.path(temp_path,  "popdens_ama.tif")
reproj_ras(infile, outfile, crs = crs(mask), ext = extent(mask), res = res(mask), method = "near")

#####################
####8. Processing####
#####################

#Check that rasters have same extent
fl <- list.files(temp_path, "ama.tif", full.names = TRUE)
r <- list()
for(i in 1:length(fl)){
  print(i)
  r[[i]] <- extent(raster(fl[[i]]))
  if(length(unique(r))!=1){
    print(i)
    break
  }
}
length(r) == length(fl)

#Synch NA and output mask
mask <- readRDS(file.path(data_path, "mask_ama.rds"))

for(i in 1:length(fl)){
  r <- raster(fl[[i]])
  mask <- raster::mask(mask, r)
  print(i)
}

saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))
removeTmpFiles(h=0)
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

rm(lu_mat)

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

data_path <- file.path(getwd(), "data", "data_ama")
dat <- readRDS(file.path(data_path, "cov.rds")) #dynamic bioclimatic variables
lu_all <- readRDS(file.path(data_path, "lu.rds"))
mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask

mask_cols <- c(grep(pattern = "l04", colnames(lu_all)), grep(pattern = "l13", colnames(lu_all)))
mask_rows <- which(rowSums(lu_all[, mask_cols]) > 0) #get rid of lu classes that have less than 10 cells with values.
dat <- dat[-mask_rows,]
lu_all <- lu_all[-mask_rows, -mask_cols]
saveRDS(lu_all, file.path(data_path, "lu.rds"))
saveRDS(dat, file.path(data_path, "cov.rds"))
mask[which(!is.na(mask[]))[mask_rows]] <- NA
saveRDS(readAll(mask), file.path(data_path, "mask_ama.rds"))

colnames(cov_mat) <- colnames
saveRDS(cov_mat, file = file.path(data_path, "cov.rds"), compress = TRUE)
