library(raster)
setwd("D:/PhD/articles/article SDM/")
macaronesia.extent <- extent(-180, 180, -79.5, -30)#size to cut
listOfiles <- c(list.files(pattern='?.tif'),list.files(pattern='?.asc'))#target rasters to cut
for (file in listOfiles) {
	env <- raster(file)
	env <- crop(env,macaronesia.extent)
	writeRaster(env, filename=file, overwrite=TRUE)
}

macaronesia.extent <- extent(-180, 180, -79.5, 80)#size to expand to
listOfiles <- c(list.files(pattern='?.tif'),list.files(pattern='?.asc'))#target rasters
for (file in listOfiles) {
	env <- raster(file)
	env <- extend(env,macaronesia.extent, value=0.0)#filling the void with 0.0
	writeRaster(env, filename=file, overwrite=TRUE)
}
