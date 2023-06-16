library(dismo)
library(ncdf4)
library(dplyr)
library(raster)


setwd("D:/PhD/articles/article SDM/")
species_locations <- function(locations,species_name){
  cranch_locations <- select(locations[locations$species==species_name,],
                             longitude,
                             latitude)
  return(cranch_locations)
}
spThin <- function(xy,dist,rep=1) {
  #dist -> distance between points in kilometers
  #rep -> replicates
  .r <- raster(resolution=0.008333333) # empty raster with 1km resolution
  .a <- raster(ext=extent(.r),resolution=res(.r)[1]*dist)
  .ac <- cellFromXY(.a,xy)
  .tbl <- data_frame(cell=.ac,x=xy[,1],y=xy[,2]) %>% group_by(cell)
  o <- list()
  for (i in 1:rep) {
    .s <- .tbl %>% sample_n(1) %>% ungroup()
    o[[i]] <- data.frame(.s %>% dplyr::select(x,y))
  }
  o
}

raster <- raster(filename)
species_name
species_points <- spThin(species_locations(read.csv("occorrences.csv"),species_name),dist=20)
presences <- extract(raster,species_points)
background <- randomPoints(env, n=1000, extf = 1.25)
colnames(background) = c('longitude', 'latitude')

absences <- extract(raster,background)
e <- evaluate(p=presences,a=absences)
t <- threshold(e,stat='prevalence')
binary <- raster>t