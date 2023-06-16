library(dismo)
library(ncdf4)
library(dplyr)
library(raster)
library(leaflet)
library(htmlwidgets)
library(sdmpredictors)

setwd("D:/PhD/articles/article SDM/")


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

TSS <- function(eval){
  return(max(eval@TPR+eval@TNR-1))
}


TSSweight <- function(eval){
  if (TSS(eval)-0.7>0){
    (TSS(eval)-0.5)^2
  }else {
    0
  }
}

weights_tsss <- function(eval){
  if (eval-0.7>0){
    (eval-0.5)^2
  }else {
    0
  }
}

species_locations <- function(locations,species_name){
  cranch_locations <- select(locations[locations$species==species_name,],
                             longitude,
                             latitude)
  return(cranch_locations)
}

plotting <- function(rasteR,species,TSS,algo,points){
  mods <- rasteR
  palete = colorNumeric(c("#5E85B8","#EDF0C0","#C13127"), c(0,1), na.color = 'transparent')
  name <- species
  add_title <- algo
  m <-  leaflet()%>%
    addProviderTiles(provider = 'Esri.OceanBasemap') %>%
    addScaleBar()%>%
    addRasterImage(colors = palete, mods, opacity = .6)%>%
    addLegend(position = 'topright',pal = palete, values(mods), title = paste(name,' - ',add_title,', TSS: ',TSS))
  return(m)
}




main <- function(species_name,
                 oegopsids,
                 env,
                 future_env1,
                 future_env2,
                 future_env3,
                 future_env4,
                 future_env5,
                 future_env6,
                 future_env7,
                 future_env8,
                 background_macaronesia_test,
                 background_macaronesia_train.env,
                 macaronesia.extent
                 ){
  cat(paste(species_name),'-',dim(species_locations(oegopsids,species_name))[1])
  presence_locations <- species_locations(oegopsids,species_name)
#  j <- na.omit(extract(env,presence_locations))
#  presence_locations <- presence_locations[-c(attr(x=j, which='na.action')),]
  presence_locations <- spThin(presence_locations,dist=20)[[1]]
  cat('-->',paste(dim(presence_locations)[1],'\n'))
  
  envs <- list(
    env,
    future_env1,
    future_env2,
    future_env3,
    future_env4,
    future_env5,
    future_env6,
    future_env7,
    future_env8
    )
  env_names = list(
    '_Present_',
    '_2050.RCP26_',
    '_2050.RCP45_',
    '_2050.RCP60_',
    '_2050.RCP85_',
    '_2100.RCP26_',
    '_2100.RCP45_',
    '_2100.RCP60_',
    '_2100.RCP85_'
    )
  predictionsmX = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c())
  tssmX <- c()#True Skill score
  ww <- 0
  cont <- 0
  threshesmX <- c()#thresholds
  for (i in 1:5){
    cat(paste(i,'\n'))
    group <- kfold(presence_locations, 4)
    presence_locations.pres_train <- presence_locations[group != 1, ]
    presence_locations.pres_test <- presence_locations[group == 1, ]
    sp.env <- extract(env,presence_locations.pres_train)#species stats
    presence_absence_env <- as.data.frame(rbind(sp.env,background_macaronesia_train.env))
    presence_absence <- c(rep(1,dim(sp.env)[1]),rep(0,dim(background_macaronesia_train.env)[1]))
    ENV_PA <- presence_absence_env
    ENV_PA$presence <- presence_absence
    

    maxEnt <- maxent(x=presence_absence_env,presence_absence,silent=FALSE)
    em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
    tssmX <- c(tssmX,TSS(em))
    w <- TSSweight(em)
    ww <- ww+w
    cont <- cont + maxEnt@results[7:14]*w
    if (i==1){#
      for (ii in 1:length(envs)){
        predictionsmX[[ii]] <- predict(envs[[ii]], maxEnt, ext=macaronesia.extent,type="response")*w
      }
      }else {
        for (ii in 1:length(envs)){
          predictionsmX[[ii]] <- predictionsmX[[ii]] + predict(envs[[ii]], maxEnt, ext=macaronesia.extent,type="response")*w
        }
      }
  }
  
  cat('MaxEnt model\n')
  wmX <- sapply(tssmX, function(x) weights_tsss(x))
  tssmX <- weighted.mean(tssmX, wmX, na.rm=TRUE)
  cont <- cont/ww
  cat(row.names(maxEnt@results)[7:14])
  cat(cont)
  add_title <- 'MaxEnt'
  for (i in 1:length(env_names)){
    predictionsmX[[i]] <- predictionsmX[[i]]/ww
    mods <- predictionsmX[[i]]
    m <-  plotting(mods,species_name,tssmX,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    writeRaster(predictionsmX[[i]], filename=paste(species_name,env_names[[i]],add_title,'.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  }

}

lste <- list.files(pattern='^Present')
env <- stack(lste)
mask_depth <- load_layers(list_layers( datasets="Bio-ORACLE" )[25,]$layer_code)
mask_depth <- mask_depth[[1]]
mask_depth[mask_depth<0] <- NA
mask_depth <-crop(mask_depth,extent(env[[4]]))
crs(env) <- crs(mask_depth)
# env <- mask(env,mask_depth)

lst1 <- list.files(pattern='^2050AOGCM.RCP26.')
future_env1 <- stack(lst1)
crs(future_env1) <- crs(mask_depth)
# future_env1 <- mask(future_env1,mask_depth)
names(future_env1) <- names(env)


lst2 <- list.files(pattern='^2050AOGCM.RCP45.')
future_env2 <- stack(lst2)
crs(future_env2) <- crs(mask_depth)
# future_env2 <- mask(future_env2,mask_depth)
names(future_env2) <- names(env)

lst3 <- list.files(pattern='^2050AOGCM.RCP60.')
future_env3 <- stack(lst3)
crs(future_env3) <- crs(mask_depth)
# future_env3 <- mask(future_env3,mask_depth)
names(future_env3) <- names(env)

lst4 <- list.files(pattern='^2050AOGCM.RCP85.')
future_env4 <- stack(lst4)
crs(future_env4) <- crs(mask_depth)
# future_env4 <- mask(future_env4,mask_depth)
names(future_env4) <- names(env)

lst5 <- list.files(pattern='^2100AOGCM.RCP26.')
future_env5 <- stack(lst5)
crs(future_env5) <- crs(mask_depth)
# future_env5 <- mask(future_env5,mask_depth)
names(future_env5) <- names(env)

lst6 <- list.files(pattern='^2100AOGCM.RCP45.')
future_env6 <- stack(lst6)
crs(future_env6) <- crs(mask_depth)
# future_env6 <- mask(future_env6,mask_depth)
names(future_env6) <- names(env)


lst7 <- list.files(pattern='^2100AOGCM.RCP60.')
future_env7 <- stack(lst7)
crs(future_env7) <- crs(mask_depth)
# future_env7 <- mask(future_env7,mask_depth)
names(future_env7) <- names(env)

lst8 <- list.files(pattern='^2100AOGCM.RCP85.')
future_env8 <- stack(lst8)
crs(future_env8) <- crs(mask_depth)
# future_env8 <- mask(future_env8,mask_depth)
names(future_env8) <- names(env)





background <- randomPoints(env, n=1000, extf = 1.25)
colnames(background) = c('longitude', 'latitude')
group <- kfold(background, 4)
background_train <- background[group != 1, ]
background_train.env <- extract(env, background_train)
background_test <- background[group == 1, ]

background_macaronesia_test <-background_test
background_macaronesia_train.env<-background_train.env



# #collinearities
# # setwd("D:/PhD/articles/article SDM/results/all_variables/")
# temp.med.suf <- raster('Present.Surface.Temperature.Mean.asc')
# temp.med.bed <- raster('Present.Benthic.Mean.Depth.Temperature.Mean.asc')
# 
# sal.med.suf <- raster('Present.Surface.Salinity.Mean.asc')
# 
# clo.med.suf <- raster('Present.Surface.Chlorophyll.Mean.asc')
# clo.min.suf <- raster('Present.Surface.Chlorophyll.Min.asc')
# 
# cur.med.suf <- raster('Present.Surface.Current.Velocity.Mean.asc.BOv2_1.asc')
# cur.med.bed <- raster('Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1.asc')
# bboxx <- extent(-180, 180, -79.5, 90)
# ice.med <- crop(raster('Present.Surface.Ice.thickness.Mean.asc'),bboxx)
# writeRaster(ice.med,overwrite=T,filename='Present.Surface.Ice.thickness.Mean.asc')
# ice.max <- crop(raster('Present.Surface.Ice.thickness.Max.asc'),bboxx)
# writeRaster(ice.max,overwrite=T,filename='Present.Surface.Ice.thickness.Max.asc')
# ice.min <- crop(raster('Present.Surface.Ice.thickness.Min.asc'),bboxx)
# writeRaster(ice.min,overwrite=T,filename='Present.Surface.Ice.thickness.Min.asc')

# env_vars <- list(as.matrix(temp.med.suf),
#                 as.matrix(temp.med.bed),
#                 as.matrix(sal.med.suf),
#                 as.matrix(clo.med.suf),
#                 as.matrix(clo.min.suf),
#                 as.matrix(ice.med),
#                 as.matrix(ice.max),
#                 as.matrix(ice.min)
#                 )
# for (i in seq(1,length(env_vars),1)){
#   for (ii in (i+1):length(env_vars)){
#     x <- c(env_vars[i][[1]])
#     y <- c(env_vars[ii][[1]])
#     result <- cor(x,y, method = "pearson", use="na.or.complete")
#     cat(result,' ')
#   }
#   cat('\n\n')
# }

species_list <- c("Abraliopsis gilchristi",
                  "Alluroteuthis antarcticus",
                  "Bathyteuthis abyssicola",
                  "Galiteuthis glacialis",
                  "Gonatus antarcticus",
                  "Histioteuthis atlantica",
                  "Histioteuthis eltaninae",
                  "Kondakovia longimana",
                  "Lycoteuthis lorigera",
                  "Martialia hyadesi",
                  "Mesonychoteuthis hamiltoni",
                  "Moroteuthis ingens",
                  "Moroteuthis robsoni",
                  "Psychroteuthis glacialis",
                  "Slosarczykovia circumantarctica",
                  "Teuthowenia pellucida",
                  "Todarodes filippovae"
)


points_collection <- read.csv('curated_occurrences.csv')
macaronesia.extent <- extent(-180, 180, -79.5, -30)
for (species_name in species_list) {
  points_set <- select(points_collection[points_collection$species==species_name,],
                        species,
                        longitude,
                        latitude)
  main(species_name,
    points_set, env,
    future_env1,
    future_env2,
    future_env3,
    future_env4,
    future_env5,
    future_env6,
    future_env7,
    future_env8,background_test,background_train.env,macaronesia.extent)
}
