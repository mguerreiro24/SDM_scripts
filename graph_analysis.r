library(leaflet)
library(raster)
library(htmlwidgets)
library(dplyr)
library(sdmpredictors)
setwd("D:/PhD/articles/article SDM/")

species_name <- c("Alluroteuthis antarcticus",
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



nameSpecies <- c()
Present <- c()
RCP852100 <- c()
RCP602100 <- c()
RCP452100 <- c()
RCP262100 <- c()
RCP852050 <- c()
RCP602050 <- c()
RCP452050 <- c()
RCP262050 <- c()
for (i in 1:length(species_name)){
  cat(paste(species_name[i],'\n'))
  nameSpecies <- c(nameSpecies,species_name[i])
  mat1.data <- raster(paste(species_name[i],'_Present_MaxEnt.tif',sep=''))
  YY <-ymax(extent(mat1.data))
  
  mat1.w <- mat1.data
  A <- area(mat1.data)
  A <- mask(A,mat1.data)

  writeRaster(mat1.w, filename=paste(species_name[i],'_present','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat1<- as.matrix(mat1.w)
  flow1 <- rowMeans(mat1,na.rm=T)#lines
  A1 <- as.matrix(A)
  A1 <- c(A1)
  mat1<- c(mat1)
  Present <- c(Present,
               sum(mat1*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))
  
  mat_2100.RCP85_.data <- raster(paste(species_name[i],'_2100.RCP85_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2100.RCP85_.data,mat1.w), filename=paste(species_name[i],'_2100.RCP85','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2100.RCP85_ <- as.matrix(mask(mat_2100.RCP85_.data,mat1.w))
  flow_2100.RCP85_ <- rowMeans(mat_2100.RCP85_,na.rm=T) - flow1
  mat_2100.RCP85_ <- c(mat_2100.RCP85_)

  RCP852100 <- c(RCP852100,sum(mat_2100.RCP85_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))

  
  mat_2100.RCP60_.data <- raster(paste(species_name[i],'_2100.RCP60_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2100.RCP60_.data,mat1.w), filename=paste(species_name[i],'_2100.RCP60','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2100.RCP60_ <- as.matrix(mask(mat_2100.RCP60_.data,mat1.w))

  flow_2100.RCP60_ <- rowMeans(mat_2100.RCP60_,na.rm=T) - flow1
  mat_2100.RCP60_ <- c(mat_2100.RCP60_)

  RCP602100 <- c(RCP602100,sum(mat_2100.RCP60_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))
  
  mat_2100.RCP45_.data <- raster(paste(species_name[i],'_2100.RCP45_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2100.RCP45_.data,mat1.w), filename=paste(species_name[i],'_2100.RCP45','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2100.RCP45_ <- as.matrix(mask(mat_2100.RCP45_.data,mat1.w))
  
  flow_2100.RCP45_ <- rowMeans(mat_2100.RCP45_,na.rm=T) - flow1
  mat_2100.RCP45_ <- c(mat_2100.RCP45_)

  RCP452100 <- c(RCP452100,sum(mat_2100.RCP45_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))

  mat_2100.RCP26_.data <- raster(paste(species_name[i],'_2100.RCP26_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2100.RCP26_.data,mat1.w), filename=paste(species_name[i],'_2100.RCP26','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2100.RCP26_ <- as.matrix(mask(mat_2100.RCP26_.data,mat1.w))

  flow_2100.RCP26_ <- rowMeans(mat_2100.RCP26_,na.rm=T) - flow1
  mat_2100.RCP26_ <- c(mat_2100.RCP26_)
  
  RCP262100 <- c(RCP262100,sum(mat_2100.RCP26_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))
  
  mat_2050.RCP85_.data <- raster(paste(species_name[i],'_2050.RCP85_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2050.RCP85_.data,mat1.w), filename=paste(species_name[i],'_2050.RCP85','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2050.RCP85_ <- as.matrix(mask(mat_2050.RCP85_.data,mat1.w))

  flow_2050.RCP85_ <- rowMeans(mat_2050.RCP85_,na.rm=T) - flow1
  mat_2050.RCP85_ <- c(mat_2050.RCP85_)
  
  RCP852050 <- c(RCP852050,sum(mat_2050.RCP85_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))
  
  mat_2050.RCP60_.data <- raster(paste(species_name[i],'_2050.RCP60_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2050.RCP60_.data,mat1.w), filename=paste(species_name[i],'_2050.RCP60','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2050.RCP60_ <- as.matrix(mask(mat_2050.RCP60_.data,mat1.w))
  
  flow_2050.RCP60_ <- rowMeans(mat_2050.RCP60_,na.rm=T) - flow1
  mat_2050.RCP60_ <- c(mat_2050.RCP60_)
  
  RCP602050 <- c(RCP602050,sum(mat_2050.RCP60_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))
  
  mat_2050.RCP45_.data <- raster(paste(species_name[i],'_2050.RCP45_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2050.RCP45_.data,mat1.w), filename=paste(species_name[i],'_2050.RCP45','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2050.RCP45_ <- as.matrix(mask(mat_2050.RCP45_.data,mat1.w))
  
  flow_2050.RCP45_ <- rowMeans(mat_2050.RCP45_,na.rm=T) - flow1
  mat_2050.RCP45_ <- c(mat_2050.RCP45_)
  
  RCP452050 <- c(RCP452050,sum(mat_2050.RCP45_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))
  
  mat_2050.RCP26_.data <- raster(paste(species_name[i],'_2050.RCP26_MaxEnt.tif',sep=''))
  writeRaster(mask(mat_2050.RCP26_.data,mat1.w), filename=paste(species_name[i],'_2050.RCP26','.tif', sep=''), options="INTERLEAVE=BAND", overwrite=TRUE)
  mat_2050.RCP26_ <- as.matrix(mask(mat_2050.RCP26_.data,mat1.w))
  
  flow_2050.RCP26_ <- rowMeans(mat_2050.RCP26_,na.rm=T) - flow1
  mat_2050.RCP26_ <- c(mat_2050.RCP26_)
  
  RCP262050 <- c(RCP262050,sum(mat_2050.RCP26_*A1,na.rm = TRUE)/sum(A1,na.rm = TRUE))
  png(file=paste(species_name[i],"_all_scenarios.png"), width=600, height=350)
  plot(c(2020,2050,2100),c(Present[i],RCP262050[i],RCP262100[i]),
       type="l",col = "blue", lwd=2.0,
       ylab = "Habitat suitability" , xlab = "year", xaxt = "n",
       ylim=c(0,
              max(c(Present[i],RCP262050[i],RCP262100[i],
                    RCP452050[i],RCP452100[i],
                    RCP602050[i],RCP602100[i],
                    RCP852050[i],RCP852100[i]))))
  axis(1, at=c(2020,2050,2100), labels=c("2020","2050","2100"),col.axis = "black")
  lines(c(2020,2050,2100),c(Present[i],RCP452050[i],RCP452100[i]),col = "green", lwd=2.0)
  lines(c(2020,2050,2100),c(Present[i],RCP602050[i],RCP602100[i]),col = "orange", lwd=2.0)
  lines(c(2020,2050,2100),c(Present[i],RCP852050[i],RCP852100[i]),col = "red", lwd=2.0)
  dev.off()
  
  png(file=paste(species_name[i],"_lines_scenarios.png"), width = 5, height = 10, units = 'in', res = 300)
  plot(flow1, -.08333333*(1:length(flow1))+YY,
       type="l",col = "black",
       ylab = "Latitude" , xlab = "Habitat suitability change", lwd=2.0,
       xlim=c(min(min(flow_2100.RCP85_,na.rm=T),min(flow_2100.RCP60_,na.rm=T),min(flow_2100.RCP45_,na.rm=T),min(flow_2100.RCP26_,na.rm=T),
                  min(flow_2050.RCP85_,na.rm=T),min(flow_2050.RCP60_,na.rm=T),min(flow_2050.RCP45_,na.rm=T),min(flow_2050.RCP26_,na.rm=T)),
              max(max(flow1,na.rm=T),max(flow_2100.RCP85_,na.rm=T),max(flow_2100.RCP60_,na.rm=T),max(flow_2100.RCP45_,na.rm=T),max(flow_2100.RCP26_,na.rm=T),
                  max(flow_2050.RCP85_,na.rm=T),max(flow_2050.RCP60_,na.rm=T),max(flow_2050.RCP45_,na.rm=T),max(flow_2050.RCP26_,na.rm=T)))
       )
  lines(flow_2050.RCP26_, -.08333333*(1:length(flow1))+YY,col = "blue", lwd=2.0)
  lines(flow_2050.RCP45_, -.08333333*(1:length(flow1))+YY,col = "green", lwd=2.0)
  lines(flow_2050.RCP60_, -.08333333*(1:length(flow1))+YY,col = "orange", lwd=2.0)
  lines(flow_2050.RCP85_, -.08333333*(1:length(flow1))+YY,col = "red", lwd=2.0)

  lines(flow_2100.RCP26_, -.08333333*(1:length(flow1))+YY,col = "blue", lty=2, lwd=2.0)
  lines(flow_2100.RCP45_, -.08333333*(1:length(flow1))+YY,col = "green", lty=2, lwd=2.0)
  lines(flow_2100.RCP60_, -.08333333*(1:length(flow1))+YY,col = "orange", lty=2, lwd=2.0)
  lines(flow_2100.RCP85_, -.08333333*(1:length(flow1))+YY,col = "red", lty=2, lwd=2.0)
  abline(v=0, col="grey", lwd=1.0)

  dev.off()
}
dat <- data.frame(nameSpecies,
                  Present,
                  RCP852050,
                  RCP602050,
                  RCP452050,
                  RCP262050,
                  RCP852100,
                  RCP602100,
                  RCP452100,
                  RCP262100)
write.csv(dat,"table1.csv", row.names = FALSE)









plotting <- function(points){
  m <-  leaflet()%>%
    addProviderTiles(provider = 'Esri.OceanBasemap') %>%
    addScaleBar()%>%
    addCircles(points$longitude,points$latitude)
  return(m)
}

species_locations <- function(locations,species_name){
  cranch_locations <- select(locations[locations$species==species_name,],
                             longitude,
                             latitude)
  return(cranch_locations)
}

oegopsids <- read.csv("occorrences.csv", na.strings = c("null", "NA"))

for (i in 1:length(species_name)){
  locations <- species_locations(oegopsids,species_name[i])
  m <-  plotting(locations)
  saveWidget(m, file=paste(species_name[i],'.html', sep=''))
}
