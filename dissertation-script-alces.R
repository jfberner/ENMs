# Sample script to allow reproduction of my dissertation data
# Chapter 1 script
# jfb
# 19 jan 2022

# For your convenience this script has been outlined.
# This script was written for a single species, H.alces, and not for the nine species in the study 
# This is an example, and some of the plot and analyses require the data results of different species, not written here
# These not runnable codes are presented, and we explain how to build the data to fit them
# This was a conscious choice

#rm(list = ls()) # clear the environment
path <- "~/Documents/Mestrado/" # set your main directory
setwd(path)
# Libraries ####
library(raster)
library(dplyr)
library(psych)
library(rgeos)
library(rgdal)
library(mapview)
library(sp)
library(sf)
library(maptools)
library(sdm)
library(tidyverse)
library(shapefiles)
library(mraster)
library(CoordinateCleaner)
library(tmap)
library(rworldmap)
library(corrplot)
library(RColorBrewer)
# Load and Prepare Data ####
## Occurrence Data ####
setwd('Occ_Data/Datasets') # path to where your occurrence data are
occ<-read.csv('occ_dataset.csv') # read occurrence data with collumns in DwC (GBIF) standard
occ <- filter(occ,specificEpithet == "alces") # this script will be written only for the species Heterophrynus alces, for simplicity. Comment this line out or include the name of other species if needed.
occ <- filter(occ,coordinateUncertaintyInMeters <=10000) # filtering occurrence data to 10km - the same resolution as our Environmental Data
occ <- as_tibble(occ) # transform to tibble
occ <- dplyr::select(occ,specificEpithet,decimalLatitude,decimalLongitude,coordinateUncertaintyInMeters) # Selecting only needed information
names(occ) <- c("species","lat","lon","uncertainty") # Renaming collumns for simplicity
sp <- as_tibble(occ) # we'll manipulate this object from now on, and leave the raw occ object in the environment in case we need it later

### Cleaning the coordinates #####
flags_spatial <- CoordinateCleaner::clean_coordinates(
  x = sp, 
  species = "species",
  lon = "lon", 
  lat = "lat",
  tests = c("capitals", # radius around capitals
            "centroids", # radius around countries and provinces centroids
            "duplicates", # duplicates
            "equal", # equal coordinates
            "gbif", # radius around GBIF headquarters
            #"institutions", # radius around biodiversity research institutions
            "seas", # points on the ocean
            #"urban", # points inside urban areas
            "validity", # points outside the coordinate system
            "zeros" # zeros e points where lat = lon
  )
)

# the result flags problematic points
#' TRUE = "clean" coordinates
#' FALSE = potentially problematic coordinates

flags_spatial %>% head
summary(flags_spatial)

# removing flagged coords
occ_data_tax_date_spa <- sp %>% 
  dplyr::filter(flags_spatial$.summary == TRUE)
occ_data_tax_date_spa

# exporting the clean coordinates as a csv
setwd(path) ; dir.create('Script_output_files') ; setwd('Script_output_files') ; dir.create('Occ') ; setwd('Occ')
readr::write_csv(occ_data_tax_date_spa, paste0("alces.csv"))

### Making SpatialPoints object and .shp #####

# Ok. Now, before we continue, we need to transform our sp dataframe into
# a Spatial Points Dataframe. To do so, we need to tell R in which rows
# are our coordinates data:

#sp <- read.csv("alces.csv") # uncomment if you're running this part of the script separately, otherwise

coordinates(sp) <- ~lon + lat     #same syntax as columns (see line 34)
proj4string(sp) <- projection(raster()) #assigning WGS84 projection
# now check class(sp), should read SpatialPointsDataFrame

#Now, we have presence only data and the species data has been loaded and cleaned.
writeOGR(obj=sp, dsn="alces.shp", layer=c("species","uncertainty"), driver="ESRI Shapefile") # to create a shapefile with cleaned occurrences

## Environmental Data ####
setwd(path)
setwd('Analises/Chapter1/0-RawEnvlayers') #path to your environmenatal datasets

setwd("BIOLCIM_wc2.1_5m_bio_download19aug2021") #bioclim layers folder
bioclimRaw <- dir(pattern = "tif$")
bioclimRaw <- raster::stack(bioclimRaw)

setwd('..') #one folder up
setwd("ENVIREM_NewWorld_current_5arcmin_generic_download19aug2021") #envirem layers folder
enviremRaw <- dir(pattern = "bil$")
enviremRaw <- raster::stack(enviremRaw)

setwd('..')
setwd("MERRA_5m_mean_00s_download19aug2021") #merraclim layers folder
merraRaw <- dir(pattern = "tif$")
merraRaw <- raster::stack(merraRaw)

# because MERRA resolution is more accurate (0.08333333) than Bioclim and ENVIREM (0.08331), we need to resample it to match the latter ones in order to compare the three

# But this will take a lot of machine time, so first, we'll crop the objects by the extents of South America, which is the region we're interested in:
plot(bioclimRaw[[1]])
extent <- drawExtent() # here I cropped to south america and these are the results:
#class: Extent xmin: -87.8773 xmax: -32.56272 ymin: -126.0932 ymax: 40.18543
# this process does not need to be perfect, we're just interested in reducing the layers for faster computation

bioclimsa <- raster::crop(bioclimRaw,extent)
enviremsa <- raster::crop(enviremRaw,extent)
merrasa <- raster::crop(merraRaw,extent)

setwd(path) ; setwd('Script_output_files') ; dir.create('Env') ; setwd('Env')

enviremsa <- raster::resample(enviremsa,bioclimsa)
merrasa <- raster::resample(merrasa,bioclimsa) # to solve the problem mentioned on line 108

### Correlation ####
#### Bioclim ####
dir.create('Bioclim') ; setwd('Bioclim')

# extracting the values
bioclim_da <- bioclimsa %>%   #bioclim_da = adjusted dimension
  raster::values() %>% 
  tibble::as_tibble() %>% 
  tidyr::drop_na()

# Correlation using Spearman method
cor_table_bioclim <- corrr::correlate(bioclim_da, method = "spearman") 

# Creating a table with the results
cor_table_bioclim_summary <- cor_table_bioclim %>% 
  corrr::shave() %>%
  corrr::fashion()

# And exporting it 
readr::write_csv(cor_table_bioclim_summary, "correlacao.csv")

#### Envirem ####
setwd('..')
dir.create('Envirem') ; setwd('Envirem')
# extracting the values
envirem_da <- enviremsa %>%   #envirem_da = adjusted dimension
  raster::values() %>% 
  tibble::as_tibble() %>% 
  tidyr::drop_na()

# Correlation using Spearman method
cor_table_envirem <- corrr::correlate(envirem_da, method = "spearman") 

# Creating a table with the results
cor_table_envirem_summary <- cor_table_envirem %>% 
  corrr::shave() %>%
  corrr::fashion()

# And exporting it 
readr::write_csv(cor_table_envirem_summary, "correlacao.csv")

#### Merraclim ####
setwd('..')
dir.create('Merra') ; setwd('Merra')
# extracting the values
envirem_da <- enviremsa %>%   #envirem_da = adjusted dimension
  raster::values() %>% 
  tibble::as_tibble() %>% 
  tidyr::drop_na()

# Correlation using Spearman method
cor_table_envirem <- corrr::correlate(envirem_da, method = "spearman") 

# Creating a table with the results
cor_table_envirem_summary <- cor_table_envirem %>% 
  corrr::shave() %>%
  corrr::fashion()

# And exporting it 
readr::write_csv(cor_table_envirem_summary, "correlacao.csv")

#### Correlation between datasets' variables ####
setwd('..')
# to check correlation between the dataset variables, I'll first select only the variables I chose based on the correlation done before, these are: bio2,3,5 and 15 for Bioclim, bio2,3,5 and 8 for Merraclim and annualPET,aridityIndex,climaticMoistureIndex and thermicityIndex for ENVIREM
bioclimsa@data@names # to see the positions of the layers I need to drop
bioclimsa <- dropLayer(bioclimsa,c(1,2,3,4,5,6,8,9,10,11,14,16,17,18,19))

enviremsa@data@names # to see the positions of the layers I need to drop
enviremsa <- dropLayer(enviremsa,c(4,5,6,7,8,9,10,11,12,13,14,15))

merrasa@data@names # to see the positions of the layers I need to drop
merrasa <- dropLayer(merrasa,c(1,2,3,4,5,6,7,8,9,10,11,14,16,17,19))

# So we don't confuse bioclim and merra layers:
names(bioclimsa) <- c('bio15-bioclim','bio2-bioclim','bio3-bioclim','bio5-bioclim')
names(merrasa) <- c('bio2-merra','bio3-merra','bio5-merra','bio8-merra')
names(enviremsa) <- c('annualPET-envirem','aridityIndex-envirem','climaticMoisture-envirem','thermicityIndex-envirem')

env <- raster::stack(bioclimsa,enviremsa,merrasa) #put them all in a single object and run the same calculations as the previous sections

env_da <- env %>%   #env_da = adjusted dimension
  raster::values() %>% 
  tibble::as_tibble() %>% 
  tidyr::drop_na()

# Correlation using Spearman method
cor_table_env <- corrr::correlate(env_da, method = "spearman") 
cor_table_env

# Creating a table with the results
cor_table_env_summary <- cor_table_env %>% 
  corrr::shave() %>%
  corrr::fashion()

# And exporting it 
readr::write_csv(cor_table_env_summary, "correlacao.csv")

# Selecting correlate variables
fi_75_env <- cor_table_env %>% 
  corrr::as_matrix() %>% 
  caret::findCorrelation(cutoff = .75, names = TRUE, verbose = TRUE)  #chose 75% correlation as a threshold

fi_75_env         #listing the correlate variables, roughness does not correlate with anything which I didn't expect

env_da_cor75 <- env_da %>% ###excluding correlate variables
  dplyr::select(-fi_75_env)
env_da_cor75

# Checking the selected variables
env_da_cor75 %>% 
  corrr::correlate(method = "spearman") %>% 
  corrr::as_matrix() %>% 
  caret::findCorrelation(cutoff = .75, names = TRUE, verbose = TRUE)

# Correlation Plot
tiff("correlacao_plot.tiff", wi = 30, he = 25, un = "cm", res = 300, comp = "lzw")
pairs.panels(x = env_da_cor75 %>% dplyr::sample_n(1e3), 
             method = "spearman",
             pch = 20, 
             ellipses = FALSE, 
             density = FALSE, 
             stars = TRUE, 
             hist.col = "gray",
             digits = 2,
             rug = FALSE,
             breaks = 10,
             ci = TRUE)
dev.off()

### Save Cropped, Selected (uncorrelated) Layers ####
# One final step: saving the selected, cropped layers on disk
setwd('Bioclim')
raster::writeRaster(x = bioclimsa, 
                    filename = paste0("", names(bioclimsa)), 
                    bylayer = TRUE, 
                    options = c("COMPRESS=DEFLATE"), 
                    format = "GTiff", 
                    overwrite = TRUE)
setwd('..')
setwd("Merra")
raster::writeRaster(x = merrasa, 
                    filename = paste0("", names(merrasa)), 
                    bylayer = TRUE, 
                    options = c("COMPRESS=DEFLATE"), 
                    format = "GTiff", 
                    overwrite = TRUE)

setwd('..')
setwd("Envirem")
raster::writeRaster(x = enviremsa, 
                    filename = paste0("", names(enviremsa)), 
                    bylayer = TRUE, 
                    options = c("COMPRESS=DEFLATE"), 
                    format = "GTiff", 
                    overwrite = TRUE)




## Make Different M sizes for Env Datasets ####

# To make a Small M (sm), a medium M (mm) and a large M (here, bm) model for each dataset and algorithm, we start by defining the M size based on the distribution of occ records, like so:

radiussm <- 0.5 #1=100km
radiusmm <- 1
radiusbm <- 1.5 # bm = big M, same as LM or large model

#generating polygons around points:
occ.buffer.sm <- gBuffer(spgeom = sp, byid = T, # sm = small M
                         width = radiussm, quadsegs = 100, 
                         capStyle = 'ROUND' , joinStyle = 'ROUND') 

occ.buffer.mm <- gBuffer(spgeom = sp, byid = T, # mm = medium M
                         width = radiusmm, quadsegs = 100, 
                         capStyle = 'ROUND' , joinStyle = 'ROUND') 

occ.buffer.bm <- gBuffer(spgeom = sp, byid = T, # bm = big M
                         width = radiusbm, quadsegs = 100, 
                         capStyle = 'ROUND' , joinStyle = 'ROUND') 

#making it a single polygon
occ.buffer.sm <- aggregate(occ.buffer.sm,dissolve=T) 
occ.buffer.mm <- aggregate(occ.buffer.mm,dissolve=T)
occ.buffer.bm <- aggregate(occ.buffer.bm,dissolve=T)

# Now, make the nine env datasets:
bio.sm <- raster::crop(x = bioclimsa, y = occ.buffer.sm)
bio.mm <- raster::crop(x = bioclimsa, y = occ.buffer.mm)
bio.bm <- raster::crop(x = bioclimsa, y = occ.buffer.bm)

merra.sm <- raster::crop(x = merrasa, y = occ.buffer.sm)
merra.mm <- raster::crop(x = merrasa, y = occ.buffer.mm)
merra.bm <- raster::crop(x = merrasa, y = occ.buffer.bm)

envirem.sm <- raster::crop(x = enviremsa, y = occ.buffer.sm)
envirem.mm <- raster::crop(x = enviremsa, y = occ.buffer.mm)
envirem.bm <- raster::crop(x = enviremsa, y = occ.buffer.bm)



# Modelling ####
setwd('..') ; setwd('..')
dir.create('M-preds') ; setwd('M-preds')

## Creating the Models + Predictions ####
## Ok now we have 3 clim datasets, with three different sizes, the next step is to apply SDMs to the three and compare them with Jimenez&Soberón's functions accum.occ and comp.accum (later section)
### Bioclim ####
# We need to remove all other collumns in the sp SpatialPointsDataFrame so the sdmData function can work
sp@data <- sp@data[1] #to remove uncertainty column

# Prepare Data for the sdm function:
d.bio.sm <- sdmData(formula = species~. , train = sp, predictors = bio.sm,
                    bg = list(n=200), method = 'gRandom')
d.bio.mm <- sdmData(formula = species~. , train = sp, predictors = bio.mm,
                    bg = list(n=200), method = 'gRandom')
d.bio.bm <- sdmData(formula = species~. , train = sp, predictors = bio.bm,
                    bg = list(n=200), method = 'gRandom')
# model
m.bio.sm <- sdm(formula = alces~., 
                data = d.bio.sm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                replication = c('boot') , n = 10) 
m.bio.mm <- sdm(formula = alces~., 
                data = d.bio.mm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                replication = c('boot') , n = 10) 
m.bio.bm <- sdm(formula = alces~., 
                data = d.bio.bm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                replication = c('boot') , n = 10) 
# prediction
p.bio.sm <- predict(m.bio.sm, bio.sm, 'predictions.bio.sm.img', mean = T, overwrite = T)
p.bio.mm <- predict(m.bio.mm, bio.mm, 'predictions.bio.mm.img', mean = T, overwrite = T)
p.bio.bm <- predict(m.bio.bm, bio.bm, 'predictions.bio.bm.img', mean = T, overwrite = T)
### Merraclim #####
## Same as above
# data
d.merra.sm <- sdmData(formula = species~. , train = sp, predictors = merra.sm,
                      bg = list(n=200), method = 'gRandom')
d.merra.mm <- sdmData(formula = species~. , train = sp, predictors = merra.mm,
                      bg = list(n=200), method = 'gRandom')
d.merra.bm <- sdmData(formula = species~. , train = sp, predictors = merra.bm,
                      bg = list(n=200), method = 'gRandom')
# model
m.merra.sm <- sdm(formula = alces~., 
                  data = d.merra.sm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                  replication = c('boot') , n = 10) 
m.merra.mm <- sdm(formula = alces~., 
                  data = d.merra.mm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                  replication = c('boot') , n = 10) 
m.merra.bm <- sdm(formula = alces~., 
                  data = d.merra.bm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                  replication = c('boot') , n = 10) 
# prediction
p.merra.sm <- predict(m.merra.sm, merra.sm, 'predictions.merra.sm.img', mean = T, overwrite = T)
p.merra.mm <- predict(m.merra.mm, merra.mm, 'predictions.merra.mm.img', mean = T, overwrite = T)
p.merra.bm <- predict(m.merra.bm, merra.bm, 'predictions.merra.bm.img', mean = T, overwrite = T)
### ENVIREM #####
# data
d.envirem.sm <- sdmData(formula = species~. , train = sp, predictors = envirem.sm,
                        bg = list(n=200), method = 'gRandom')
d.envirem.mm <- sdmData(formula = species~. , train = sp, predictors = envirem.mm,
                        bg = list(n=200), method = 'gRandom')
d.envirem.bm <- sdmData(formula = species~. , train = sp, predictors = envirem.bm,
                        bg = list(n=200), method = 'gRandom')
# model
m.envirem.sm <- sdm(formula = alces~., 
                    data = d.envirem.sm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                    replication = c('boot') , n = 10) 
m.envirem.mm <- sdm(formula = alces~., 
                    data = d.envirem.mm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                    replication = c('boot') , n = 10) 
m.envirem.bm <- sdm(formula = alces~., 
                    data = d.envirem.bm, methods = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'),
                    replication = c('boot') , n = 10) 
# prediction
p.envirem.sm <- predict(m.envirem.sm, envirem.sm, 'predictions.envirem.sm.img', mean = T, overwrite = T)
p.envirem.mm <- predict(m.envirem.mm, envirem.mm, 'predictions.envirem.mm.img', mean = T, overwrite = T)
p.envirem.bm <- predict(m.envirem.bm, envirem.bm, 'predictions.envirem.bm.img', mean = T, overwrite = T)
## Evaluating the Models - AOcCs ####
# AUC, TSS and other statistics related to ENMs are all available in the M objects, easily accessible with the gui(m.bio.sm) function (for example)

# However, we want to analize the Accumulation of Occurrences Curve (AOcC), proposed by Jimenez & Soberon in 2020. <https://doi.org/10.1111/2041-210X.13479>

# For that, I modified the accum.occ function slightly to show method names on the plot titles. My version can be found at <https://github.com/jfberner/ENMs/blob/main/accum_curve_jfb.RData>

load("~/Documents/Mestrado/Analises/SDM-hyperTest-master/SDM-Hyper-Test-Master-Functions-Jul212021.RData") #original function from Jimenez's github
load("~/Documents/Mestrado/Analises/Model Comparison/accum_curve_jfb.RData") #my modified Jimenez&Soberon function with plot titles, loaded on top of the former to mask the accum.occ function

# We have to manipulate the data a bit to use the accum.occ and comp.accum functions from Laura Jimenez and Jorge Soberon

### Bioclim ####
#### sm ####
bd.bio.sm <- as.data.frame(bio.sm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.bio.sm <- as.data.frame(p.bio.sm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.bio.sm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.bio.sm.glm<-dplyr::select(pd.bio.sm,long,lat,GLM)
output.mod.bio.sm.glm<-cbind(pd.bio.sm.glm,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.glm <- raster::extract(p.bio.sm[[1]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.glm <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.glm,occ.p_env.bio.sm)
occ.pnts.bio.sm.glm <- as.data.frame(occ.pnts.bio.sm.glm)
occ.pnts.bio.sm.glm <-drop_na(occ.pnts.bio.sm.glm)


acc.bio.sm.glm<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.sm.glm,
                          occ.pnts = occ.pnts.bio.sm.glm,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.bio.sm.svm<-dplyr::select(pd.bio.sm,long,lat,SVM)
output.mod.bio.sm.svm<-cbind(pd.bio.sm.svm,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.svm <- raster::extract(p.bio.sm[[2]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.svm <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.svm,occ.p_env.bio.sm)
occ.pnts.bio.sm.svm <- as.data.frame(occ.pnts.bio.sm.svm)
occ.pnts.bio.sm.svm <-drop_na(occ.pnts.bio.sm.svm)


acc.bio.sm.svm<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.sm.svm,
                          occ.pnts = occ.pnts.bio.sm.svm,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.bio.sm.rf<-dplyr::select(pd.bio.sm,long,lat,RF)
output.mod.bio.sm.rf<-cbind(pd.bio.sm.rf,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.rf <- raster::extract(p.bio.sm[[3]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.rf <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.rf,occ.p_env.bio.sm)
occ.pnts.bio.sm.rf <- as.data.frame(occ.pnts.bio.sm.rf)
occ.pnts.bio.sm.rf <-drop_na(occ.pnts.bio.sm.rf)


acc.bio.sm.rf<-accum.occ(sp.name='Heterophrynus',
                         output.mod = output.mod.bio.sm.rf,
                         occ.pnts = occ.pnts.bio.sm.rf,
                         null.mod = "hypergeom",
                         conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.bio.sm.brt<-dplyr::select(pd.bio.sm,long,lat,BRT)
output.mod.bio.sm.brt<-cbind(pd.bio.sm.brt,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.brt <- raster::extract(p.bio.sm[[4]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.brt <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.brt,occ.p_env.bio.sm)
occ.pnts.bio.sm.brt <- as.data.frame(occ.pnts.bio.sm.brt)
occ.pnts.bio.sm.brt <-drop_na(occ.pnts.bio.sm.brt)


acc.bio.sm.brt<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.sm.brt,
                          occ.pnts = occ.pnts.bio.sm.brt,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.bio.sm.mars<-dplyr::select(pd.bio.sm,long,lat,MARS)
output.mod.bio.sm.mars<-cbind(pd.bio.sm.mars,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.mars <- raster::extract(p.bio.sm[[5]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.mars <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.mars,occ.p_env.bio.sm)
occ.pnts.bio.sm.mars <- as.data.frame(occ.pnts.bio.sm.mars)
occ.pnts.bio.sm.mars <-drop_na(occ.pnts.bio.sm.mars)


acc.bio.sm.mars<-accum.occ(sp.name='Heterophrynus',
                           output.mod = output.mod.bio.sm.mars,
                           occ.pnts = occ.pnts.bio.sm.mars,
                           null.mod = "hypergeom",
                           conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.bio.sm.maxent<-dplyr::select(pd.bio.sm,long,lat,MAXENT)
output.mod.bio.sm.maxent<-cbind(pd.bio.sm.maxent,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.maxent <- raster::extract(p.bio.sm[[6]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.maxent <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.maxent,occ.p_env.bio.sm)
occ.pnts.bio.sm.maxent <- as.data.frame(occ.pnts.bio.sm.maxent)
occ.pnts.bio.sm.maxent <-drop_na(occ.pnts.bio.sm.maxent)


acc.bio.sm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                             output.mod = output.mod.bio.sm.maxent,
                             occ.pnts = occ.pnts.bio.sm.maxent,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.bio.sm.maxlike<-dplyr::select(pd.bio.sm,long,lat,MAXLIKE)
output.mod.bio.sm.maxlike<-cbind(pd.bio.sm.maxlike,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.maxlike <- raster::extract(p.bio.sm[[7]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.maxlike <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.maxlike,occ.p_env.bio.sm)
occ.pnts.bio.sm.maxlike <- as.data.frame(occ.pnts.bio.sm.maxlike)
occ.pnts.bio.sm.maxlike <-drop_na(occ.pnts.bio.sm.maxlike)


acc.bio.sm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                              output.mod = output.mod.bio.sm.maxlike,
                              occ.pnts = occ.pnts.bio.sm.maxlike,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.bio.sm.glmnet<-dplyr::select(pd.bio.sm,long,lat,GLMNET)
output.mod.bio.sm.glmnet<-cbind(pd.bio.sm.glmnet,bd.bio.sm)

occ.p_env.bio.sm <- raster::extract(bio.sm,sp,as.data.frame=T)
occ.p_preds.bio.sm.glmnet <- raster::extract(p.bio.sm[[8]],sp,as.data.frame=T)
occ.coords.bio.sm <- sp@coords
occ.pnts.bio.sm.glmnet <- cbind(occ.coords.bio.sm,occ.p_preds.bio.sm.glmnet,occ.p_env.bio.sm)
occ.pnts.bio.sm.glmnet <- as.data.frame(occ.pnts.bio.sm.glmnet)
occ.pnts.bio.sm.glmnet <-drop_na(occ.pnts.bio.sm.glmnet)


acc.bio.sm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                             output.mod = output.mod.bio.sm.glmnet,
                             occ.pnts = occ.pnts.bio.sm.glmnet,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()
#### mm ####
bd.bio.mm <- as.data.frame(bio.mm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.bio.mm <- as.data.frame(p.bio.mm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.bio.mm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.bio.mm.glm<-dplyr::select(pd.bio.mm,long,lat,GLM)
output.mod.bio.mm.glm<-cbind(pd.bio.mm.glm,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.glm <- raster::extract(p.bio.mm[[1]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.glm <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.glm,occ.p_env.bio.mm)
occ.pnts.bio.mm.glm <- as.data.frame(occ.pnts.bio.mm.glm)
occ.pnts.bio.mm.glm <-drop_na(occ.pnts.bio.mm.glm)


acc.bio.mm.glm<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.mm.glm,
                          occ.pnts = occ.pnts.bio.mm.glm,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.bio.mm.svm<-dplyr::select(pd.bio.mm,long,lat,SVM)
output.mod.bio.mm.svm<-cbind(pd.bio.mm.svm,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.svm <- raster::extract(p.bio.mm[[2]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.svm <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.svm,occ.p_env.bio.mm)
occ.pnts.bio.mm.svm <- as.data.frame(occ.pnts.bio.mm.svm)
occ.pnts.bio.mm.svm <-drop_na(occ.pnts.bio.mm.svm)


acc.bio.mm.svm<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.mm.svm,
                          occ.pnts = occ.pnts.bio.mm.svm,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.bio.mm.rf<-dplyr::select(pd.bio.mm,long,lat,RF)
output.mod.bio.mm.rf<-cbind(pd.bio.mm.rf,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.rf <- raster::extract(p.bio.mm[[3]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.rf <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.rf,occ.p_env.bio.mm)
occ.pnts.bio.mm.rf <- as.data.frame(occ.pnts.bio.mm.rf)
occ.pnts.bio.mm.rf <-drop_na(occ.pnts.bio.mm.rf)


acc.bio.mm.rf<-accum.occ(sp.name='Heterophrynus',
                         output.mod = output.mod.bio.mm.rf,
                         occ.pnts = occ.pnts.bio.mm.rf,
                         null.mod = "hypergeom",
                         conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.bio.mm.brt<-dplyr::select(pd.bio.mm,long,lat,BRT)
output.mod.bio.mm.brt<-cbind(pd.bio.mm.brt,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.brt <- raster::extract(p.bio.mm[[4]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.brt <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.brt,occ.p_env.bio.mm)
occ.pnts.bio.mm.brt <- as.data.frame(occ.pnts.bio.mm.brt)
occ.pnts.bio.mm.brt <-drop_na(occ.pnts.bio.mm.brt)


acc.bio.mm.brt<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.mm.brt,
                          occ.pnts = occ.pnts.bio.mm.brt,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.bio.mm.mars<-dplyr::select(pd.bio.mm,long,lat,MARS)
output.mod.bio.mm.mars<-cbind(pd.bio.mm.mars,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.mars <- raster::extract(p.bio.mm[[5]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.mars <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.mars,occ.p_env.bio.mm)
occ.pnts.bio.mm.mars <- as.data.frame(occ.pnts.bio.mm.mars)
occ.pnts.bio.mm.mars <-drop_na(occ.pnts.bio.mm.mars)


acc.bio.mm.mars<-accum.occ(sp.name='Heterophrynus',
                           output.mod = output.mod.bio.mm.mars,
                           occ.pnts = occ.pnts.bio.mm.mars,
                           null.mod = "hypergeom",
                           conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.bio.mm.maxent<-dplyr::select(pd.bio.mm,long,lat,MAXENT)
output.mod.bio.mm.maxent<-cbind(pd.bio.mm.maxent,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.maxent <- raster::extract(p.bio.mm[[6]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.maxent <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.maxent,occ.p_env.bio.mm)
occ.pnts.bio.mm.maxent <- as.data.frame(occ.pnts.bio.mm.maxent)
occ.pnts.bio.mm.maxent <-drop_na(occ.pnts.bio.mm.maxent)


acc.bio.mm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                             output.mod = output.mod.bio.mm.maxent,
                             occ.pnts = occ.pnts.bio.mm.maxent,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.bio.mm.maxlike<-dplyr::select(pd.bio.mm,long,lat,MAXLIKE)
output.mod.bio.mm.maxlike<-cbind(pd.bio.mm.maxlike,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.maxlike <- raster::extract(p.bio.mm[[7]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.maxlike <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.maxlike,occ.p_env.bio.mm)
occ.pnts.bio.mm.maxlike <- as.data.frame(occ.pnts.bio.mm.maxlike)
occ.pnts.bio.mm.maxlike <-drop_na(occ.pnts.bio.mm.maxlike)


acc.bio.mm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                              output.mod = output.mod.bio.mm.maxlike,
                              occ.pnts = occ.pnts.bio.mm.maxlike,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.bio.mm.glmnet<-dplyr::select(pd.bio.mm,long,lat,GLMNET)
output.mod.bio.mm.glmnet<-cbind(pd.bio.mm.glmnet,bd.bio.mm)

occ.p_env.bio.mm <- raster::extract(bio.mm,sp,as.data.frame=T)
occ.p_preds.bio.mm.glmnet <- raster::extract(p.bio.mm[[8]],sp,as.data.frame=T)
occ.coords.bio.mm <- sp@coords
occ.pnts.bio.mm.glmnet <- cbind(occ.coords.bio.mm,occ.p_preds.bio.mm.glmnet,occ.p_env.bio.mm)
occ.pnts.bio.mm.glmnet <- as.data.frame(occ.pnts.bio.mm.glmnet)
occ.pnts.bio.mm.glmnet <-drop_na(occ.pnts.bio.mm.glmnet)


acc.bio.mm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                             output.mod = output.mod.bio.mm.glmnet,
                             occ.pnts = occ.pnts.bio.mm.glmnet,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()
#### bm ####
bd.bio.bm <- as.data.frame(bio.bm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.bio.bm <- as.data.frame(p.bio.bm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.bio.bm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.bio.bm.glm<-dplyr::select(pd.bio.bm,long,lat,GLM)
output.mod.bio.bm.glm<-cbind(pd.bio.bm.glm,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.glm <- raster::extract(p.bio.bm[[1]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.glm <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.glm,occ.p_env.bio.bm)
occ.pnts.bio.bm.glm <- as.data.frame(occ.pnts.bio.bm.glm)
occ.pnts.bio.bm.glm <-drop_na(occ.pnts.bio.bm.glm)


acc.bio.bm.glm<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.bm.glm,
                          occ.pnts = occ.pnts.bio.bm.glm,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.bio.bm.svm<-dplyr::select(pd.bio.bm,long,lat,SVM)
output.mod.bio.bm.svm<-cbind(pd.bio.bm.svm,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.svm <- raster::extract(p.bio.bm[[2]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.svm <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.svm,occ.p_env.bio.bm)
occ.pnts.bio.bm.svm <- as.data.frame(occ.pnts.bio.bm.svm)
occ.pnts.bio.bm.svm <-drop_na(occ.pnts.bio.bm.svm)


acc.bio.bm.svm<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.bm.svm,
                          occ.pnts = occ.pnts.bio.bm.svm,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.bio.bm.rf<-dplyr::select(pd.bio.bm,long,lat,RF)
output.mod.bio.bm.rf<-cbind(pd.bio.bm.rf,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.rf <- raster::extract(p.bio.bm[[3]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.rf <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.rf,occ.p_env.bio.bm)
occ.pnts.bio.bm.rf <- as.data.frame(occ.pnts.bio.bm.rf)
occ.pnts.bio.bm.rf <-drop_na(occ.pnts.bio.bm.rf)


acc.bio.bm.rf<-accum.occ(sp.name='Heterophrynus',
                         output.mod = output.mod.bio.bm.rf,
                         occ.pnts = occ.pnts.bio.bm.rf,
                         null.mod = "hypergeom",
                         conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.bio.bm.brt<-dplyr::select(pd.bio.bm,long,lat,BRT)
output.mod.bio.bm.brt<-cbind(pd.bio.bm.brt,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.brt <- raster::extract(p.bio.bm[[4]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.brt <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.brt,occ.p_env.bio.bm)
occ.pnts.bio.bm.brt <- as.data.frame(occ.pnts.bio.bm.brt)
occ.pnts.bio.bm.brt <-drop_na(occ.pnts.bio.bm.brt)


acc.bio.bm.brt<-accum.occ(sp.name='Heterophrynus',
                          output.mod = output.mod.bio.bm.brt,
                          occ.pnts = occ.pnts.bio.bm.brt,
                          null.mod = "hypergeom",
                          conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.bio.bm.mars<-dplyr::select(pd.bio.bm,long,lat,MARS)
output.mod.bio.bm.mars<-cbind(pd.bio.bm.mars,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.mars <- raster::extract(p.bio.bm[[5]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.mars <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.mars,occ.p_env.bio.bm)
occ.pnts.bio.bm.mars <- as.data.frame(occ.pnts.bio.bm.mars)
occ.pnts.bio.bm.mars <-drop_na(occ.pnts.bio.bm.mars)


acc.bio.bm.mars<-accum.occ(sp.name='Heterophrynus',
                           output.mod = output.mod.bio.bm.mars,
                           occ.pnts = occ.pnts.bio.bm.mars,
                           null.mod = "hypergeom",
                           conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.bio.bm.maxent<-dplyr::select(pd.bio.bm,long,lat,MAXENT)
output.mod.bio.bm.maxent<-cbind(pd.bio.bm.maxent,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.maxent <- raster::extract(p.bio.bm[[6]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.maxent <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.maxent,occ.p_env.bio.bm)
occ.pnts.bio.bm.maxent <- as.data.frame(occ.pnts.bio.bm.maxent)
occ.pnts.bio.bm.maxent <-drop_na(occ.pnts.bio.bm.maxent)


acc.bio.bm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                             output.mod = output.mod.bio.bm.maxent,
                             occ.pnts = occ.pnts.bio.bm.maxent,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.bio.bm.maxlike<-dplyr::select(pd.bio.bm,long,lat,MAXLIKE)
output.mod.bio.bm.maxlike<-cbind(pd.bio.bm.maxlike,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.maxlike <- raster::extract(p.bio.bm[[7]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.maxlike <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.maxlike,occ.p_env.bio.bm)
occ.pnts.bio.bm.maxlike <- as.data.frame(occ.pnts.bio.bm.maxlike)
occ.pnts.bio.bm.maxlike <-drop_na(occ.pnts.bio.bm.maxlike)


acc.bio.bm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                              output.mod = output.mod.bio.bm.maxlike,
                              occ.pnts = occ.pnts.bio.bm.maxlike,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.bio.bm.glmnet<-dplyr::select(pd.bio.bm,long,lat,GLMNET)
output.mod.bio.bm.glmnet<-cbind(pd.bio.bm.glmnet,bd.bio.bm)

occ.p_env.bio.bm <- raster::extract(bio.bm,sp,as.data.frame=T)
occ.p_preds.bio.bm.glmnet <- raster::extract(p.bio.bm[[8]],sp,as.data.frame=T)
occ.coords.bio.bm <- sp@coords
occ.pnts.bio.bm.glmnet <- cbind(occ.coords.bio.bm,occ.p_preds.bio.bm.glmnet,occ.p_env.bio.bm)
occ.pnts.bio.bm.glmnet <- as.data.frame(occ.pnts.bio.bm.glmnet)
occ.pnts.bio.bm.glmnet <-drop_na(occ.pnts.bio.bm.glmnet)


acc.bio.bm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                             output.mod = output.mod.bio.bm.glmnet,
                             occ.pnts = occ.pnts.bio.bm.glmnet,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#### Generating the AOcCs - comp.accplot ####
# Load the function from Jimenez&Soberon
load("~/Documents/Mestrado/Analises/SDM-hyperTest-master/SDM-Hyper-Test-Master-Functions-Jul212021.RData") 
load("~/Documents/Mestrado/Analises/Model Comparison/accum_curve_jfb.RData") #my modified Jimenez&Soberon function with plot titles, load it again to overwrite the accum.occ function

# sm 

bio.sm.comp <- list(acc.bio.sm.glm,acc.bio.sm.svm,acc.bio.sm.rf,acc.bio.sm.brt,acc.bio.sm.mars,acc.bio.sm.maxent,acc.bio.sm.maxlike,acc.bio.sm.glmnet)

comp.accplot(mods=bio.sm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.bio.sm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Small Model Comparison")
dev.off()
# mm
bio.mm.comp <- list(acc.bio.mm.glm,acc.bio.mm.svm,acc.bio.mm.rf,acc.bio.mm.brt,acc.bio.mm.mars,acc.bio.mm.maxent,acc.bio.mm.maxlike,acc.bio.mm.glmnet)

comp.accplot(mods=bio.mm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.bio.mm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Medium Model Comparison")
dev.off()
# bm
bio.bm.comp <- list(acc.bio.bm.glm,acc.bio.bm.svm,acc.bio.bm.rf,acc.bio.bm.brt,acc.bio.bm.mars,acc.bio.bm.maxent,acc.bio.bm.maxlike,acc.bio.bm.glmnet)

comp.accplot(mods=bio.bm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.bio.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Big Model Comparison")
dev.off()
### Merraclim ####
#### sm ####
bd.merra.sm <- as.data.frame(merra.sm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.merra.sm <- as.data.frame(p.merra.sm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.merra.sm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.merra.sm.glm<-dplyr::select(pd.merra.sm,long,lat,GLM)
output.mod.merra.sm.glm<-cbind(pd.merra.sm.glm,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.glm <- raster::extract(p.merra.sm[[1]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.glm <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.glm,occ.p_env.merra.sm)
occ.pnts.merra.sm.glm <- as.data.frame(occ.pnts.merra.sm.glm)
occ.pnts.merra.sm.glm <-drop_na(occ.pnts.merra.sm.glm)


acc.merra.sm.glm<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.sm.glm,
                            occ.pnts = occ.pnts.merra.sm.glm,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.merra.sm.svm<-dplyr::select(pd.merra.sm,long,lat,SVM)
output.mod.merra.sm.svm<-cbind(pd.merra.sm.svm,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.svm <- raster::extract(p.merra.sm[[2]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.svm <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.svm,occ.p_env.merra.sm)
occ.pnts.merra.sm.svm <- as.data.frame(occ.pnts.merra.sm.svm)
occ.pnts.merra.sm.svm <-drop_na(occ.pnts.merra.sm.svm)


acc.merra.sm.svm<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.sm.svm,
                            occ.pnts = occ.pnts.merra.sm.svm,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.merra.sm.rf<-dplyr::select(pd.merra.sm,long,lat,RF)
output.mod.merra.sm.rf<-cbind(pd.merra.sm.rf,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.rf <- raster::extract(p.merra.sm[[3]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.rf <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.rf,occ.p_env.merra.sm)
occ.pnts.merra.sm.rf <- as.data.frame(occ.pnts.merra.sm.rf)
occ.pnts.merra.sm.rf <-drop_na(occ.pnts.merra.sm.rf)


acc.merra.sm.rf<-accum.occ(sp.name='Heterophrynus',
                           output.mod = output.mod.merra.sm.rf,
                           occ.pnts = occ.pnts.merra.sm.rf,
                           null.mod = "hypergeom",
                           conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.merra.sm.brt<-dplyr::select(pd.merra.sm,long,lat,BRT)
output.mod.merra.sm.brt<-cbind(pd.merra.sm.brt,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.brt <- raster::extract(p.merra.sm[[4]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.brt <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.brt,occ.p_env.merra.sm)
occ.pnts.merra.sm.brt <- as.data.frame(occ.pnts.merra.sm.brt)
occ.pnts.merra.sm.brt <-drop_na(occ.pnts.merra.sm.brt)


acc.merra.sm.brt<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.sm.brt,
                            occ.pnts = occ.pnts.merra.sm.brt,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.merra.sm.mars<-dplyr::select(pd.merra.sm,long,lat,MARS)
output.mod.merra.sm.mars<-cbind(pd.merra.sm.mars,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.mars <- raster::extract(p.merra.sm[[5]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.mars <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.mars,occ.p_env.merra.sm)
occ.pnts.merra.sm.mars <- as.data.frame(occ.pnts.merra.sm.mars)
occ.pnts.merra.sm.mars <-drop_na(occ.pnts.merra.sm.mars)


acc.merra.sm.mars<-accum.occ(sp.name='Heterophrynus',
                             output.mod = output.mod.merra.sm.mars,
                             occ.pnts = occ.pnts.merra.sm.mars,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.merra.sm.maxent<-dplyr::select(pd.merra.sm,long,lat,MAXENT)
output.mod.merra.sm.maxent<-cbind(pd.merra.sm.maxent,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.maxent <- raster::extract(p.merra.sm[[6]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.maxent <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.maxent,occ.p_env.merra.sm)
occ.pnts.merra.sm.maxent <- as.data.frame(occ.pnts.merra.sm.maxent)
occ.pnts.merra.sm.maxent <-drop_na(occ.pnts.merra.sm.maxent)


acc.merra.sm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                               output.mod = output.mod.merra.sm.maxent,
                               occ.pnts = occ.pnts.merra.sm.maxent,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.merra.sm.maxlike<-dplyr::select(pd.merra.sm,long,lat,MAXLIKE)
output.mod.merra.sm.maxlike<-cbind(pd.merra.sm.maxlike,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.maxlike <- raster::extract(p.merra.sm[[7]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.maxlike <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.maxlike,occ.p_env.merra.sm)
occ.pnts.merra.sm.maxlike <- as.data.frame(occ.pnts.merra.sm.maxlike)
occ.pnts.merra.sm.maxlike <-drop_na(occ.pnts.merra.sm.maxlike)


acc.merra.sm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                                output.mod = output.mod.merra.sm.maxlike,
                                occ.pnts = occ.pnts.merra.sm.maxlike,
                                null.mod = "hypergeom",
                                conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.merra.sm.glmnet<-dplyr::select(pd.merra.sm,long,lat,GLMNET)
output.mod.merra.sm.glmnet<-cbind(pd.merra.sm.glmnet,bd.merra.sm)

occ.p_env.merra.sm <- raster::extract(merra.sm,sp,as.data.frame=T)
occ.p_preds.merra.sm.glmnet <- raster::extract(p.merra.sm[[8]],sp,as.data.frame=T)
occ.coords.merra.sm <- sp@coords
occ.pnts.merra.sm.glmnet <- cbind(occ.coords.merra.sm,occ.p_preds.merra.sm.glmnet,occ.p_env.merra.sm)
occ.pnts.merra.sm.glmnet <- as.data.frame(occ.pnts.merra.sm.glmnet)
occ.pnts.merra.sm.glmnet <-drop_na(occ.pnts.merra.sm.glmnet)


acc.merra.sm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                               output.mod = output.mod.merra.sm.glmnet,
                               occ.pnts = occ.pnts.merra.sm.glmnet,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()
#### mm ####
bd.merra.mm <- as.data.frame(merra.mm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.merra.mm <- as.data.frame(p.merra.mm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.merra.mm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.merra.mm.glm<-dplyr::select(pd.merra.mm,long,lat,GLM)
output.mod.merra.mm.glm<-cbind(pd.merra.mm.glm,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.glm <- raster::extract(p.merra.mm[[1]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.glm <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.glm,occ.p_env.merra.mm)
occ.pnts.merra.mm.glm <- as.data.frame(occ.pnts.merra.mm.glm)
occ.pnts.merra.mm.glm <-drop_na(occ.pnts.merra.mm.glm)


acc.merra.mm.glm<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.mm.glm,
                            occ.pnts = occ.pnts.merra.mm.glm,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.merra.mm.svm<-dplyr::select(pd.merra.mm,long,lat,SVM)
output.mod.merra.mm.svm<-cbind(pd.merra.mm.svm,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.svm <- raster::extract(p.merra.mm[[2]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.svm <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.svm,occ.p_env.merra.mm)
occ.pnts.merra.mm.svm <- as.data.frame(occ.pnts.merra.mm.svm)
occ.pnts.merra.mm.svm <-drop_na(occ.pnts.merra.mm.svm)


acc.merra.mm.svm<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.mm.svm,
                            occ.pnts = occ.pnts.merra.mm.svm,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.merra.mm.rf<-dplyr::select(pd.merra.mm,long,lat,RF)
output.mod.merra.mm.rf<-cbind(pd.merra.mm.rf,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.rf <- raster::extract(p.merra.mm[[3]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.rf <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.rf,occ.p_env.merra.mm)
occ.pnts.merra.mm.rf <- as.data.frame(occ.pnts.merra.mm.rf)
occ.pnts.merra.mm.rf <-drop_na(occ.pnts.merra.mm.rf)


acc.merra.mm.rf<-accum.occ(sp.name='Heterophrynus',
                           output.mod = output.mod.merra.mm.rf,
                           occ.pnts = occ.pnts.merra.mm.rf,
                           null.mod = "hypergeom",
                           conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.merra.mm.brt<-dplyr::select(pd.merra.mm,long,lat,BRT)
output.mod.merra.mm.brt<-cbind(pd.merra.mm.brt,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.brt <- raster::extract(p.merra.mm[[4]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.brt <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.brt,occ.p_env.merra.mm)
occ.pnts.merra.mm.brt <- as.data.frame(occ.pnts.merra.mm.brt)
occ.pnts.merra.mm.brt <-drop_na(occ.pnts.merra.mm.brt)


acc.merra.mm.brt<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.mm.brt,
                            occ.pnts = occ.pnts.merra.mm.brt,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.merra.mm.mars<-dplyr::select(pd.merra.mm,long,lat,MARS)
output.mod.merra.mm.mars<-cbind(pd.merra.mm.mars,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.mars <- raster::extract(p.merra.mm[[5]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.mars <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.mars,occ.p_env.merra.mm)
occ.pnts.merra.mm.mars <- as.data.frame(occ.pnts.merra.mm.mars)
occ.pnts.merra.mm.mars <-drop_na(occ.pnts.merra.mm.mars)


acc.merra.mm.mars<-accum.occ(sp.name='Heterophrynus',
                             output.mod = output.mod.merra.mm.mars,
                             occ.pnts = occ.pnts.merra.mm.mars,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.merra.mm.maxent<-dplyr::select(pd.merra.mm,long,lat,MAXENT)
output.mod.merra.mm.maxent<-cbind(pd.merra.mm.maxent,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.maxent <- raster::extract(p.merra.mm[[6]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.maxent <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.maxent,occ.p_env.merra.mm)
occ.pnts.merra.mm.maxent <- as.data.frame(occ.pnts.merra.mm.maxent)
occ.pnts.merra.mm.maxent <-drop_na(occ.pnts.merra.mm.maxent)


acc.merra.mm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                               output.mod = output.mod.merra.mm.maxent,
                               occ.pnts = occ.pnts.merra.mm.maxent,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.merra.mm.maxlike<-dplyr::select(pd.merra.mm,long,lat,MAXLIKE)
output.mod.merra.mm.maxlike<-cbind(pd.merra.mm.maxlike,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.maxlike <- raster::extract(p.merra.mm[[7]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.maxlike <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.maxlike,occ.p_env.merra.mm)
occ.pnts.merra.mm.maxlike <- as.data.frame(occ.pnts.merra.mm.maxlike)
occ.pnts.merra.mm.maxlike <-drop_na(occ.pnts.merra.mm.maxlike)


acc.merra.mm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                                output.mod = output.mod.merra.mm.maxlike,
                                occ.pnts = occ.pnts.merra.mm.maxlike,
                                null.mod = "hypergeom",
                                conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.merra.mm.glmnet<-dplyr::select(pd.merra.mm,long,lat,GLMNET)
output.mod.merra.mm.glmnet<-cbind(pd.merra.mm.glmnet,bd.merra.mm)

occ.p_env.merra.mm <- raster::extract(merra.mm,sp,as.data.frame=T)
occ.p_preds.merra.mm.glmnet <- raster::extract(p.merra.mm[[8]],sp,as.data.frame=T)
occ.coords.merra.mm <- sp@coords
occ.pnts.merra.mm.glmnet <- cbind(occ.coords.merra.mm,occ.p_preds.merra.mm.glmnet,occ.p_env.merra.mm)
occ.pnts.merra.mm.glmnet <- as.data.frame(occ.pnts.merra.mm.glmnet)
occ.pnts.merra.mm.glmnet <-drop_na(occ.pnts.merra.mm.glmnet)


acc.merra.mm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                               output.mod = output.mod.merra.mm.glmnet,
                               occ.pnts = occ.pnts.merra.mm.glmnet,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()
#### bm ####
bd.merra.bm <- as.data.frame(merra.bm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.merra.bm <- as.data.frame(p.merra.bm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.merra.bm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.merra.bm.glm<-dplyr::select(pd.merra.bm,long,lat,GLM)
output.mod.merra.bm.glm<-cbind(pd.merra.bm.glm,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.glm <- raster::extract(p.merra.bm[[1]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.glm <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.glm,occ.p_env.merra.bm)
occ.pnts.merra.bm.glm <- as.data.frame(occ.pnts.merra.bm.glm)
occ.pnts.merra.bm.glm <-drop_na(occ.pnts.merra.bm.glm)


acc.merra.bm.glm<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.bm.glm,
                            occ.pnts = occ.pnts.merra.bm.glm,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.merra.bm.svm<-dplyr::select(pd.merra.bm,long,lat,SVM)
output.mod.merra.bm.svm<-cbind(pd.merra.bm.svm,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.svm <- raster::extract(p.merra.bm[[2]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.svm <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.svm,occ.p_env.merra.bm)
occ.pnts.merra.bm.svm <- as.data.frame(occ.pnts.merra.bm.svm)
occ.pnts.merra.bm.svm <-drop_na(occ.pnts.merra.bm.svm)


acc.merra.bm.svm<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.bm.svm,
                            occ.pnts = occ.pnts.merra.bm.svm,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.merra.bm.rf<-dplyr::select(pd.merra.bm,long,lat,RF)
output.mod.merra.bm.rf<-cbind(pd.merra.bm.rf,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.rf <- raster::extract(p.merra.bm[[3]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.rf <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.rf,occ.p_env.merra.bm)
occ.pnts.merra.bm.rf <- as.data.frame(occ.pnts.merra.bm.rf)
occ.pnts.merra.bm.rf <-drop_na(occ.pnts.merra.bm.rf)


acc.merra.bm.rf<-accum.occ(sp.name='Heterophrynus',
                           output.mod = output.mod.merra.bm.rf,
                           occ.pnts = occ.pnts.merra.bm.rf,
                           null.mod = "hypergeom",
                           conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.merra.bm.brt<-dplyr::select(pd.merra.bm,long,lat,BRT)
output.mod.merra.bm.brt<-cbind(pd.merra.bm.brt,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.brt <- raster::extract(p.merra.bm[[4]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.brt <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.brt,occ.p_env.merra.bm)
occ.pnts.merra.bm.brt <- as.data.frame(occ.pnts.merra.bm.brt)
occ.pnts.merra.bm.brt <-drop_na(occ.pnts.merra.bm.brt)


acc.merra.bm.brt<-accum.occ(sp.name='Heterophrynus',
                            output.mod = output.mod.merra.bm.brt,
                            occ.pnts = occ.pnts.merra.bm.brt,
                            null.mod = "hypergeom",
                            conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.merra.bm.mars<-dplyr::select(pd.merra.bm,long,lat,MARS)
output.mod.merra.bm.mars<-cbind(pd.merra.bm.mars,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.mars <- raster::extract(p.merra.bm[[5]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.mars <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.mars,occ.p_env.merra.bm)
occ.pnts.merra.bm.mars <- as.data.frame(occ.pnts.merra.bm.mars)
occ.pnts.merra.bm.mars <-drop_na(occ.pnts.merra.bm.mars)


acc.merra.bm.mars<-accum.occ(sp.name='Heterophrynus',
                             output.mod = output.mod.merra.bm.mars,
                             occ.pnts = occ.pnts.merra.bm.mars,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.merra.bm.maxent<-dplyr::select(pd.merra.bm,long,lat,MAXENT)
output.mod.merra.bm.maxent<-cbind(pd.merra.bm.maxent,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.maxent <- raster::extract(p.merra.bm[[6]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.maxent <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.maxent,occ.p_env.merra.bm)
occ.pnts.merra.bm.maxent <- as.data.frame(occ.pnts.merra.bm.maxent)
occ.pnts.merra.bm.maxent <-drop_na(occ.pnts.merra.bm.maxent)


acc.merra.bm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                               output.mod = output.mod.merra.bm.maxent,
                               occ.pnts = occ.pnts.merra.bm.maxent,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.merra.bm.maxlike<-dplyr::select(pd.merra.bm,long,lat,MAXLIKE)
output.mod.merra.bm.maxlike<-cbind(pd.merra.bm.maxlike,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.maxlike <- raster::extract(p.merra.bm[[7]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.maxlike <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.maxlike,occ.p_env.merra.bm)
occ.pnts.merra.bm.maxlike <- as.data.frame(occ.pnts.merra.bm.maxlike)
occ.pnts.merra.bm.maxlike <-drop_na(occ.pnts.merra.bm.maxlike)


acc.merra.bm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                                output.mod = output.mod.merra.bm.maxlike,
                                occ.pnts = occ.pnts.merra.bm.maxlike,
                                null.mod = "hypergeom",
                                conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.merra.bm.glmnet<-dplyr::select(pd.merra.bm,long,lat,GLMNET)
output.mod.merra.bm.glmnet<-cbind(pd.merra.bm.glmnet,bd.merra.bm)

occ.p_env.merra.bm <- raster::extract(merra.bm,sp,as.data.frame=T)
occ.p_preds.merra.bm.glmnet <- raster::extract(p.merra.bm[[8]],sp,as.data.frame=T)
occ.coords.merra.bm <- sp@coords
occ.pnts.merra.bm.glmnet <- cbind(occ.coords.merra.bm,occ.p_preds.merra.bm.glmnet,occ.p_env.merra.bm)
occ.pnts.merra.bm.glmnet <- as.data.frame(occ.pnts.merra.bm.glmnet)
occ.pnts.merra.bm.glmnet <-drop_na(occ.pnts.merra.bm.glmnet)


acc.merra.bm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                               output.mod = output.mod.merra.bm.glmnet,
                               occ.pnts = occ.pnts.merra.bm.glmnet,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#### Generating the AOcCs - comp.accplot ####
# sm 
merra.sm.comp <- list(acc.merra.sm.glm,acc.merra.sm.svm,acc.merra.sm.rf,acc.merra.sm.brt,acc.merra.sm.mars,acc.merra.sm.maxent,acc.merra.sm.maxlike,acc.merra.sm.glmnet)

comp.accplot(mods=merra.sm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.merra.sm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Small Model Comparison")
dev.off()

# mm
merra.mm.comp <- list(acc.merra.mm.glm,acc.merra.mm.svm,acc.merra.mm.rf,acc.merra.mm.brt,acc.merra.mm.mars,acc.merra.mm.maxent,acc.merra.mm.maxlike,acc.merra.mm.glmnet)

comp.accplot(mods=merra.mm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.merra.mm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Medium Model Comparison")
dev.off()

# bm
merra.bm.comp <- list(acc.merra.bm.glm,acc.merra.bm.svm,acc.merra.bm.rf,acc.merra.bm.brt,acc.merra.bm.mars,acc.merra.bm.maxent,acc.merra.bm.maxlike,acc.merra.bm.glmnet)

comp.accplot(mods=merra.bm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.merra.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Big Model Comparison")
dev.off()

### Envirem ####
#### sm ####
bd.envirem.sm <- as.data.frame(envirem.sm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.envirem.sm <- as.data.frame(p.envirem.sm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.envirem.sm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.envirem.sm.glm<-dplyr::select(pd.envirem.sm,long,lat,GLM)
output.mod.envirem.sm.glm<-cbind(pd.envirem.sm.glm,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.glm <- raster::extract(p.envirem.sm[[1]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.glm <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.glm,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.glm <- as.data.frame(occ.pnts.envirem.sm.glm)
occ.pnts.envirem.sm.glm <-drop_na(occ.pnts.envirem.sm.glm)


acc.envirem.sm.glm<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.sm.glm,
                              occ.pnts = occ.pnts.envirem.sm.glm,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.envirem.sm.svm<-dplyr::select(pd.envirem.sm,long,lat,SVM)
output.mod.envirem.sm.svm<-cbind(pd.envirem.sm.svm,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.svm <- raster::extract(p.envirem.sm[[2]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.svm <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.svm,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.svm <- as.data.frame(occ.pnts.envirem.sm.svm)
occ.pnts.envirem.sm.svm <-drop_na(occ.pnts.envirem.sm.svm)


acc.envirem.sm.svm<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.sm.svm,
                              occ.pnts = occ.pnts.envirem.sm.svm,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.envirem.sm.rf<-dplyr::select(pd.envirem.sm,long,lat,RF)
output.mod.envirem.sm.rf<-cbind(pd.envirem.sm.rf,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.rf <- raster::extract(p.envirem.sm[[3]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.rf <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.rf,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.rf <- as.data.frame(occ.pnts.envirem.sm.rf)
occ.pnts.envirem.sm.rf <-drop_na(occ.pnts.envirem.sm.rf)


acc.envirem.sm.rf<-accum.occ(sp.name='Heterophrynus',
                             output.mod = output.mod.envirem.sm.rf,
                             occ.pnts = occ.pnts.envirem.sm.rf,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.envirem.sm.brt<-dplyr::select(pd.envirem.sm,long,lat,BRT)
output.mod.envirem.sm.brt<-cbind(pd.envirem.sm.brt,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.brt <- raster::extract(p.envirem.sm[[4]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.brt <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.brt,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.brt <- as.data.frame(occ.pnts.envirem.sm.brt)
occ.pnts.envirem.sm.brt <-drop_na(occ.pnts.envirem.sm.brt)


acc.envirem.sm.brt<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.sm.brt,
                              occ.pnts = occ.pnts.envirem.sm.brt,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.envirem.sm.mars<-dplyr::select(pd.envirem.sm,long,lat,MARS)
output.mod.envirem.sm.mars<-cbind(pd.envirem.sm.mars,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.mars <- raster::extract(p.envirem.sm[[5]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.mars <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.mars,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.mars <- as.data.frame(occ.pnts.envirem.sm.mars)
occ.pnts.envirem.sm.mars <-drop_na(occ.pnts.envirem.sm.mars)


acc.envirem.sm.mars<-accum.occ(sp.name='Heterophrynus',
                               output.mod = output.mod.envirem.sm.mars,
                               occ.pnts = occ.pnts.envirem.sm.mars,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.envirem.sm.maxent<-dplyr::select(pd.envirem.sm,long,lat,MAXENT)
output.mod.envirem.sm.maxent<-cbind(pd.envirem.sm.maxent,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.maxent <- raster::extract(p.envirem.sm[[6]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.maxent <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.maxent,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.maxent <- as.data.frame(occ.pnts.envirem.sm.maxent)
occ.pnts.envirem.sm.maxent <-drop_na(occ.pnts.envirem.sm.maxent)


acc.envirem.sm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                                 output.mod = output.mod.envirem.sm.maxent,
                                 occ.pnts = occ.pnts.envirem.sm.maxent,
                                 null.mod = "hypergeom",
                                 conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.envirem.sm.maxlike<-dplyr::select(pd.envirem.sm,long,lat,MAXLIKE)
output.mod.envirem.sm.maxlike<-cbind(pd.envirem.sm.maxlike,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.maxlike <- raster::extract(p.envirem.sm[[7]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.maxlike <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.maxlike,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.maxlike <- as.data.frame(occ.pnts.envirem.sm.maxlike)
occ.pnts.envirem.sm.maxlike <-drop_na(occ.pnts.envirem.sm.maxlike)


acc.envirem.sm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                                  output.mod = output.mod.envirem.sm.maxlike,
                                  occ.pnts = occ.pnts.envirem.sm.maxlike,
                                  null.mod = "hypergeom",
                                  conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.envirem.sm.glmnet<-dplyr::select(pd.envirem.sm,long,lat,GLMNET)
output.mod.envirem.sm.glmnet<-cbind(pd.envirem.sm.glmnet,bd.envirem.sm)

occ.p_env.envirem.sm <- raster::extract(envirem.sm,sp,as.data.frame=T)
occ.p_preds.envirem.sm.glmnet <- raster::extract(p.envirem.sm[[8]],sp,as.data.frame=T)
occ.coords.envirem.sm <- sp@coords
occ.pnts.envirem.sm.glmnet <- cbind(occ.coords.envirem.sm,occ.p_preds.envirem.sm.glmnet,occ.p_env.envirem.sm)
occ.pnts.envirem.sm.glmnet <- as.data.frame(occ.pnts.envirem.sm.glmnet)
occ.pnts.envirem.sm.glmnet <-drop_na(occ.pnts.envirem.sm.glmnet)


acc.envirem.sm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                                 output.mod = output.mod.envirem.sm.glmnet,
                                 occ.pnts = occ.pnts.envirem.sm.glmnet,
                                 null.mod = "hypergeom",
                                 conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()
#### mm ####
bd.envirem.mm <- as.data.frame(envirem.mm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.envirem.mm <- as.data.frame(p.envirem.mm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.envirem.mm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.envirem.mm.glm<-dplyr::select(pd.envirem.mm,long,lat,GLM)
output.mod.envirem.mm.glm<-cbind(pd.envirem.mm.glm,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.glm <- raster::extract(p.envirem.mm[[1]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.glm <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.glm,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.glm <- as.data.frame(occ.pnts.envirem.mm.glm)
occ.pnts.envirem.mm.glm <-drop_na(occ.pnts.envirem.mm.glm)


acc.envirem.mm.glm<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.mm.glm,
                              occ.pnts = occ.pnts.envirem.mm.glm,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.envirem.mm.svm<-dplyr::select(pd.envirem.mm,long,lat,SVM)
output.mod.envirem.mm.svm<-cbind(pd.envirem.mm.svm,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.svm <- raster::extract(p.envirem.mm[[2]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.svm <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.svm,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.svm <- as.data.frame(occ.pnts.envirem.mm.svm)
occ.pnts.envirem.mm.svm <-drop_na(occ.pnts.envirem.mm.svm)


acc.envirem.mm.svm<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.mm.svm,
                              occ.pnts = occ.pnts.envirem.mm.svm,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.envirem.mm.rf<-dplyr::select(pd.envirem.mm,long,lat,RF)
output.mod.envirem.mm.rf<-cbind(pd.envirem.mm.rf,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.rf <- raster::extract(p.envirem.mm[[3]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.rf <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.rf,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.rf <- as.data.frame(occ.pnts.envirem.mm.rf)
occ.pnts.envirem.mm.rf <-drop_na(occ.pnts.envirem.mm.rf)


acc.envirem.mm.rf<-accum.occ(sp.name='Heterophrynus',
                             output.mod = output.mod.envirem.mm.rf,
                             occ.pnts = occ.pnts.envirem.mm.rf,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.envirem.mm.brt<-dplyr::select(pd.envirem.mm,long,lat,BRT)
output.mod.envirem.mm.brt<-cbind(pd.envirem.mm.brt,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.brt <- raster::extract(p.envirem.mm[[4]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.brt <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.brt,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.brt <- as.data.frame(occ.pnts.envirem.mm.brt)
occ.pnts.envirem.mm.brt <-drop_na(occ.pnts.envirem.mm.brt)


acc.envirem.mm.brt<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.mm.brt,
                              occ.pnts = occ.pnts.envirem.mm.brt,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.envirem.mm.mars<-dplyr::select(pd.envirem.mm,long,lat,MARS)
output.mod.envirem.mm.mars<-cbind(pd.envirem.mm.mars,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.mars <- raster::extract(p.envirem.mm[[5]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.mars <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.mars,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.mars <- as.data.frame(occ.pnts.envirem.mm.mars)
occ.pnts.envirem.mm.mars <-drop_na(occ.pnts.envirem.mm.mars)


acc.envirem.mm.mars<-accum.occ(sp.name='Heterophrynus',
                               output.mod = output.mod.envirem.mm.mars,
                               occ.pnts = occ.pnts.envirem.mm.mars,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.envirem.mm.maxent<-dplyr::select(pd.envirem.mm,long,lat,MAXENT)
output.mod.envirem.mm.maxent<-cbind(pd.envirem.mm.maxent,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.maxent <- raster::extract(p.envirem.mm[[6]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.maxent <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.maxent,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.maxent <- as.data.frame(occ.pnts.envirem.mm.maxent)
occ.pnts.envirem.mm.maxent <-drop_na(occ.pnts.envirem.mm.maxent)


acc.envirem.mm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                                 output.mod = output.mod.envirem.mm.maxent,
                                 occ.pnts = occ.pnts.envirem.mm.maxent,
                                 null.mod = "hypergeom",
                                 conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.envirem.mm.maxlike<-dplyr::select(pd.envirem.mm,long,lat,MAXLIKE)
output.mod.envirem.mm.maxlike<-cbind(pd.envirem.mm.maxlike,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.maxlike <- raster::extract(p.envirem.mm[[7]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.maxlike <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.maxlike,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.maxlike <- as.data.frame(occ.pnts.envirem.mm.maxlike)
occ.pnts.envirem.mm.maxlike <-drop_na(occ.pnts.envirem.mm.maxlike)


acc.envirem.mm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                                  output.mod = output.mod.envirem.mm.maxlike,
                                  occ.pnts = occ.pnts.envirem.mm.maxlike,
                                  null.mod = "hypergeom",
                                  conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.envirem.mm.glmnet<-dplyr::select(pd.envirem.mm,long,lat,GLMNET)
output.mod.envirem.mm.glmnet<-cbind(pd.envirem.mm.glmnet,bd.envirem.mm)

occ.p_env.envirem.mm <- raster::extract(envirem.mm,sp,as.data.frame=T)
occ.p_preds.envirem.mm.glmnet <- raster::extract(p.envirem.mm[[8]],sp,as.data.frame=T)
occ.coords.envirem.mm <- sp@coords
occ.pnts.envirem.mm.glmnet <- cbind(occ.coords.envirem.mm,occ.p_preds.envirem.mm.glmnet,occ.p_env.envirem.mm)
occ.pnts.envirem.mm.glmnet <- as.data.frame(occ.pnts.envirem.mm.glmnet)
occ.pnts.envirem.mm.glmnet <-drop_na(occ.pnts.envirem.mm.glmnet)


acc.envirem.mm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                                 output.mod = output.mod.envirem.mm.glmnet,
                                 occ.pnts = occ.pnts.envirem.mm.glmnet,
                                 null.mod = "hypergeom",
                                 conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()
#### bm ####
bd.envirem.bm <- as.data.frame(envirem.bm, row.names=NULL, na.rm=F,xy=F,long=F)
pd.envirem.bm <- as.data.frame(p.envirem.bm, row.names=NULL, na.rm=F,xy=T,long=F)
names(pd.envirem.bm)<- c("long","lat","GLM","SVM","RF","BRT","MARS","MAXENT","MAXLIKE","GLMNET")

#GLM
pd.envirem.bm.glm<-dplyr::select(pd.envirem.bm,long,lat,GLM)
output.mod.envirem.bm.glm<-cbind(pd.envirem.bm.glm,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.glm <- raster::extract(p.envirem.bm[[1]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.glm <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.glm,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.glm <- as.data.frame(occ.pnts.envirem.bm.glm)
occ.pnts.envirem.bm.glm <-drop_na(occ.pnts.envirem.bm.glm)


acc.envirem.bm.glm<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.bm.glm,
                              occ.pnts = occ.pnts.envirem.bm.glm,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
# SVM
pd.envirem.bm.svm<-dplyr::select(pd.envirem.bm,long,lat,SVM)
output.mod.envirem.bm.svm<-cbind(pd.envirem.bm.svm,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.svm <- raster::extract(p.envirem.bm[[2]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.svm <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.svm,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.svm <- as.data.frame(occ.pnts.envirem.bm.svm)
occ.pnts.envirem.bm.svm <-drop_na(occ.pnts.envirem.bm.svm)


acc.envirem.bm.svm<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.bm.svm,
                              occ.pnts = occ.pnts.envirem.bm.svm,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#RF
pd.envirem.bm.rf<-dplyr::select(pd.envirem.bm,long,lat,RF)
output.mod.envirem.bm.rf<-cbind(pd.envirem.bm.rf,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.rf <- raster::extract(p.envirem.bm[[3]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.rf <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.rf,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.rf <- as.data.frame(occ.pnts.envirem.bm.rf)
occ.pnts.envirem.bm.rf <-drop_na(occ.pnts.envirem.bm.rf)


acc.envirem.bm.rf<-accum.occ(sp.name='Heterophrynus',
                             output.mod = output.mod.envirem.bm.rf,
                             occ.pnts = occ.pnts.envirem.bm.rf,
                             null.mod = "hypergeom",
                             conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#BRT
pd.envirem.bm.brt<-dplyr::select(pd.envirem.bm,long,lat,BRT)
output.mod.envirem.bm.brt<-cbind(pd.envirem.bm.brt,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.brt <- raster::extract(p.envirem.bm[[4]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.brt <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.brt,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.brt <- as.data.frame(occ.pnts.envirem.bm.brt)
occ.pnts.envirem.bm.brt <-drop_na(occ.pnts.envirem.bm.brt)


acc.envirem.bm.brt<-accum.occ(sp.name='Heterophrynus',
                              output.mod = output.mod.envirem.bm.brt,
                              occ.pnts = occ.pnts.envirem.bm.brt,
                              null.mod = "hypergeom",
                              conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()
#MARS
pd.envirem.bm.mars<-dplyr::select(pd.envirem.bm,long,lat,MARS)
output.mod.envirem.bm.mars<-cbind(pd.envirem.bm.mars,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.mars <- raster::extract(p.envirem.bm[[5]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.mars <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.mars,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.mars <- as.data.frame(occ.pnts.envirem.bm.mars)
occ.pnts.envirem.bm.mars <-drop_na(occ.pnts.envirem.bm.mars)


acc.envirem.bm.mars<-accum.occ(sp.name='Heterophrynus',
                               output.mod = output.mod.envirem.bm.mars,
                               occ.pnts = occ.pnts.envirem.bm.mars,
                               null.mod = "hypergeom",
                               conlev = 0.05, bios = 0)
dev.off () ; dev.off() ; dev.off()

#MAXENT
pd.envirem.bm.maxent<-dplyr::select(pd.envirem.bm,long,lat,MAXENT)
output.mod.envirem.bm.maxent<-cbind(pd.envirem.bm.maxent,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.maxent <- raster::extract(p.envirem.bm[[6]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.maxent <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.maxent,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.maxent <- as.data.frame(occ.pnts.envirem.bm.maxent)
occ.pnts.envirem.bm.maxent <-drop_na(occ.pnts.envirem.bm.maxent)


acc.envirem.bm.maxent<-accum.occ(sp.name='Heterophrynus alces',
                                 output.mod = output.mod.envirem.bm.maxent,
                                 occ.pnts = occ.pnts.envirem.bm.maxent,
                                 null.mod = "hypergeom",
                                 conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#maxlike
pd.envirem.bm.maxlike<-dplyr::select(pd.envirem.bm,long,lat,MAXLIKE)
output.mod.envirem.bm.maxlike<-cbind(pd.envirem.bm.maxlike,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.maxlike <- raster::extract(p.envirem.bm[[7]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.maxlike <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.maxlike,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.maxlike <- as.data.frame(occ.pnts.envirem.bm.maxlike)
occ.pnts.envirem.bm.maxlike <-drop_na(occ.pnts.envirem.bm.maxlike)


acc.envirem.bm.maxlike<-accum.occ(sp.name='Heterophrynus alces',
                                  output.mod = output.mod.envirem.bm.maxlike,
                                  occ.pnts = occ.pnts.envirem.bm.maxlike,
                                  null.mod = "hypergeom",
                                  conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#GLMNet
pd.envirem.bm.glmnet<-dplyr::select(pd.envirem.bm,long,lat,GLMNET)
output.mod.envirem.bm.glmnet<-cbind(pd.envirem.bm.glmnet,bd.envirem.bm)

occ.p_env.envirem.bm <- raster::extract(envirem.bm,sp,as.data.frame=T)
occ.p_preds.envirem.bm.glmnet <- raster::extract(p.envirem.bm[[8]],sp,as.data.frame=T)
occ.coords.envirem.bm <- sp@coords
occ.pnts.envirem.bm.glmnet <- cbind(occ.coords.envirem.bm,occ.p_preds.envirem.bm.glmnet,occ.p_env.envirem.bm)
occ.pnts.envirem.bm.glmnet <- as.data.frame(occ.pnts.envirem.bm.glmnet)
occ.pnts.envirem.bm.glmnet <-drop_na(occ.pnts.envirem.bm.glmnet)


acc.envirem.bm.glmnet<-accum.occ(sp.name='Heterophrynus alces',
                                 output.mod = output.mod.envirem.bm.glmnet,
                                 occ.pnts = occ.pnts.envirem.bm.glmnet,
                                 null.mod = "hypergeom",
                                 conlev = 0.05, bios = 0)
dev.off() ; dev.off() ; dev.off()

#### Generating the AOcCs - comp.accplot ####
# sm 
envirem.sm.comp <- list(acc.envirem.sm.glm,acc.envirem.sm.svm,acc.envirem.sm.rf,acc.envirem.sm.brt,acc.envirem.sm.mars,acc.envirem.sm.maxent,acc.envirem.sm.maxlike,acc.envirem.sm.glmnet)

comp.accplot(mods=envirem.sm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.sm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Small Model Comparison")
dev.off()

# mm
envirem.mm.comp <- list(acc.envirem.mm.glm,acc.envirem.mm.svm,acc.envirem.mm.rf,acc.envirem.mm.brt,acc.envirem.mm.mars,acc.envirem.mm.maxent,acc.envirem.mm.maxlike,acc.envirem.mm.glmnet)

comp.accplot(mods=envirem.mm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.mm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Medium Model Comparison")
dev.off()

# bm
envirem.bm.comp <- list(acc.envirem.bm.glm,acc.envirem.bm.svm,acc.envirem.bm.rf,acc.envirem.bm.brt,acc.envirem.bm.mars,acc.envirem.bm.maxent,acc.envirem.bm.maxlike,acc.envirem.bm.glmnet)

comp.accplot(mods=envirem.bm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('glm','svm','rf','brt','mars','maxent','maxlike','glmnet'), 
             alpha = 0.05)
title("Big Model Comparison")
dev.off()


## Algorithm Comparison among Envs & Ms #####
# To generate the AOcCs in Fig. S2 (by algorithm instead of by clim/M) # simply run the following
#GLM
glm.comp <- list(acc.bio.sm.glm,acc.bio.mm.glm,acc.bio.bm.glm,acc.merra.sm.glm,acc.merra.mm.glm,acc.merra.bm.glm,acc.envirem.sm.glm,acc.envirem.mm.glm,acc.envirem.bm.glm)

comp.accplot(mods=glm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("GLM - Model Comparison")
dev.off()

#SVM
svm.comp <- list(acc.bio.sm.svm,acc.bio.mm.svm,acc.bio.bm.svm,acc.merra.sm.svm,acc.merra.mm.svm,acc.merra.bm.svm,acc.envirem.sm.svm,acc.envirem.mm.svm,acc.envirem.bm.svm)

comp.accplot(mods=svm.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("SVM - Model Comparison")
dev.off()

#RF
rf.comp <- list(acc.bio.sm.rf,acc.bio.mm.rf,acc.bio.bm.rf,acc.merra.sm.rf,acc.merra.mm.rf,acc.merra.bm.rf,acc.envirem.sm.rf,acc.envirem.mm.rf,acc.envirem.bm.rf)

comp.accplot(mods=rf.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("RF - Model Comparison")
dev.off()

#BRT
brt.comp <- list(acc.bio.sm.brt,acc.bio.mm.brt,acc.bio.bm.brt,acc.merra.sm.brt,acc.merra.mm.brt,acc.merra.bm.brt,acc.envirem.sm.brt,acc.envirem.mm.brt,acc.envirem.bm.brt)

comp.accplot(mods=brt.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("BRT - Model Comparison")
dev.off()

#MARS
mars.comp <- list(acc.bio.sm.mars,acc.bio.mm.mars,acc.bio.bm.mars,acc.merra.sm.mars,acc.merra.mm.mars,acc.merra.bm.mars,acc.envirem.sm.mars,acc.envirem.mm.mars,acc.envirem.bm.mars)

comp.accplot(mods=mars.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("MARS - Model Comparison")
dev.off()

#MaxEnt
maxent.comp <- list(acc.bio.sm.maxent,acc.bio.mm.maxent,acc.bio.bm.maxent,acc.merra.sm.maxent,acc.merra.mm.maxent,acc.merra.bm.maxent,acc.envirem.sm.maxent,acc.envirem.mm.maxent,acc.envirem.bm.maxent)

comp.accplot(mods=maxent.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("MaxEnt - Model Comparison")
dev.off()

#MaxLike
maxlike.comp <- list(acc.bio.sm.maxlike,acc.bio.mm.maxlike,acc.bio.bm.maxlike,acc.merra.sm.maxlike,acc.merra.mm.maxlike,acc.merra.bm.maxlike,acc.envirem.sm.maxlike,acc.envirem.mm.maxlike,acc.envirem.bm.maxlike)

comp.accplot(mods=maxlike.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("MaxLike - Model Comparison")
dev.off()

#GLMNet
glmnet.comp <- list(acc.bio.sm.glmnet,acc.bio.mm.glmnet,acc.bio.bm.glmnet,acc.merra.sm.glmnet,acc.merra.mm.glmnet,acc.merra.bm.glmnet,acc.envirem.sm.glmnet,acc.envirem.mm.glmnet,acc.envirem.bm.glmnet)

comp.accplot(mods=glmnet.comp,
             nocc = length(sp),
             ncells = raster::ncell(p.envirem.bm),
             sp.name = 'Heterophrynus alces', 
             mods.names = c('BioClim-SM','BioClim-MM','BioClim-BM','MERRAClim-SM','MERRAClim-MM','MERRAClim-BM','ENVIREM-SM','ENVIREM-MM','ENVIREM-BM'), 
             alpha = 0.05)
title("GLMNet - Model Comparison")
dev.off()
## Saving the models and aocc results ####
setwd(path) ; setwd('Script_output_files')

dir.create("Results")
setwd("Results")
#Save the models
write.sdm(m.bio.sm,"m.bio.sm")
write.sdm(m.bio.mm,"m.bio.mm")
write.sdm(m.bio.bm,"m.bio.bm")
write.sdm(m.merra.sm,"m.merra.sm")
write.sdm(m.merra.mm,"m.merra.mm")
write.sdm(m.merra.bm,"m.merra.bm")
write.sdm(m.envirem.sm,"m.envirem.sm")
write.sdm(m.envirem.mm,"m.envirem.mm")
write.sdm(m.envirem.bm,"m.envirem.bm")
#Save the accum.occ results
bio.sm.comp.output<-plyr::adply(bio.sm.comp,1,unlist,.id = NA); write_csv(bio.sm.comp.output,"bio.sm.results.csv")
bio.mm.comp.output<-plyr::adply(bio.mm.comp,1,unlist,.id = NA); write_csv(bio.mm.comp.output,"bio.mm.results.csv")
bio.bm.comp.output<-plyr::adply(bio.bm.comp,1,unlist,.id = NA); write_csv(bio.bm.comp.output,"bio.bm.results.csv")

merra.sm.comp.output<-plyr::adply(merra.sm.comp,1,unlist,.id = NA); write_csv(merra.sm.comp.output,"merra.sm.results.csv")
merra.mm.comp.output<-plyr::adply(merra.mm.comp,1,unlist,.id = NA); write_csv(merra.mm.comp.output,"merra.mm.results.csv")
merra.bm.comp.output<-plyr::adply(merra.bm.comp,1,unlist,.id = NA); write_csv(merra.bm.comp.output,"merra.bm.results.csv")

envirem.sm.comp.output<-plyr::adply(envirem.sm.comp,1,unlist,.id = NA); write_csv(envirem.sm.comp.output,"envirem.sm.results.csv")
envirem.mm.comp.output<-plyr::adply(envirem.mm.comp,1,unlist,.id = NA); write_csv(envirem.mm.comp.output,"envirem.mm.results.csv")
envirem.bm.comp.output<-plyr::adply(envirem.bm.comp,1,unlist,.id = NA); write_csv(envirem.bm.comp.output,"envirem.bm.results.csv")

# Model Comparison - Schoener's D-stat ####
## Compare Environmental Datasets ####
setwd('..')
setwd('M-preds')
### Load Rasters' Layers individually #####
alces.bio.sm.raster.glm <- raster::raster("predictions.bio.sm.img",band=1) ; alces.bio.sm.raster.glm[alces.bio.sm.raster.glm < 0 ] <- 0
alces.bio.sm.raster.svm <- raster::raster("predictions.bio.sm.img",band=2) ; alces.bio.sm.raster.svm[alces.bio.sm.raster.svm < 0 ] <- 0
alces.bio.sm.raster.rf <- raster::raster("predictions.bio.sm.img",band=3) ; alces.bio.sm.raster.rf[alces.bio.sm.raster.rf < 0 ] <- 0
alces.bio.sm.raster.brt <- raster::raster("predictions.bio.sm.img",band=4) ; alces.bio.sm.raster.brt[alces.bio.sm.raster.brt < 0 ] <- 0
alces.bio.sm.raster.mars <- raster::raster("predictions.bio.sm.img",band=5) ; alces.bio.sm.raster.mars[alces.bio.sm.raster.mars < 0 ] <- 0
alces.bio.sm.raster.maxent <- raster::raster("predictions.bio.sm.img",band=6) ; alces.bio.sm.raster.maxent[alces.bio.sm.raster.maxent < 0 ] <- 0
alces.bio.sm.raster.maxlike <- raster::raster("predictions.bio.sm.img",band=7) ; alces.bio.sm.raster.maxlike[alces.bio.sm.raster.maxlike < 0 ] <- 0
alces.bio.sm.raster.glmnet <- raster::raster("predictions.bio.sm.img",band=8) ;  alces.bio.sm.raster.glmnet <- raster::raster("predictions.bio.sm.img",band=8) ; alces.bio.sm.raster.glmnet <- alces.bio.sm.raster.glmnet + abs(x = alces.bio.sm.raster.glmnet@data@min); alces.bio.sm.raster.glmnet <- alces.bio.sm.raster.glmnet + abs(x = alces.bio.sm.raster.glmnet@data@min)

alces.bio.mm.raster.glm <- raster::raster("predictions.bio.mm.img",band=1) ; alces.bio.mm.raster.glm[alces.bio.mm.raster.glm < 0 ] <- 0
alces.bio.mm.raster.svm <- raster::raster("predictions.bio.mm.img",band=2) ; alces.bio.mm.raster.svm[alces.bio.mm.raster.svm < 0 ] <- 0
alces.bio.mm.raster.rf <- raster::raster("predictions.bio.mm.img",band=3) ; alces.bio.mm.raster.rf[alces.bio.mm.raster.rf < 0 ] <- 0
alces.bio.mm.raster.brt <- raster::raster("predictions.bio.mm.img",band=4) ; alces.bio.mm.raster.brt[alces.bio.mm.raster.brt < 0 ] <- 0
alces.bio.mm.raster.mars <- raster::raster("predictions.bio.mm.img",band=5) ; alces.bio.mm.raster.mars[alces.bio.mm.raster.mars < 0 ] <- 0
alces.bio.mm.raster.maxent <- raster::raster("predictions.bio.mm.img",band=6) ; alces.bio.mm.raster.maxent[alces.bio.mm.raster.maxent < 0 ] <- 0
alces.bio.mm.raster.maxlike <- raster::raster("predictions.bio.mm.img",band=7) ; alces.bio.mm.raster.maxlike[alces.bio.mm.raster.maxlike < 0 ] <- 0
alces.bio.mm.raster.glmnet <- raster::raster("predictions.bio.mm.img",band=8) ; alces.bio.mm.raster.glmnet <- raster::raster("predictions.bio.mm.img",band=8) ; alces.bio.mm.raster.glmnet <- alces.bio.mm.raster.glmnet + abs(x = alces.bio.mm.raster.glmnet@data@min); alces.bio.mm.raster.glmnet <- alces.bio.mm.raster.glmnet + abs(x = alces.bio.mm.raster.glmnet@data@min)

alces.bio.bm.raster.glm <- raster::raster("predictions.bio.bm.img",band=1) ; alces.bio.bm.raster.glm[alces.bio.bm.raster.glm < 0 ] <- 0
alces.bio.bm.raster.svm <- raster::raster("predictions.bio.bm.img",band=2) ; alces.bio.bm.raster.svm[alces.bio.bm.raster.svm < 0 ] <- 0
alces.bio.bm.raster.rf <- raster::raster("predictions.bio.bm.img",band=3) ; alces.bio.bm.raster.rf[alces.bio.bm.raster.rf < 0 ] <- 0
alces.bio.bm.raster.brt <- raster::raster("predictions.bio.bm.img",band=4) ; alces.bio.bm.raster.brt[alces.bio.bm.raster.brt < 0 ] <- 0
alces.bio.bm.raster.mars <- raster::raster("predictions.bio.bm.img",band=5) ; alces.bio.bm.raster.mars[alces.bio.bm.raster.mars < 0 ] <- 0
alces.bio.bm.raster.maxent <- raster::raster("predictions.bio.bm.img",band=6) ; alces.bio.bm.raster.maxent[alces.bio.bm.raster.maxent < 0 ] <- 0
alces.bio.bm.raster.maxlike <- raster::raster("predictions.bio.bm.img",band=7) ; alces.bio.bm.raster.maxlike[alces.bio.bm.raster.maxlike < 0 ] <- 0
alces.bio.bm.raster.glmnet <- raster::raster("predictions.bio.bm.img",band=8) ; alces.bio.bm.raster.glmnet <- raster::raster("predictions.bio.bm.img",band=8) ; alces.bio.bm.raster.glmnet <- alces.bio.bm.raster.glmnet + abs(x = alces.bio.bm.raster.glmnet@data@min); alces.bio.bm.raster.glmnet <- alces.bio.bm.raster.glmnet + abs(x = alces.bio.bm.raster.glmnet@data@min)


alces.merra.sm.raster.glm <- raster::raster("predictions.merra.sm.img",band=1) ; alces.merra.sm.raster.glm[alces.merra.sm.raster.glm < 0 ] <- 0 ; raster::crop(alces.merra.sm.raster.glm,alces.bio.sm.raster.glm); alces.merra.sm.raster.glm@extent<- alces.bio.sm.raster.glm@extent
alces.merra.sm.raster.svm <- raster::raster("predictions.merra.sm.img",band=2) ; alces.merra.sm.raster.svm[alces.merra.sm.raster.svm < 0 ] <- 0 ; raster::crop(alces.merra.sm.raster.svm,alces.bio.sm.raster.svm) ;  alces.merra.sm.raster.svm@extent<- alces.bio.sm.raster.svm@extent
alces.merra.sm.raster.rf <- raster::raster("predictions.merra.sm.img",band=3) ; alces.merra.sm.raster.rf[alces.merra.sm.raster.rf < 0 ] <- 0 ; raster::crop(alces.merra.sm.raster.rf,alces.bio.sm.raster.rf) ; alces.merra.sm.raster.rf@extent<- alces.bio.sm.raster.rf@extent
alces.merra.sm.raster.brt <- raster::raster("predictions.merra.sm.img",band=4) ; alces.merra.sm.raster.brt[alces.merra.sm.raster.brt < 0 ] <- 0 ; raster::crop(alces.merra.sm.raster.brt,alces.bio.sm.raster.brt) ; alces.merra.sm.raster.brt@extent<- alces.bio.sm.raster.brt@extent
alces.merra.sm.raster.mars <- raster::raster("predictions.merra.sm.img",band=5) ; alces.merra.sm.raster.mars[alces.merra.sm.raster.mars < 0 ] <- 0 ; raster::crop(alces.merra.sm.raster.mars,alces.bio.sm.raster.mars) ;  alces.merra.sm.raster.mars@extent<- alces.bio.sm.raster.mars@extent
alces.merra.sm.raster.maxent <- raster::raster("predictions.merra.sm.img",band=6) ; alces.merra.sm.raster.maxent[alces.merra.sm.raster.maxent < 0 ] <- 0 ; raster::crop(alces.merra.sm.raster.maxent,alces.bio.sm.raster.maxent) ; alces.merra.sm.raster.maxent@extent<- alces.bio.sm.raster.maxent@extent
alces.merra.sm.raster.maxlike <- raster::raster("predictions.merra.sm.img",band=7) ; alces.merra.sm.raster.maxlike[alces.merra.sm.raster.maxlike < 0 ] <- 0 ; raster::crop(alces.merra.sm.raster.maxlike,alces.bio.sm.raster.maxlike) ; alces.merra.sm.raster.maxlike@extent<- alces.bio.sm.raster.maxlike@extent
alces.merra.sm.raster.glmnet <- raster::raster("predictions.merra.sm.img",band=8) ; alces.merra.sm.raster.glmnet <- alces.merra.sm.raster.glmnet + abs(x = alces.merra.sm.raster.glmnet@data@min); alces.merra.sm.raster.glmnet <- alces.merra.sm.raster.glmnet + abs(x = alces.merra.sm.raster.glmnet@data@min) ; raster::crop(alces.merra.sm.raster.glmnet,alces.bio.sm.raster.glmnet) ; alces.merra.sm.raster.glmnet@extent<- alces.bio.sm.raster.glmnet@extent

alces.merra.mm.raster.glm <- raster::raster("predictions.merra.mm.img",band=1) ; alces.merra.mm.raster.glm[alces.merra.mm.raster.glm < 0 ] <- 0 ; raster::crop(alces.merra.mm.raster.glm,alces.bio.mm.raster.glm); alces.merra.mm.raster.glm@extent<- alces.bio.mm.raster.glm@extent
alces.merra.mm.raster.svm <- raster::raster("predictions.merra.mm.img",band=2) ; alces.merra.mm.raster.svm[alces.merra.mm.raster.svm < 0 ] <- 0 ; raster::crop(alces.merra.mm.raster.svm,alces.bio.mm.raster.svm) ;  alces.merra.mm.raster.svm@extent<- alces.bio.mm.raster.svm@extent
alces.merra.mm.raster.rf <- raster::raster("predictions.merra.mm.img",band=3) ; alces.merra.mm.raster.rf[alces.merra.mm.raster.rf < 0 ] <- 0 ; raster::crop(alces.merra.mm.raster.rf,alces.bio.mm.raster.rf) ; alces.merra.mm.raster.rf@extent<- alces.bio.mm.raster.rf@extent
alces.merra.mm.raster.brt <- raster::raster("predictions.merra.mm.img",band=4) ; alces.merra.mm.raster.brt[alces.merra.mm.raster.brt < 0 ] <- 0 ; raster::crop(alces.merra.mm.raster.brt,alces.bio.mm.raster.brt) ; alces.merra.mm.raster.brt@extent<- alces.bio.mm.raster.brt@extent
alces.merra.mm.raster.mars <- raster::raster("predictions.merra.mm.img",band=5) ; alces.merra.mm.raster.mars[alces.merra.mm.raster.mars < 0 ] <- 0 ; raster::crop(alces.merra.mm.raster.mars,alces.bio.mm.raster.mars) ;  alces.merra.mm.raster.mars@extent<- alces.bio.mm.raster.mars@extent
alces.merra.mm.raster.maxent <- raster::raster("predictions.merra.mm.img",band=6) ; alces.merra.mm.raster.maxent[alces.merra.mm.raster.maxent < 0 ] <- 0 ; raster::crop(alces.merra.mm.raster.maxent,alces.bio.mm.raster.maxent) ; alces.merra.mm.raster.maxent@extent<- alces.bio.mm.raster.maxent@extent
alces.merra.mm.raster.maxlike <- raster::raster("predictions.merra.mm.img",band=7) ; alces.merra.mm.raster.maxlike[alces.merra.mm.raster.maxlike < 0 ] <- 0 ; raster::crop(alces.merra.mm.raster.maxlike,alces.bio.mm.raster.maxlike) ; alces.merra.mm.raster.maxlike@extent<- alces.bio.mm.raster.maxlike@extent
alces.merra.mm.raster.glmnet <- raster::raster("predictions.merra.mm.img",band=8) ; alces.merra.mm.raster.glmnet <- alces.merra.mm.raster.glmnet + abs(x = alces.merra.mm.raster.glmnet@data@min); alces.merra.mm.raster.glmnet <- alces.merra.mm.raster.glmnet + abs(x = alces.merra.mm.raster.glmnet@data@min) ; raster::crop(alces.merra.mm.raster.glmnet,alces.bio.mm.raster.glmnet) ; alces.merra.mm.raster.glmnet@extent<- alces.bio.mm.raster.glmnet@extent

alces.merra.bm.raster.glm <- raster::raster("predictions.merra.bm.img",band=1) ; alces.merra.bm.raster.glm[alces.merra.bm.raster.glm < 0 ] <- 0 ; raster::crop(alces.merra.bm.raster.glm,alces.bio.bm.raster.glm); alces.merra.bm.raster.glm@extent<- alces.bio.bm.raster.glm@extent
alces.merra.bm.raster.svm <- raster::raster("predictions.merra.bm.img",band=2) ; alces.merra.bm.raster.svm[alces.merra.bm.raster.svm < 0 ] <- 0 ; raster::crop(alces.merra.bm.raster.svm,alces.bio.bm.raster.svm) ;  alces.merra.bm.raster.svm@extent<- alces.bio.bm.raster.svm@extent
alces.merra.bm.raster.rf <- raster::raster("predictions.merra.bm.img",band=3) ; alces.merra.bm.raster.rf[alces.merra.bm.raster.rf < 0 ] <- 0 ; raster::crop(alces.merra.bm.raster.rf,alces.bio.bm.raster.rf) ; alces.merra.bm.raster.rf@extent<- alces.bio.bm.raster.rf@extent
alces.merra.bm.raster.brt <- raster::raster("predictions.merra.bm.img",band=4) ; alces.merra.bm.raster.brt[alces.merra.bm.raster.brt < 0 ] <- 0 ; raster::crop(alces.merra.bm.raster.brt,alces.bio.bm.raster.brt) ; alces.merra.bm.raster.brt@extent<- alces.bio.bm.raster.brt@extent
alces.merra.bm.raster.mars <- raster::raster("predictions.merra.bm.img",band=5) ; alces.merra.bm.raster.mars[alces.merra.bm.raster.mars < 0 ] <- 0 ; raster::crop(alces.merra.bm.raster.mars,alces.bio.bm.raster.mars) ;  alces.merra.bm.raster.mars@extent<- alces.bio.bm.raster.mars@extent
alces.merra.bm.raster.maxent <- raster::raster("predictions.merra.bm.img",band=6) ; alces.merra.bm.raster.maxent[alces.merra.bm.raster.maxent < 0 ] <- 0 ; raster::crop(alces.merra.bm.raster.maxent,alces.bio.bm.raster.maxent) ; alces.merra.bm.raster.maxent@extent<- alces.bio.bm.raster.maxent@extent
alces.merra.bm.raster.maxlike <- raster::raster("predictions.merra.bm.img",band=7) ; alces.merra.bm.raster.maxlike[alces.merra.bm.raster.maxlike < 0 ] <- 0 ; raster::crop(alces.merra.bm.raster.maxlike,alces.bio.bm.raster.maxlike) ; alces.merra.bm.raster.maxlike@extent<- alces.bio.bm.raster.maxlike@extent
alces.merra.bm.raster.glmnet <- raster::raster("predictions.merra.bm.img",band=8) ; alces.merra.bm.raster.glmnet <- alces.merra.bm.raster.glmnet + abs(x = alces.merra.bm.raster.glmnet@data@min); alces.merra.bm.raster.glmnet <- alces.merra.bm.raster.glmnet + abs(x = alces.merra.bm.raster.glmnet@data@min) ; raster::crop(alces.merra.bm.raster.glmnet,alces.bio.bm.raster.glmnet) ; alces.merra.bm.raster.glmnet@extent<- alces.bio.bm.raster.glmnet@extent

alces.envirem.sm.raster.glm <- raster::raster("predictions.envirem.sm.img",band=1) ; alces.envirem.sm.raster.glm[alces.envirem.sm.raster.glm < 0 ] <- 0
alces.envirem.sm.raster.svm <- raster::raster("predictions.envirem.sm.img",band=2) ; alces.envirem.sm.raster.svm[alces.envirem.sm.raster.svm < 0 ] <- 0
alces.envirem.sm.raster.rf <- raster::raster("predictions.envirem.sm.img",band=3) ; alces.envirem.sm.raster.rf[alces.envirem.sm.raster.rf < 0 ] <- 0
alces.envirem.sm.raster.brt <- raster::raster("predictions.envirem.sm.img",band=4) ; alces.envirem.sm.raster.brt[alces.envirem.sm.raster.brt < 0 ] <- 0
alces.envirem.sm.raster.mars <- raster::raster("predictions.envirem.sm.img",band=5) ; alces.envirem.sm.raster.mars[alces.envirem.sm.raster.mars < 0 ] <- 0
alces.envirem.sm.raster.maxent <- raster::raster("predictions.envirem.sm.img",band=6) ; alces.envirem.sm.raster.maxent[alces.envirem.sm.raster.maxent < 0 ] <- 0
alces.envirem.sm.raster.maxlike <- raster::raster("predictions.envirem.sm.img",band=7) ; alces.envirem.sm.raster.maxlike[alces.envirem.sm.raster.maxlike < 0 ] <- 0
alces.envirem.sm.raster.glmnet <- raster::raster("predictions.envirem.sm.img",band=8) ; alces.envirem.sm.raster.glmnet <- alces.envirem.sm.raster.glmnet + abs(x = alces.envirem.sm.raster.glmnet@data@min); alces.envirem.sm.raster.glmnet <- alces.envirem.sm.raster.glmnet + abs(x = alces.envirem.sm.raster.glmnet@data@min)

alces.envirem.mm.raster.glm <- raster::raster("predictions.envirem.mm.img",band=1) ; alces.envirem.mm.raster.glm[alces.envirem.mm.raster.glm < 0 ] <- 0
alces.envirem.mm.raster.svm <- raster::raster("predictions.envirem.mm.img",band=2) ; alces.envirem.mm.raster.svm[alces.envirem.mm.raster.svm < 0 ] <- 0
alces.envirem.mm.raster.rf <- raster::raster("predictions.envirem.mm.img",band=3) ; alces.envirem.mm.raster.rf[alces.envirem.mm.raster.rf < 0 ] <- 0
alces.envirem.mm.raster.brt <- raster::raster("predictions.envirem.mm.img",band=4) ; alces.envirem.mm.raster.brt[alces.envirem.mm.raster.brt < 0 ] <- 0
alces.envirem.mm.raster.mars <- raster::raster("predictions.envirem.mm.img",band=5) ; alces.envirem.mm.raster.mars[alces.envirem.mm.raster.mars < 0 ] <- 0
alces.envirem.mm.raster.maxent <- raster::raster("predictions.envirem.mm.img",band=6) ; alces.envirem.mm.raster.maxent[alces.envirem.mm.raster.maxent < 0 ] <- 0
alces.envirem.mm.raster.maxlike <- raster::raster("predictions.envirem.mm.img",band=7) ; alces.envirem.mm.raster.maxlike[alces.envirem.mm.raster.maxlike < 0 ] <- 0
alces.envirem.mm.raster.glmnet <- raster::raster("predictions.envirem.mm.img",band=8) ; alces.envirem.mm.raster.glmnet <- alces.envirem.mm.raster.glmnet + abs(x = alces.envirem.mm.raster.glmnet@data@min); alces.envirem.mm.raster.glmnet <- alces.envirem.mm.raster.glmnet + abs(x = alces.envirem.mm.raster.glmnet@data@min)

alces.envirem.bm.raster.glm <- raster::raster("predictions.envirem.bm.img",band=1) ; alces.envirem.bm.raster.glm[alces.envirem.bm.raster.glm < 0 ] <- 0
alces.envirem.bm.raster.svm <- raster::raster("predictions.envirem.bm.img",band=2) ; alces.envirem.bm.raster.svm[alces.envirem.bm.raster.svm < 0 ] <- 0
alces.envirem.bm.raster.rf <- raster::raster("predictions.envirem.bm.img",band=3) ; alces.envirem.bm.raster.rf[alces.envirem.bm.raster.rf < 0 ] <- 0
alces.envirem.bm.raster.brt <- raster::raster("predictions.envirem.bm.img",band=4) ; alces.envirem.bm.raster.brt[alces.envirem.bm.raster.brt < 0 ] <- 0
alces.envirem.bm.raster.mars <- raster::raster("predictions.envirem.bm.img",band=5) ; alces.envirem.bm.raster.mars[alces.envirem.bm.raster.mars < 0 ] <- 0
alces.envirem.bm.raster.maxent <- raster::raster("predictions.envirem.bm.img",band=6) ; alces.envirem.bm.raster.maxent[alces.envirem.bm.raster.maxent < 0 ] <- 0
alces.envirem.bm.raster.maxlike <- raster::raster("predictions.envirem.bm.img",band=7) ; alces.envirem.bm.raster.maxlike[alces.envirem.bm.raster.maxlike < 0 ] <- 0
alces.envirem.bm.raster.glmnet <- raster::raster("predictions.envirem.bm.img",band=8) ; alces.envirem.bm.raster.glmnet <- alces.envirem.bm.raster.glmnet + abs(x = alces.envirem.bm.raster.glmnet@data@min); alces.envirem.bm.raster.glmnet <- alces.envirem.bm.raster.glmnet + abs(x = alces.envirem.bm.raster.glmnet@data@min)

#### Analyze - nicheOverlap #####
#GLM
alces.bioxmerra.sm.glm <- nicheOverlap(alces.bio.sm.raster.glm,alces.merra.sm.raster.glm,stat = 'D')
alces.bioxenvirem.sm.glm <- nicheOverlap(alces.bio.sm.raster.glm,alces.envirem.sm.raster.glm,stat = 'D')
alces.enviremxmerra.sm.glm <- nicheOverlap(alces.envirem.sm.raster.glm,alces.merra.sm.raster.glm,stat = 'D')

alces.bioxmerra.mm.glm <- nicheOverlap(alces.bio.mm.raster.glm,alces.merra.mm.raster.glm,stat = 'D')
alces.bioxenvirem.mm.glm <- nicheOverlap(alces.bio.mm.raster.glm,alces.envirem.mm.raster.glm,stat = 'D')
alces.enviremxmerra.mm.glm <- nicheOverlap(alces.envirem.mm.raster.glm,alces.merra.mm.raster.glm,stat = 'D')

alces.bioxmerra.bm.glm <- nicheOverlap(alces.bio.bm.raster.glm,alces.merra.bm.raster.glm,stat = 'D')
alces.bioxenvirem.bm.glm <- nicheOverlap(alces.bio.bm.raster.glm,alces.envirem.bm.raster.glm,stat = 'D')
alces.enviremxmerra.bm.glm <- nicheOverlap(alces.envirem.bm.raster.glm,alces.merra.bm.raster.glm,stat = 'D')

#SVM
alces.bioxmerra.sm.svm <- nicheOverlap(alces.bio.sm.raster.svm,alces.merra.sm.raster.svm,stat = 'D')
alces.bioxenvirem.sm.svm <- nicheOverlap(alces.bio.sm.raster.svm,alces.envirem.sm.raster.svm,stat = 'D')
alces.enviremxmerra.sm.svm <- nicheOverlap(alces.envirem.sm.raster.svm,alces.merra.sm.raster.svm,stat = 'D')

alces.bioxmerra.mm.svm <- nicheOverlap(alces.bio.mm.raster.svm,alces.merra.mm.raster.svm,stat = 'D')
alces.bioxenvirem.mm.svm <- nicheOverlap(alces.bio.mm.raster.svm,alces.envirem.mm.raster.svm,stat = 'D')
alces.enviremxmerra.mm.svm <- nicheOverlap(alces.envirem.mm.raster.svm,alces.merra.mm.raster.svm,stat = 'D')

alces.bioxmerra.bm.svm <- nicheOverlap(alces.bio.bm.raster.svm,alces.merra.bm.raster.svm,stat = 'D')
alces.bioxenvirem.bm.svm <- nicheOverlap(alces.bio.bm.raster.svm,alces.envirem.bm.raster.svm,stat = 'D')
alces.enviremxmerra.bm.svm <- nicheOverlap(alces.envirem.bm.raster.svm,alces.merra.bm.raster.svm,stat = 'D')

#RF
alces.bioxmerra.sm.rf <- nicheOverlap(alces.bio.sm.raster.rf,alces.merra.sm.raster.rf,stat = 'D')
alces.bioxenvirem.sm.rf <- nicheOverlap(alces.bio.sm.raster.rf,alces.envirem.sm.raster.rf,stat = 'D')
alces.enviremxmerra.sm.rf <- nicheOverlap(alces.envirem.sm.raster.rf,alces.merra.sm.raster.rf,stat = 'D')

alces.bioxmerra.mm.rf <- nicheOverlap(alces.bio.mm.raster.rf,alces.merra.mm.raster.rf,stat = 'D')
alces.bioxenvirem.mm.rf <- nicheOverlap(alces.bio.mm.raster.rf,alces.envirem.mm.raster.rf,stat = 'D')
alces.enviremxmerra.mm.rf <- nicheOverlap(alces.envirem.mm.raster.rf,alces.merra.mm.raster.rf,stat = 'D')

alces.bioxmerra.bm.rf <- nicheOverlap(alces.bio.bm.raster.rf,alces.merra.bm.raster.rf,stat = 'D')
alces.bioxenvirem.bm.rf <- nicheOverlap(alces.bio.bm.raster.rf,alces.envirem.bm.raster.rf,stat = 'D')
alces.enviremxmerra.bm.rf <- nicheOverlap(alces.envirem.bm.raster.rf,alces.merra.bm.raster.rf,stat = 'D')

#BRT
alces.bioxmerra.sm.brt <- nicheOverlap(alces.bio.sm.raster.brt,alces.merra.sm.raster.brt,stat = 'D')
alces.bioxenvirem.sm.brt <- nicheOverlap(alces.bio.sm.raster.brt,alces.envirem.sm.raster.brt,stat = 'D')
alces.enviremxmerra.sm.brt <- nicheOverlap(alces.envirem.sm.raster.brt,alces.merra.sm.raster.brt,stat = 'D')

alces.bioxmerra.mm.brt <- nicheOverlap(alces.bio.mm.raster.brt,alces.merra.mm.raster.brt,stat = 'D')
alces.bioxenvirem.mm.brt <- nicheOverlap(alces.bio.mm.raster.brt,alces.envirem.mm.raster.brt,stat = 'D')
alces.enviremxmerra.mm.brt <- nicheOverlap(alces.envirem.mm.raster.brt,alces.merra.mm.raster.brt,stat = 'D')

alces.bioxmerra.bm.brt <- nicheOverlap(alces.bio.bm.raster.brt,alces.merra.bm.raster.brt,stat = 'D')
alces.bioxenvirem.bm.brt <- nicheOverlap(alces.bio.bm.raster.brt,alces.envirem.bm.raster.brt,stat = 'D')
alces.enviremxmerra.bm.brt <- nicheOverlap(alces.envirem.bm.raster.brt,alces.merra.bm.raster.brt,stat = 'D')

#MARS
alces.bioxmerra.sm.mars <- nicheOverlap(alces.bio.sm.raster.mars,alces.merra.sm.raster.mars,stat = 'D')
alces.bioxenvirem.sm.mars <- nicheOverlap(alces.bio.sm.raster.mars,alces.envirem.sm.raster.mars,stat = 'D')
alces.enviremxmerra.sm.mars <- nicheOverlap(alces.envirem.sm.raster.mars,alces.merra.sm.raster.mars,stat = 'D')

alces.bioxmerra.mm.mars <- nicheOverlap(alces.bio.mm.raster.mars,alces.merra.mm.raster.mars,stat = 'D')
alces.bioxenvirem.mm.mars <- nicheOverlap(alces.bio.mm.raster.mars,alces.envirem.mm.raster.mars,stat = 'D')
alces.enviremxmerra.mm.mars <- nicheOverlap(alces.envirem.mm.raster.mars,alces.merra.mm.raster.mars,stat = 'D')

alces.bioxmerra.bm.mars <- nicheOverlap(alces.bio.bm.raster.mars,alces.merra.bm.raster.mars,stat = 'D')
alces.bioxenvirem.bm.mars <- nicheOverlap(alces.bio.bm.raster.mars,alces.envirem.bm.raster.mars,stat = 'D')
alces.enviremxmerra.bm.mars <- nicheOverlap(alces.envirem.bm.raster.mars,alces.merra.bm.raster.mars,stat = 'D')

#MAXENT
alces.bioxmerra.sm.maxent <- nicheOverlap(alces.bio.sm.raster.maxent,alces.merra.sm.raster.maxent,stat = 'D')
alces.bioxenvirem.sm.maxent <- nicheOverlap(alces.bio.sm.raster.maxent,alces.envirem.sm.raster.maxent,stat = 'D')
alces.enviremxmerra.sm.maxent <- nicheOverlap(alces.envirem.sm.raster.maxent,alces.merra.sm.raster.maxent,stat = 'D')

alces.bioxmerra.mm.maxent <- nicheOverlap(alces.bio.mm.raster.maxent,alces.merra.mm.raster.maxent,stat = 'D')
alces.bioxenvirem.mm.maxent <- nicheOverlap(alces.bio.mm.raster.maxent,alces.envirem.mm.raster.maxent,stat = 'D')
alces.enviremxmerra.mm.maxent <- nicheOverlap(alces.envirem.mm.raster.maxent,alces.merra.mm.raster.maxent,stat = 'D')

alces.bioxmerra.bm.maxent <- nicheOverlap(alces.bio.bm.raster.maxent,alces.merra.bm.raster.maxent,stat = 'D')
alces.bioxenvirem.bm.maxent <- nicheOverlap(alces.bio.bm.raster.maxent,alces.envirem.bm.raster.maxent,stat = 'D')
alces.enviremxmerra.bm.maxent <- nicheOverlap(alces.envirem.bm.raster.maxent,alces.merra.bm.raster.maxent,stat = 'D')

#MaxLike
alces.bioxmerra.sm.maxlike <- nicheOverlap(alces.bio.sm.raster.maxlike,alces.merra.sm.raster.maxlike,stat = 'D')
alces.bioxenvirem.sm.maxlike <- nicheOverlap(alces.bio.sm.raster.maxlike,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.enviremxmerra.sm.maxlike <- nicheOverlap(alces.envirem.sm.raster.maxlike,alces.merra.sm.raster.maxlike,stat = 'D')

alces.bioxmerra.mm.maxlike <- nicheOverlap(alces.bio.mm.raster.maxlike,alces.merra.mm.raster.maxlike,stat = 'D')
alces.bioxenvirem.mm.maxlike <- nicheOverlap(alces.bio.mm.raster.maxlike,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.enviremxmerra.mm.maxlike <- nicheOverlap(alces.envirem.mm.raster.maxlike,alces.merra.mm.raster.maxlike,stat = 'D')

alces.bioxmerra.bm.maxlike <- nicheOverlap(alces.bio.bm.raster.maxlike,alces.merra.bm.raster.maxlike,stat = 'D')
alces.bioxenvirem.bm.maxlike <- nicheOverlap(alces.bio.bm.raster.maxlike,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.enviremxmerra.bm.maxlike <- nicheOverlap(alces.envirem.bm.raster.maxlike,alces.merra.bm.raster.maxlike,stat = 'D')

#GLMnet
alces.bioxmerra.sm.glmnet <- nicheOverlap(alces.bio.sm.raster.glmnet,alces.merra.sm.raster.glmnet,stat = 'D')
alces.bioxenvirem.sm.glmnet <- nicheOverlap(alces.bio.sm.raster.glmnet,alces.envirem.sm.raster.glmnet,stat = 'D')
alces.enviremxmerra.sm.glmnet <- nicheOverlap(alces.envirem.sm.raster.glmnet,alces.merra.sm.raster.glmnet,stat = 'D')

alces.bioxmerra.mm.glmnet <- nicheOverlap(alces.bio.mm.raster.glmnet,alces.merra.mm.raster.glmnet,stat = 'D')
alces.bioxenvirem.mm.glmnet <- nicheOverlap(alces.bio.mm.raster.glmnet,alces.envirem.mm.raster.glmnet,stat = 'D')
alces.enviremxmerra.mm.glmnet <- nicheOverlap(alces.envirem.mm.raster.glmnet,alces.merra.mm.raster.glmnet,stat = 'D')

alces.bioxmerra.bm.glmnet <- nicheOverlap(alces.bio.bm.raster.glmnet,alces.merra.bm.raster.glmnet,stat = 'D')
alces.bioxenvirem.bm.glmnet <- nicheOverlap(alces.bio.bm.raster.glmnet,alces.envirem.bm.raster.glmnet,stat = 'D')
alces.enviremxmerra.bm.glmnet <- nicheOverlap(alces.envirem.bm.raster.glmnet,alces.merra.bm.raster.glmnet,stat = 'D')

#### Organize the results #####
# here, we'll create three vectors, one for SM, one for MM and one for BM. 
# vectors will be standardized, so we can build a dataframe from it
# BE = Bioclim x Envirem ; EM = Envirem x Merraclim ; BM = Bioclim x Merraclim

sm.alces <- c(alces.bioxenvirem.sm.glm,alces.bioxmerra.sm.glm,alces.enviremxmerra.sm.glm,alces.bioxenvirem.sm.svm,alces.bioxmerra.sm.svm,alces.enviremxmerra.sm.svm,alces.bioxenvirem.sm.rf,alces.bioxmerra.sm.rf,alces.enviremxmerra.sm.rf,alces.bioxenvirem.sm.brt,alces.bioxmerra.sm.brt,alces.enviremxmerra.sm.brt,alces.bioxenvirem.sm.mars,alces.bioxmerra.sm.mars,alces.enviremxmerra.sm.mars,alces.bioxenvirem.sm.maxent,alces.bioxmerra.sm.maxent,alces.enviremxmerra.sm.maxent,alces.bioxenvirem.sm.maxlike,alces.bioxmerra.sm.maxlike,alces.enviremxmerra.sm.maxlike,alces.bioxenvirem.sm.glmnet,alces.bioxmerra.sm.glmnet,alces.enviremxmerra.sm.glmnet)

mm.alces <- c(alces.bioxenvirem.mm.glm,alces.bioxmerra.mm.glm,alces.enviremxmerra.mm.glm,alces.bioxenvirem.mm.svm,alces.bioxmerra.mm.svm,alces.enviremxmerra.mm.svm,alces.bioxenvirem.mm.rf,alces.bioxmerra.mm.rf,alces.enviremxmerra.mm.rf,alces.bioxenvirem.mm.brt,alces.bioxmerra.mm.brt,alces.enviremxmerra.mm.brt,alces.bioxenvirem.mm.mars,alces.bioxmerra.mm.mars,alces.enviremxmerra.mm.mars,alces.bioxenvirem.mm.maxent,alces.bioxmerra.mm.maxent,alces.enviremxmerra.mm.maxent,alces.bioxenvirem.mm.maxlike,alces.bioxmerra.mm.maxlike,alces.enviremxmerra.mm.maxlike,alces.bioxenvirem.mm.glmnet,alces.bioxmerra.mm.glmnet,alces.enviremxmerra.mm.glmnet)

bm.alces <- c(alces.bioxenvirem.bm.glm,alces.bioxmerra.bm.glm,alces.enviremxmerra.bm.glm,alces.bioxenvirem.bm.svm,alces.bioxmerra.bm.svm,alces.enviremxmerra.bm.svm,alces.bioxenvirem.bm.rf,alces.bioxmerra.bm.rf,alces.enviremxmerra.bm.rf,alces.bioxenvirem.bm.brt,alces.bioxmerra.bm.brt,alces.enviremxmerra.bm.brt,alces.bioxenvirem.bm.mars,alces.bioxmerra.bm.mars,alces.enviremxmerra.bm.mars,alces.bioxenvirem.bm.maxent,alces.bioxmerra.bm.maxent,alces.enviremxmerra.bm.maxent,alces.bioxenvirem.bm.maxlike,alces.bioxmerra.bm.maxlike,alces.enviremxmerra.bm.maxlike,alces.bioxenvirem.bm.glmnet,alces.bioxmerra.bm.glmnet,alces.enviremxmerra.bm.glmnet)

setwd(path)

### Put it all together #####
# bind s, m and b Ms

colnamesdf <- c('BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim','BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim','BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim','BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim','BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim','BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim','BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim','BioClim x Envirem','BioClim x MERRAclim','Envirem x MERRAclim')
methodrow <- c('GLM','GLM','GLM','SVM','SVM','SVM','RF','RF','RF','BRT','BRT','BRT','MARS','MARS','MARS','MAXENT','MAXENT','MAXENT','MAXLIKE','MAXLIKE','MAXLIKE','GLMNET','GLMNET','GLMNET')

sm.alces1 <- rbind(colnamesdf,methodrow,sm.alces) ; sm.alces1 <- t(sm.alces1) ; sm.alces1 <- as.data.frame(sm.alces1) ; sm.alces1$species <- 'H. alces (n=21)'; colnames(sm.alces1) <- c('Niche Overlap','algorithm','values','species')

# In my paper, this is where all species came together. You can follow the same process that created sm.alces1 for the other species, and bind them in a single sm object:

# sm <- dplyr::bind_rows(sm.cervinus1,sm.armiger1,sm.alces1,sm.cheiracanthus1,sm.elaphus1,sm.boterorum1,sm.vesanicus1,sm.batesii1,sm.longicornis1)

mm.alces1 <- rbind(colnamesdf,methodrow,mm.alces) ; mm.alces1 <- t(mm.alces1) ; mm.alces1 <- as.data.frame(mm.alces1) ; mm.alces1$species <- 'H. alces (n=21)'; colnames(mm.alces1) <- c('Niche Overlap','algorithm','values','species')

#same:
#mm <- dplyr::bind_rows(mm.cervinus1,mm.armiger1,mm.alces1,mm.cheiracanthus1,mm.elaphus1,mm.boterorum1,mm.vesanicus1,mm.batesii1,mm.longicornis1)

bm.alces1 <- rbind(colnamesdf,methodrow,bm.alces) ; bm.alces1 <- t(bm.alces1) ; bm.alces1 <- as.data.frame(bm.alces1) ; bm.alces1$species <- 'H. alces (n=21)'; colnames(bm.alces1) <- c('Niche Overlap','algorithm','values','species')

#same:
#bm <- dplyr::bind_rows(bm.cervinus1,bm.armiger1,bm.alces1,bm.cheiracanthus1,bm.elaphus1,bm.boterorum1,bm.vesanicus1,bm.batesii1,bm.longicornis1)

##### BoxPlot Similarities #####
# Note: the code in the next few lines won't work, because we did not create an M object (see line 2253), which was created by me including nine species, so the next lines are just an example of how to recreate Figure S8

# Figure S8 was generated with the following lines, but it can't be done here for this script covers a single species:
# 
# sm <- as_tibble(sm) ; sm$`Niche Overlap` <- as.factor(sm$`Niche Overlap`) ; sm$algorithm <- as.factor(sm$algorithm) ; sm$species <- factor(sm$species,levels = unique(sm$species)) ; sm$values <- as.numeric(sm$values) ; sm$msize <- as.factor('Small M')
# mm <- as_tibble(mm) ; mm$`Niche Overlap` <- as.factor(mm$`Niche Overlap`) ; mm$algorithm <- as.factor(mm$algorithm) ; mm$species <- factor(mm$species,levels = unique(mm$species)) ; mm$values <- as.numeric(mm$values) ; mm$msize <- as.factor('Medium M')
# bm <- as_tibble(bm) ; bm$`Niche Overlap` <- as.factor(bm$`Niche Overlap`) ; bm$algorithm <- as.factor(bm$algorithm)  ; bm$species <- factor(bm$species,levels = unique(bm$species)) ;  bm$values <- as.numeric(bm$values);  bm$msize <- as.factor('Large M')
# 
# M <- dplyr::bind_rows(sm,mm,bm)
# names(M) <- c('nicheOverlap','algorithm','values','species','msize')
# M %>%
#   ggplot( aes(x=nicheOverlap, y=values, fill=nicheOverlap)) +
#   geom_violin(width=0.4,alpha=0.2) +
#   geom_boxplot(width=0.1, color="black", alpha=0.8)+
#   geom_jitter(color="black", size=0.5, alpha=0.15, width=0.15) +
#   scale_fill_viridis(discrete = TRUE,alpha = 0.7) +
#   scale_x_discrete(name='Environmental Datasets') +
#   scale_y_continuous(name='Schoener´s D Statistic')+
#   theme_bw()+
#   theme(legend.position = "none",
#         legend.text = element_text(size=15),
#         legend.title = element_text(size= 16),
#         plot.title = element_text(size=20),
#         axis.title.y = element_text(vjust=2,hjust=0.5,size=20)
#   )#+
# #facet_grid(msize~algorithm,scales = 'free_y') #comment and Uncomment this line to generate different plots. FigS8 = uncommented

setwd(path) ; setwd('Script-Output-Files') ; dir.create('D-stat-ENVs') ; setwd('D-stat-ENVs')
# ggsave(filename='d-stat-boxplot_byalgM.png',device = 'png',dpi='print',limitsize = F,width = 10, height = 10,bg = 'white')
write.csv(M,'d-stat_env-results.csv')
## Compare Algorithms ####
### BioClim ####
# Bio - SM
alces.rfxbrt.sm.bio <- nicheOverlap(alces.bio.sm.raster.rf,alces.bio.sm.raster.brt,stat = 'D')
alces.rfxglm.sm.bio <- nicheOverlap(alces.bio.sm.raster.rf,alces.envirem.sm.raster.glm,stat = 'D')
alces.rfxglmnet.sm.bio <- nicheOverlap(alces.bio.sm.raster.rf,alces.bio.sm.raster.glmnet,stat = 'D')
alces.rfxmars.sm.bio <- nicheOverlap(alces.bio.sm.raster.rf,alces.envirem.sm.raster.mars,stat = 'D')
alces.rfxmaxent.sm.bio <- nicheOverlap(alces.bio.sm.raster.rf,alces.bio.sm.raster.maxent,stat = 'D')
alces.rfxmaxlike.sm.bio <- nicheOverlap(alces.bio.sm.raster.rf,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.rfxsvm.sm.bio <- nicheOverlap(alces.bio.sm.raster.rf,alces.bio.sm.raster.svm,stat = 'D')

alces.brtxglm.sm.bio <- nicheOverlap(alces.bio.sm.raster.brt,alces.envirem.sm.raster.glm,stat = 'D')
alces.brtxglmnet.sm.bio <- nicheOverlap(alces.bio.sm.raster.brt,alces.bio.sm.raster.glmnet,stat = 'D')
alces.brtxmars.sm.bio <- nicheOverlap(alces.bio.sm.raster.brt,alces.envirem.sm.raster.mars,stat = 'D')
alces.brtxmaxent.sm.bio <- nicheOverlap(alces.bio.sm.raster.brt,alces.bio.sm.raster.maxent,stat = 'D')
alces.brtxmaxlike.sm.bio <- nicheOverlap(alces.bio.sm.raster.brt,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.brtxsvm.sm.bio <- nicheOverlap(alces.bio.sm.raster.brt,alces.bio.sm.raster.svm,stat = 'D')

alces.glmxglmnet.sm.bio <- nicheOverlap(alces.bio.sm.raster.glm,alces.bio.sm.raster.glmnet,stat = 'D')
alces.glmxmars.sm.bio <- nicheOverlap(alces.bio.sm.raster.glm,alces.envirem.sm.raster.mars,stat = 'D')
alces.glmxmaxent.sm.bio <- nicheOverlap(alces.bio.sm.raster.glm,alces.bio.sm.raster.maxent,stat = 'D')
alces.glmxmaxlike.sm.bio <- nicheOverlap(alces.bio.sm.raster.glm,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.glmxsvm.sm.bio <- nicheOverlap(alces.bio.sm.raster.glm,alces.bio.sm.raster.svm,stat = 'D')

alces.glmnetxmars.sm.bio <- nicheOverlap(alces.bio.sm.raster.glmnet,alces.envirem.sm.raster.mars,stat = 'D')
alces.glmnetxmaxent.sm.bio <- nicheOverlap(alces.bio.sm.raster.glmnet,alces.bio.sm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.sm.bio <- nicheOverlap(alces.bio.sm.raster.glmnet,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.sm.bio <- nicheOverlap(alces.bio.sm.raster.glmnet,alces.bio.sm.raster.svm,stat = 'D')

alces.marsxmaxent.sm.bio <- nicheOverlap(alces.bio.sm.raster.mars,alces.bio.sm.raster.maxent,stat = 'D')
alces.marsxmaxlike.sm.bio <- nicheOverlap(alces.bio.sm.raster.mars,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.marsxsvm.sm.bio <- nicheOverlap(alces.bio.sm.raster.mars,alces.bio.sm.raster.svm,stat = 'D')

alces.maxentxmaxlike.sm.bio <- nicheOverlap(alces.bio.sm.raster.maxent,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.maxentxsvm.sm.bio <- nicheOverlap(alces.bio.sm.raster.maxent,alces.bio.sm.raster.svm,stat = 'D')

alces.maxlikexsvm.sm.bio <- nicheOverlap(alces.bio.sm.raster.maxent,alces.bio.sm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.sm.bio,alces.rfxglm.sm.bio,alces.rfxglmnet.sm.bio,alces.rfxmars.sm.bio,alces.rfxmaxent.sm.bio,alces.rfxmaxlike.sm.bio,alces.rfxsvm.sm.bio)

alces.algD.bio.sm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.bio.sm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.bio.sm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.bio.sm[1,1] <- 1 ; alces.algD.bio.sm[2,1] <- alces.rfxbrt.sm.bio ; alces.algD.bio.sm[3,1] <- alces.rfxglm.sm.bio ; alces.algD.bio.sm[4,1] <- alces.rfxglmnet.sm.bio ; alces.algD.bio.sm[5,1] <- alces.rfxmars.sm.bio ; alces.algD.bio.sm[6,1] <- alces.rfxmaxent.sm.bio ; alces.algD.bio.sm[7,1] <- alces.rfxmaxlike.sm.bio ; alces.algD.bio.sm[8,1] <- alces.rfxsvm.sm.bio; alces.algD.bio.sm[2,2] <- 1 ; alces.algD.bio.sm[3,2] <- alces.brtxglm.sm.bio ; alces.algD.bio.sm[4,2] <- alces.brtxglmnet.sm.bio ; alces.algD.bio.sm[5,2] <- alces.brtxmars.sm.bio ; alces.algD.bio.sm[6,2] <- alces.brtxmaxent.sm.bio ; alces.algD.bio.sm[7,2] <- alces.brtxmaxlike.sm.bio ; alces.algD.bio.sm[8,2] <- alces.brtxsvm.sm.bio ; alces.algD.bio.sm[3,3] <- 1 ; alces.algD.bio.sm[4,3] <- alces.glmxglmnet.sm.bio ; alces.algD.bio.sm[5,3] <- alces.glmxmars.sm.bio ; alces.algD.bio.sm[6,3] <- alces.glmxmaxent.sm.bio ; alces.algD.bio.sm[7,3] <- alces.glmnetxmaxlike.sm.bio ; alces.algD.bio.sm[8,3] <- alces.glmnetxsvm.sm.bio ; alces.algD.bio.sm[4,4] <- 1 ; alces.algD.bio.sm[5,4] <- alces.glmnetxmars.sm.bio ; alces.algD.bio.sm[6,4] <- alces.glmnetxmaxent.sm.bio ; alces.algD.bio.sm[7,4] <- alces.glmnetxmaxlike.sm.bio ; alces.algD.bio.sm[8,4] <- alces.glmnetxsvm.sm.bio ; alces.algD.bio.sm[5,5] <- 1 ; alces.algD.bio.sm[6,5] <- alces.marsxmaxent.sm.bio ; alces.algD.bio.sm[7,5] <- alces.marsxmaxlike.sm.bio ; alces.algD.bio.sm[8,5] <- alces.marsxsvm.sm.bio ; alces.algD.bio.sm[6,6] <- 1 ; alces.algD.bio.sm[7,6] <- alces.maxentxmaxlike.sm.bio ; alces.algD.bio.sm[8,6] <- alces.maxentxsvm.sm.bio ; alces.algD.bio.sm[7,7] <- 1 ; alces.algD.bio.sm[8,7] <- alces.maxlikexsvm.sm.bio ; alces.algD.bio.sm[8,8] <- 1 ;

# Bio - mm
alces.rfxbrt.mm.bio <- nicheOverlap(alces.bio.mm.raster.rf,alces.bio.mm.raster.brt,stat = 'D')
alces.rfxglm.mm.bio <- nicheOverlap(alces.bio.mm.raster.rf,alces.envirem.mm.raster.glm,stat = 'D')
alces.rfxglmnet.mm.bio <- nicheOverlap(alces.bio.mm.raster.rf,alces.bio.mm.raster.glmnet,stat = 'D')
alces.rfxmars.mm.bio <- nicheOverlap(alces.bio.mm.raster.rf,alces.envirem.mm.raster.mars,stat = 'D')
alces.rfxmaxent.mm.bio <- nicheOverlap(alces.bio.mm.raster.rf,alces.bio.mm.raster.maxent,stat = 'D')
alces.rfxmaxlike.mm.bio <- nicheOverlap(alces.bio.mm.raster.rf,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.rfxsvm.mm.bio <- nicheOverlap(alces.bio.mm.raster.rf,alces.bio.mm.raster.svm,stat = 'D')

alces.brtxglm.mm.bio <- nicheOverlap(alces.bio.mm.raster.brt,alces.envirem.mm.raster.glm,stat = 'D')
alces.brtxglmnet.mm.bio <- nicheOverlap(alces.bio.mm.raster.brt,alces.bio.mm.raster.glmnet,stat = 'D')
alces.brtxmars.mm.bio <- nicheOverlap(alces.bio.mm.raster.brt,alces.envirem.mm.raster.mars,stat = 'D')
alces.brtxmaxent.mm.bio <- nicheOverlap(alces.bio.mm.raster.brt,alces.bio.mm.raster.maxent,stat = 'D')
alces.brtxmaxlike.mm.bio <- nicheOverlap(alces.bio.mm.raster.brt,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.brtxsvm.mm.bio <- nicheOverlap(alces.bio.mm.raster.brt,alces.bio.mm.raster.svm,stat = 'D')

alces.glmxglmnet.mm.bio <- nicheOverlap(alces.bio.mm.raster.glm,alces.bio.mm.raster.glmnet,stat = 'D')
alces.glmxmars.mm.bio <- nicheOverlap(alces.bio.mm.raster.glm,alces.envirem.mm.raster.mars,stat = 'D')
alces.glmxmaxent.mm.bio <- nicheOverlap(alces.bio.mm.raster.glm,alces.bio.mm.raster.maxent,stat = 'D')
alces.glmxmaxlike.mm.bio <- nicheOverlap(alces.bio.mm.raster.glm,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.glmxsvm.mm.bio <- nicheOverlap(alces.bio.mm.raster.glm,alces.bio.mm.raster.svm,stat = 'D')

alces.glmnetxmars.mm.bio <- nicheOverlap(alces.bio.mm.raster.glmnet,alces.envirem.mm.raster.mars,stat = 'D')
alces.glmnetxmaxent.mm.bio <- nicheOverlap(alces.bio.mm.raster.glmnet,alces.bio.mm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.mm.bio <- nicheOverlap(alces.bio.mm.raster.glmnet,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.mm.bio <- nicheOverlap(alces.bio.mm.raster.glmnet,alces.bio.mm.raster.svm,stat = 'D')

alces.marsxmaxent.mm.bio <- nicheOverlap(alces.bio.mm.raster.mars,alces.bio.mm.raster.maxent,stat = 'D')
alces.marsxmaxlike.mm.bio <- nicheOverlap(alces.bio.mm.raster.mars,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.marsxsvm.mm.bio <- nicheOverlap(alces.bio.mm.raster.mars,alces.bio.mm.raster.svm,stat = 'D')

alces.maxentxmaxlike.mm.bio <- nicheOverlap(alces.bio.mm.raster.maxent,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.maxentxsvm.mm.bio <- nicheOverlap(alces.bio.mm.raster.maxent,alces.bio.mm.raster.svm,stat = 'D')

alces.maxlikexsvm.mm.bio <- nicheOverlap(alces.bio.mm.raster.maxent,alces.bio.mm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.mm.bio,alces.rfxglm.mm.bio,alces.rfxglmnet.mm.bio,alces.rfxmars.mm.bio,alces.rfxmaxent.mm.bio,alces.rfxmaxlike.mm.bio,alces.rfxsvm.mm.bio)

alces.algD.bio.mm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.bio.mm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.bio.mm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.bio.mm[1,1] <- 1 ; alces.algD.bio.mm[2,1] <- alces.rfxbrt.mm.bio ; alces.algD.bio.mm[3,1] <- alces.rfxglm.mm.bio ; alces.algD.bio.mm[4,1] <- alces.rfxglmnet.mm.bio ; alces.algD.bio.mm[5,1] <- alces.rfxmars.mm.bio ; alces.algD.bio.mm[6,1] <- alces.rfxmaxent.mm.bio ; alces.algD.bio.mm[7,1] <- alces.rfxmaxlike.mm.bio ; alces.algD.bio.mm[8,1] <- alces.rfxsvm.mm.bio; alces.algD.bio.mm[2,2] <- 1 ; alces.algD.bio.mm[3,2] <- alces.brtxglm.mm.bio ; alces.algD.bio.mm[4,2] <- alces.brtxglmnet.mm.bio ; alces.algD.bio.mm[5,2] <- alces.brtxmars.mm.bio ; alces.algD.bio.mm[6,2] <- alces.brtxmaxent.mm.bio ; alces.algD.bio.mm[7,2] <- alces.brtxmaxlike.mm.bio ; alces.algD.bio.mm[8,2] <- alces.brtxsvm.mm.bio ; alces.algD.bio.mm[3,3] <- 1 ; alces.algD.bio.mm[4,3] <- alces.glmxglmnet.mm.bio ; alces.algD.bio.mm[5,3] <- alces.glmxmars.mm.bio ; alces.algD.bio.mm[6,3] <- alces.glmxmaxent.mm.bio ; alces.algD.bio.mm[7,3] <- alces.glmnetxmaxlike.mm.bio ; alces.algD.bio.mm[8,3] <- alces.glmnetxsvm.mm.bio ; alces.algD.bio.mm[4,4] <- 1 ; alces.algD.bio.mm[5,4] <- alces.glmnetxmars.mm.bio ; alces.algD.bio.mm[6,4] <- alces.glmnetxmaxent.mm.bio ; alces.algD.bio.mm[7,4] <- alces.glmnetxmaxlike.mm.bio ; alces.algD.bio.mm[8,4] <- alces.glmnetxsvm.mm.bio ; alces.algD.bio.mm[5,5] <- 1 ; alces.algD.bio.mm[6,5] <- alces.marsxmaxent.mm.bio ; alces.algD.bio.mm[7,5] <- alces.marsxmaxlike.mm.bio ; alces.algD.bio.mm[8,5] <- alces.marsxsvm.mm.bio ; alces.algD.bio.mm[6,6] <- 1 ; alces.algD.bio.mm[7,6] <- alces.maxentxmaxlike.mm.bio ; alces.algD.bio.mm[8,6] <- alces.maxentxsvm.mm.bio ; alces.algD.bio.mm[7,7] <- 1 ; alces.algD.bio.mm[8,7] <- alces.maxlikexsvm.mm.bio ; alces.algD.bio.mm[8,8] <- 1 ;


# Bio - bm
alces.rfxbrt.bm.bio <- nicheOverlap(alces.bio.bm.raster.rf,alces.bio.bm.raster.brt,stat = 'D')
alces.rfxglm.bm.bio <- nicheOverlap(alces.bio.bm.raster.rf,alces.envirem.bm.raster.glm,stat = 'D')
alces.rfxglmnet.bm.bio <- nicheOverlap(alces.bio.bm.raster.rf,alces.bio.bm.raster.glmnet,stat = 'D')
alces.rfxmars.bm.bio <- nicheOverlap(alces.bio.bm.raster.rf,alces.envirem.bm.raster.mars,stat = 'D')
alces.rfxmaxent.bm.bio <- nicheOverlap(alces.bio.bm.raster.rf,alces.bio.bm.raster.maxent,stat = 'D')
alces.rfxmaxlike.bm.bio <- nicheOverlap(alces.bio.bm.raster.rf,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.rfxsvm.bm.bio <- nicheOverlap(alces.bio.bm.raster.rf,alces.bio.bm.raster.svm,stat = 'D')

alces.brtxglm.bm.bio <- nicheOverlap(alces.bio.bm.raster.brt,alces.envirem.bm.raster.glm,stat = 'D')
alces.brtxglmnet.bm.bio <- nicheOverlap(alces.bio.bm.raster.brt,alces.bio.bm.raster.glmnet,stat = 'D')
alces.brtxmars.bm.bio <- nicheOverlap(alces.bio.bm.raster.brt,alces.envirem.bm.raster.mars,stat = 'D')
alces.brtxmaxent.bm.bio <- nicheOverlap(alces.bio.bm.raster.brt,alces.bio.bm.raster.maxent,stat = 'D')
alces.brtxmaxlike.bm.bio <- nicheOverlap(alces.bio.bm.raster.brt,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.brtxsvm.bm.bio <- nicheOverlap(alces.bio.bm.raster.brt,alces.bio.bm.raster.svm,stat = 'D')

alces.glmxglmnet.bm.bio <- nicheOverlap(alces.bio.bm.raster.glm,alces.bio.bm.raster.glmnet,stat = 'D')
alces.glmxmars.bm.bio <- nicheOverlap(alces.bio.bm.raster.glm,alces.envirem.bm.raster.mars,stat = 'D')
alces.glmxmaxent.bm.bio <- nicheOverlap(alces.bio.bm.raster.glm,alces.bio.bm.raster.maxent,stat = 'D')
alces.glmxmaxlike.bm.bio <- nicheOverlap(alces.bio.bm.raster.glm,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.glmxsvm.bm.bio <- nicheOverlap(alces.bio.bm.raster.glm,alces.bio.bm.raster.svm,stat = 'D')

alces.glmnetxmars.bm.bio <- nicheOverlap(alces.bio.bm.raster.glmnet,alces.envirem.bm.raster.mars,stat = 'D')
alces.glmnetxmaxent.bm.bio <- nicheOverlap(alces.bio.bm.raster.glmnet,alces.bio.bm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.bm.bio <- nicheOverlap(alces.bio.bm.raster.glmnet,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.bm.bio <- nicheOverlap(alces.bio.bm.raster.glmnet,alces.bio.bm.raster.svm,stat = 'D')

alces.marsxmaxent.bm.bio <- nicheOverlap(alces.bio.bm.raster.mars,alces.bio.bm.raster.maxent,stat = 'D')
alces.marsxmaxlike.bm.bio <- nicheOverlap(alces.bio.bm.raster.mars,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.marsxsvm.bm.bio <- nicheOverlap(alces.bio.bm.raster.mars,alces.bio.bm.raster.svm,stat = 'D')

alces.maxentxmaxlike.bm.bio <- nicheOverlap(alces.bio.bm.raster.maxent,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.maxentxsvm.bm.bio <- nicheOverlap(alces.bio.bm.raster.maxent,alces.bio.bm.raster.svm,stat = 'D')

alces.maxlikexsvm.bm.bio <- nicheOverlap(alces.bio.bm.raster.maxent,alces.bio.bm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.bm.bio,alces.rfxglm.bm.bio,alces.rfxglmnet.bm.bio,alces.rfxmars.bm.bio,alces.rfxmaxent.bm.bio,alces.rfxmaxlike.bm.bio,alces.rfxsvm.bm.bio)

alces.algD.bio.bm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.bio.bm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.bio.bm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.bio.bm[1,1] <- 1 ; alces.algD.bio.bm[2,1] <- alces.rfxbrt.bm.bio ; alces.algD.bio.bm[3,1] <- alces.rfxglm.bm.bio ; alces.algD.bio.bm[4,1] <- alces.rfxglmnet.bm.bio ; alces.algD.bio.bm[5,1] <- alces.rfxmars.bm.bio ; alces.algD.bio.bm[6,1] <- alces.rfxmaxent.bm.bio ; alces.algD.bio.bm[7,1] <- alces.rfxmaxlike.bm.bio ; alces.algD.bio.bm[8,1] <- alces.rfxsvm.bm.bio; alces.algD.bio.bm[2,2] <- 1 ; alces.algD.bio.bm[3,2] <- alces.brtxglm.bm.bio ; alces.algD.bio.bm[4,2] <- alces.brtxglmnet.bm.bio ; alces.algD.bio.bm[5,2] <- alces.brtxmars.bm.bio ; alces.algD.bio.bm[6,2] <- alces.brtxmaxent.bm.bio ; alces.algD.bio.bm[7,2] <- alces.brtxmaxlike.bm.bio ; alces.algD.bio.bm[8,2] <- alces.brtxsvm.bm.bio ; alces.algD.bio.bm[3,3] <- 1 ; alces.algD.bio.bm[4,3] <- alces.glmxglmnet.bm.bio ; alces.algD.bio.bm[5,3] <- alces.glmxmars.bm.bio ; alces.algD.bio.bm[6,3] <- alces.glmxmaxent.bm.bio ; alces.algD.bio.bm[7,3] <- alces.glmnetxmaxlike.bm.bio ; alces.algD.bio.bm[8,3] <- alces.glmnetxsvm.bm.bio ; alces.algD.bio.bm[4,4] <- 1 ; alces.algD.bio.bm[5,4] <- alces.glmnetxmars.bm.bio ; alces.algD.bio.bm[6,4] <- alces.glmnetxmaxent.bm.bio ; alces.algD.bio.bm[7,4] <- alces.glmnetxmaxlike.bm.bio ; alces.algD.bio.bm[8,4] <- alces.glmnetxsvm.bm.bio ; alces.algD.bio.bm[5,5] <- 1 ; alces.algD.bio.bm[6,5] <- alces.marsxmaxent.bm.bio ; alces.algD.bio.bm[7,5] <- alces.marsxmaxlike.bm.bio ; alces.algD.bio.bm[8,5] <- alces.marsxsvm.bm.bio ; alces.algD.bio.bm[6,6] <- 1 ; alces.algD.bio.bm[7,6] <- alces.maxentxmaxlike.bm.bio ; alces.algD.bio.bm[8,6] <- alces.maxentxsvm.bm.bio ; alces.algD.bio.bm[7,7] <- 1 ; alces.algD.bio.bm[8,7] <- alces.maxlikexsvm.bm.bio ; alces.algD.bio.bm[8,8] <- 1 ;

### Envirem ####
# envirem - SM
alces.rfxbrt.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.rf,alces.envirem.sm.raster.brt,stat = 'D')
alces.rfxglm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.rf,alces.envirem.sm.raster.glm,stat = 'D')
alces.rfxglmnet.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.rf,alces.envirem.sm.raster.glmnet,stat = 'D')
alces.rfxmars.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.rf,alces.envirem.sm.raster.mars,stat = 'D')
alces.rfxmaxent.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.rf,alces.envirem.sm.raster.maxent,stat = 'D')
alces.rfxmaxlike.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.rf,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.rfxsvm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.rf,alces.envirem.sm.raster.svm,stat = 'D')

alces.brtxglm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.brt,alces.envirem.sm.raster.glm,stat = 'D')
alces.brtxglmnet.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.brt,alces.envirem.sm.raster.glmnet,stat = 'D')
alces.brtxmars.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.brt,alces.envirem.sm.raster.mars,stat = 'D')
alces.brtxmaxent.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.brt,alces.envirem.sm.raster.maxent,stat = 'D')
alces.brtxmaxlike.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.brt,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.brtxsvm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.brt,alces.envirem.sm.raster.svm,stat = 'D')

alces.glmxglmnet.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glm,alces.envirem.sm.raster.glmnet,stat = 'D')
alces.glmxmars.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glm,alces.envirem.sm.raster.mars,stat = 'D')
alces.glmxmaxent.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glm,alces.envirem.sm.raster.maxent,stat = 'D')
alces.glmxmaxlike.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glm,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.glmxsvm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glm,alces.envirem.sm.raster.svm,stat = 'D')

alces.glmnetxmars.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glmnet,alces.envirem.sm.raster.mars,stat = 'D')
alces.glmnetxmaxent.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glmnet,alces.envirem.sm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glmnet,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.glmnet,alces.envirem.sm.raster.svm,stat = 'D')

alces.marsxmaxent.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.mars,alces.envirem.sm.raster.maxent,stat = 'D')
alces.marsxmaxlike.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.mars,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.marsxsvm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.mars,alces.envirem.sm.raster.svm,stat = 'D')

alces.maxentxmaxlike.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.maxent,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.maxentxsvm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.maxent,alces.envirem.sm.raster.svm,stat = 'D')

alces.maxlikexsvm.sm.envirem <- nicheOverlap(alces.envirem.sm.raster.maxent,alces.envirem.sm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.sm.envirem,alces.rfxglm.sm.envirem,alces.rfxglmnet.sm.envirem,alces.rfxmars.sm.envirem,alces.rfxmaxent.sm.envirem,alces.rfxmaxlike.sm.envirem,alces.rfxsvm.sm.envirem)

alces.algD.envirem.sm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.envirem.sm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.envirem.sm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.envirem.sm[1,1] <- 1 ; alces.algD.envirem.sm[2,1] <- alces.rfxbrt.sm.envirem ; alces.algD.envirem.sm[3,1] <- alces.rfxglm.sm.envirem ; alces.algD.envirem.sm[4,1] <- alces.rfxglmnet.sm.envirem ; alces.algD.envirem.sm[5,1] <- alces.rfxmars.sm.envirem ; alces.algD.envirem.sm[6,1] <- alces.rfxmaxent.sm.envirem ; alces.algD.envirem.sm[7,1] <- alces.rfxmaxlike.sm.envirem ; alces.algD.envirem.sm[8,1] <- alces.rfxsvm.sm.envirem; alces.algD.envirem.sm[2,2] <- 1 ; alces.algD.envirem.sm[3,2] <- alces.brtxglm.sm.envirem ; alces.algD.envirem.sm[4,2] <- alces.brtxglmnet.sm.envirem ; alces.algD.envirem.sm[5,2] <- alces.brtxmars.sm.envirem ; alces.algD.envirem.sm[6,2] <- alces.brtxmaxent.sm.envirem ; alces.algD.envirem.sm[7,2] <- alces.brtxmaxlike.sm.envirem ; alces.algD.envirem.sm[8,2] <- alces.brtxsvm.sm.envirem ; alces.algD.envirem.sm[3,3] <- 1 ; alces.algD.envirem.sm[4,3] <- alces.glmxglmnet.sm.envirem ; alces.algD.envirem.sm[5,3] <- alces.glmxmars.sm.envirem ; alces.algD.envirem.sm[6,3] <- alces.glmxmaxent.sm.envirem ; alces.algD.envirem.sm[7,3] <- alces.glmnetxmaxlike.sm.envirem ; alces.algD.envirem.sm[8,3] <- alces.glmnetxsvm.sm.envirem ; alces.algD.envirem.sm[4,4] <- 1 ; alces.algD.envirem.sm[5,4] <- alces.glmnetxmars.sm.envirem ; alces.algD.envirem.sm[6,4] <- alces.glmnetxmaxent.sm.envirem ; alces.algD.envirem.sm[7,4] <- alces.glmnetxmaxlike.sm.envirem ; alces.algD.envirem.sm[8,4] <- alces.glmnetxsvm.sm.envirem ; alces.algD.envirem.sm[5,5] <- 1 ; alces.algD.envirem.sm[6,5] <- alces.marsxmaxent.sm.envirem ; alces.algD.envirem.sm[7,5] <- alces.marsxmaxlike.sm.envirem ; alces.algD.envirem.sm[8,5] <- alces.marsxsvm.sm.envirem ; alces.algD.envirem.sm[6,6] <- 1 ; alces.algD.envirem.sm[7,6] <- alces.maxentxmaxlike.sm.envirem ; alces.algD.envirem.sm[8,6] <- alces.maxentxsvm.sm.envirem ; alces.algD.envirem.sm[7,7] <- 1 ; alces.algD.envirem.sm[8,7] <- alces.maxlikexsvm.sm.envirem ; alces.algD.envirem.sm[8,8] <- 1 ;

# envirem - mm
alces.rfxbrt.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.rf,alces.envirem.mm.raster.brt,stat = 'D')
alces.rfxglm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.rf,alces.envirem.mm.raster.glm,stat = 'D')
alces.rfxglmnet.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.rf,alces.envirem.mm.raster.glmnet,stat = 'D')
alces.rfxmars.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.rf,alces.envirem.mm.raster.mars,stat = 'D')
alces.rfxmaxent.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.rf,alces.envirem.mm.raster.maxent,stat = 'D')
alces.rfxmaxlike.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.rf,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.rfxsvm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.rf,alces.envirem.mm.raster.svm,stat = 'D')

alces.brtxglm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.brt,alces.envirem.mm.raster.glm,stat = 'D')
alces.brtxglmnet.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.brt,alces.envirem.mm.raster.glmnet,stat = 'D')
alces.brtxmars.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.brt,alces.envirem.mm.raster.mars,stat = 'D')
alces.brtxmaxent.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.brt,alces.envirem.mm.raster.maxent,stat = 'D')
alces.brtxmaxlike.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.brt,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.brtxsvm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.brt,alces.envirem.mm.raster.svm,stat = 'D')

alces.glmxglmnet.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glm,alces.envirem.mm.raster.glmnet,stat = 'D')
alces.glmxmars.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glm,alces.envirem.mm.raster.mars,stat = 'D')
alces.glmxmaxent.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glm,alces.envirem.mm.raster.maxent,stat = 'D')
alces.glmxmaxlike.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glm,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.glmxsvm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glm,alces.envirem.mm.raster.svm,stat = 'D')

alces.glmnetxmars.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glmnet,alces.envirem.mm.raster.mars,stat = 'D')
alces.glmnetxmaxent.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glmnet,alces.envirem.mm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glmnet,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.glmnet,alces.envirem.mm.raster.svm,stat = 'D')

alces.marsxmaxent.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.mars,alces.envirem.mm.raster.maxent,stat = 'D')
alces.marsxmaxlike.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.mars,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.marsxsvm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.mars,alces.envirem.mm.raster.svm,stat = 'D')

alces.maxentxmaxlike.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.maxent,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.maxentxsvm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.maxent,alces.envirem.mm.raster.svm,stat = 'D')

alces.maxlikexsvm.mm.envirem <- nicheOverlap(alces.envirem.mm.raster.maxent,alces.envirem.mm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.mm.envirem,alces.rfxglm.mm.envirem,alces.rfxglmnet.mm.envirem,alces.rfxmars.mm.envirem,alces.rfxmaxent.mm.envirem,alces.rfxmaxlike.mm.envirem,alces.rfxsvm.mm.envirem)

alces.algD.envirem.mm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.envirem.mm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.envirem.mm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.envirem.mm[1,1] <- 1 ; alces.algD.envirem.mm[2,1] <- alces.rfxbrt.mm.envirem ; alces.algD.envirem.mm[3,1] <- alces.rfxglm.mm.envirem ; alces.algD.envirem.mm[4,1] <- alces.rfxglmnet.mm.envirem ; alces.algD.envirem.mm[5,1] <- alces.rfxmars.mm.envirem ; alces.algD.envirem.mm[6,1] <- alces.rfxmaxent.mm.envirem ; alces.algD.envirem.mm[7,1] <- alces.rfxmaxlike.mm.envirem ; alces.algD.envirem.mm[8,1] <- alces.rfxsvm.mm.envirem; alces.algD.envirem.mm[2,2] <- 1 ; alces.algD.envirem.mm[3,2] <- alces.brtxglm.mm.envirem ; alces.algD.envirem.mm[4,2] <- alces.brtxglmnet.mm.envirem ; alces.algD.envirem.mm[5,2] <- alces.brtxmars.mm.envirem ; alces.algD.envirem.mm[6,2] <- alces.brtxmaxent.mm.envirem ; alces.algD.envirem.mm[7,2] <- alces.brtxmaxlike.mm.envirem ; alces.algD.envirem.mm[8,2] <- alces.brtxsvm.mm.envirem ; alces.algD.envirem.mm[3,3] <- 1 ; alces.algD.envirem.mm[4,3] <- alces.glmxglmnet.mm.envirem ; alces.algD.envirem.mm[5,3] <- alces.glmxmars.mm.envirem ; alces.algD.envirem.mm[6,3] <- alces.glmxmaxent.mm.envirem ; alces.algD.envirem.mm[7,3] <- alces.glmnetxmaxlike.mm.envirem ; alces.algD.envirem.mm[8,3] <- alces.glmnetxsvm.mm.envirem ; alces.algD.envirem.mm[4,4] <- 1 ; alces.algD.envirem.mm[5,4] <- alces.glmnetxmars.mm.envirem ; alces.algD.envirem.mm[6,4] <- alces.glmnetxmaxent.mm.envirem ; alces.algD.envirem.mm[7,4] <- alces.glmnetxmaxlike.mm.envirem ; alces.algD.envirem.mm[8,4] <- alces.glmnetxsvm.mm.envirem ; alces.algD.envirem.mm[5,5] <- 1 ; alces.algD.envirem.mm[6,5] <- alces.marsxmaxent.mm.envirem ; alces.algD.envirem.mm[7,5] <- alces.marsxmaxlike.mm.envirem ; alces.algD.envirem.mm[8,5] <- alces.marsxsvm.mm.envirem ; alces.algD.envirem.mm[6,6] <- 1 ; alces.algD.envirem.mm[7,6] <- alces.maxentxmaxlike.mm.envirem ; alces.algD.envirem.mm[8,6] <- alces.maxentxsvm.mm.envirem ; alces.algD.envirem.mm[7,7] <- 1 ; alces.algD.envirem.mm[8,7] <- alces.maxlikexsvm.mm.envirem ; alces.algD.envirem.mm[8,8] <- 1 ;


# envirem - bm
alces.rfxbrt.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.rf,alces.envirem.bm.raster.brt,stat = 'D')
alces.rfxglm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.rf,alces.envirem.bm.raster.glm,stat = 'D')
alces.rfxglmnet.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.rf,alces.envirem.bm.raster.glmnet,stat = 'D')
alces.rfxmars.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.rf,alces.envirem.bm.raster.mars,stat = 'D')
alces.rfxmaxent.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.rf,alces.envirem.bm.raster.maxent,stat = 'D')
alces.rfxmaxlike.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.rf,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.rfxsvm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.rf,alces.envirem.bm.raster.svm,stat = 'D')

alces.brtxglm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.brt,alces.envirem.bm.raster.glm,stat = 'D')
alces.brtxglmnet.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.brt,alces.envirem.bm.raster.glmnet,stat = 'D')
alces.brtxmars.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.brt,alces.envirem.bm.raster.mars,stat = 'D')
alces.brtxmaxent.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.brt,alces.envirem.bm.raster.maxent,stat = 'D')
alces.brtxmaxlike.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.brt,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.brtxsvm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.brt,alces.envirem.bm.raster.svm,stat = 'D')

alces.glmxglmnet.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glm,alces.envirem.bm.raster.glmnet,stat = 'D')
alces.glmxmars.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glm,alces.envirem.bm.raster.mars,stat = 'D')
alces.glmxmaxent.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glm,alces.envirem.bm.raster.maxent,stat = 'D')
alces.glmxmaxlike.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glm,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.glmxsvm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glm,alces.envirem.bm.raster.svm,stat = 'D')

alces.glmnetxmars.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glmnet,alces.envirem.bm.raster.mars,stat = 'D')
alces.glmnetxmaxent.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glmnet,alces.envirem.bm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glmnet,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.glmnet,alces.envirem.bm.raster.svm,stat = 'D')

alces.marsxmaxent.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.mars,alces.envirem.bm.raster.maxent,stat = 'D')
alces.marsxmaxlike.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.mars,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.marsxsvm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.mars,alces.envirem.bm.raster.svm,stat = 'D')

alces.maxentxmaxlike.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.maxent,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.maxentxsvm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.maxent,alces.envirem.bm.raster.svm,stat = 'D')

alces.maxlikexsvm.bm.envirem <- nicheOverlap(alces.envirem.bm.raster.maxent,alces.envirem.bm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.bm.envirem,alces.rfxglm.bm.envirem,alces.rfxglmnet.bm.envirem,alces.rfxmars.bm.envirem,alces.rfxmaxent.bm.envirem,alces.rfxmaxlike.bm.envirem,alces.rfxsvm.bm.envirem)

alces.algD.envirem.bm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.envirem.bm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.envirem.bm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.envirem.bm[1,1] <- 1 ; alces.algD.envirem.bm[2,1] <- alces.rfxbrt.bm.envirem ; alces.algD.envirem.bm[3,1] <- alces.rfxglm.bm.envirem ; alces.algD.envirem.bm[4,1] <- alces.rfxglmnet.bm.envirem ; alces.algD.envirem.bm[5,1] <- alces.rfxmars.bm.envirem ; alces.algD.envirem.bm[6,1] <- alces.rfxmaxent.bm.envirem ; alces.algD.envirem.bm[7,1] <- alces.rfxmaxlike.bm.envirem ; alces.algD.envirem.bm[8,1] <- alces.rfxsvm.bm.envirem; alces.algD.envirem.bm[2,2] <- 1 ; alces.algD.envirem.bm[3,2] <- alces.brtxglm.bm.envirem ; alces.algD.envirem.bm[4,2] <- alces.brtxglmnet.bm.envirem ; alces.algD.envirem.bm[5,2] <- alces.brtxmars.bm.envirem ; alces.algD.envirem.bm[6,2] <- alces.brtxmaxent.bm.envirem ; alces.algD.envirem.bm[7,2] <- alces.brtxmaxlike.bm.envirem ; alces.algD.envirem.bm[8,2] <- alces.brtxsvm.bm.envirem ; alces.algD.envirem.bm[3,3] <- 1 ; alces.algD.envirem.bm[4,3] <- alces.glmxglmnet.bm.envirem ; alces.algD.envirem.bm[5,3] <- alces.glmxmars.bm.envirem ; alces.algD.envirem.bm[6,3] <- alces.glmxmaxent.bm.envirem ; alces.algD.envirem.bm[7,3] <- alces.glmnetxmaxlike.bm.envirem ; alces.algD.envirem.bm[8,3] <- alces.glmnetxsvm.bm.envirem ; alces.algD.envirem.bm[4,4] <- 1 ; alces.algD.envirem.bm[5,4] <- alces.glmnetxmars.bm.envirem ; alces.algD.envirem.bm[6,4] <- alces.glmnetxmaxent.bm.envirem ; alces.algD.envirem.bm[7,4] <- alces.glmnetxmaxlike.bm.envirem ; alces.algD.envirem.bm[8,4] <- alces.glmnetxsvm.bm.envirem ; alces.algD.envirem.bm[5,5] <- 1 ; alces.algD.envirem.bm[6,5] <- alces.marsxmaxent.bm.envirem ; alces.algD.envirem.bm[7,5] <- alces.marsxmaxlike.bm.envirem ; alces.algD.envirem.bm[8,5] <- alces.marsxsvm.bm.envirem ; alces.algD.envirem.bm[6,6] <- 1 ; alces.algD.envirem.bm[7,6] <- alces.maxentxmaxlike.bm.envirem ; alces.algD.envirem.bm[8,6] <- alces.maxentxsvm.bm.envirem ; alces.algD.envirem.bm[7,7] <- 1 ; alces.algD.envirem.bm[8,7] <- alces.maxlikexsvm.bm.envirem ; alces.algD.envirem.bm[8,8] <- 1 ;

### Merraclim #####
# merra - SM
alces.rfxbrt.sm.merra <- nicheOverlap(alces.merra.sm.raster.rf,alces.merra.sm.raster.brt,stat = 'D')
alces.rfxglm.sm.merra <- nicheOverlap(alces.merra.sm.raster.rf,alces.envirem.sm.raster.glm,stat = 'D')
alces.rfxglmnet.sm.merra <- nicheOverlap(alces.merra.sm.raster.rf,alces.merra.sm.raster.glmnet,stat = 'D')
alces.rfxmars.sm.merra <- nicheOverlap(alces.merra.sm.raster.rf,alces.envirem.sm.raster.mars,stat = 'D')
alces.rfxmaxent.sm.merra <- nicheOverlap(alces.merra.sm.raster.rf,alces.merra.sm.raster.maxent,stat = 'D')
alces.rfxmaxlike.sm.merra <- nicheOverlap(alces.merra.sm.raster.rf,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.rfxsvm.sm.merra <- nicheOverlap(alces.merra.sm.raster.rf,alces.merra.sm.raster.svm,stat = 'D')

alces.brtxglm.sm.merra <- nicheOverlap(alces.merra.sm.raster.brt,alces.envirem.sm.raster.glm,stat = 'D')
alces.brtxglmnet.sm.merra <- nicheOverlap(alces.merra.sm.raster.brt,alces.merra.sm.raster.glmnet,stat = 'D')
alces.brtxmars.sm.merra <- nicheOverlap(alces.merra.sm.raster.brt,alces.envirem.sm.raster.mars,stat = 'D')
alces.brtxmaxent.sm.merra <- nicheOverlap(alces.merra.sm.raster.brt,alces.merra.sm.raster.maxent,stat = 'D')
alces.brtxmaxlike.sm.merra <- nicheOverlap(alces.merra.sm.raster.brt,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.brtxsvm.sm.merra <- nicheOverlap(alces.merra.sm.raster.brt,alces.merra.sm.raster.svm,stat = 'D')

alces.glmxglmnet.sm.merra <- nicheOverlap(alces.merra.sm.raster.glm,alces.merra.sm.raster.glmnet,stat = 'D')
alces.glmxmars.sm.merra <- nicheOverlap(alces.merra.sm.raster.glm,alces.envirem.sm.raster.mars,stat = 'D')
alces.glmxmaxent.sm.merra <- nicheOverlap(alces.merra.sm.raster.glm,alces.merra.sm.raster.maxent,stat = 'D')
alces.glmxmaxlike.sm.merra <- nicheOverlap(alces.merra.sm.raster.glm,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.glmxsvm.sm.merra <- nicheOverlap(alces.merra.sm.raster.glm,alces.merra.sm.raster.svm,stat = 'D')

alces.glmnetxmars.sm.merra <- nicheOverlap(alces.merra.sm.raster.glmnet,alces.envirem.sm.raster.mars,stat = 'D')
alces.glmnetxmaxent.sm.merra <- nicheOverlap(alces.merra.sm.raster.glmnet,alces.merra.sm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.sm.merra <- nicheOverlap(alces.merra.sm.raster.glmnet,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.sm.merra <- nicheOverlap(alces.merra.sm.raster.glmnet,alces.merra.sm.raster.svm,stat = 'D')

alces.marsxmaxent.sm.merra <- nicheOverlap(alces.merra.sm.raster.mars,alces.merra.sm.raster.maxent,stat = 'D')
alces.marsxmaxlike.sm.merra <- nicheOverlap(alces.merra.sm.raster.mars,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.marsxsvm.sm.merra <- nicheOverlap(alces.merra.sm.raster.mars,alces.merra.sm.raster.svm,stat = 'D')

alces.maxentxmaxlike.sm.merra <- nicheOverlap(alces.merra.sm.raster.maxent,alces.envirem.sm.raster.maxlike,stat = 'D')
alces.maxentxsvm.sm.merra <- nicheOverlap(alces.merra.sm.raster.maxent,alces.merra.sm.raster.svm,stat = 'D')

alces.maxlikexsvm.sm.merra <- nicheOverlap(alces.merra.sm.raster.maxent,alces.merra.sm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.sm.merra,alces.rfxglm.sm.merra,alces.rfxglmnet.sm.merra,alces.rfxmars.sm.merra,alces.rfxmaxent.sm.merra,alces.rfxmaxlike.sm.merra,alces.rfxsvm.sm.merra)

alces.algD.merra.sm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.merra.sm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.merra.sm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.merra.sm[1,1] <- 1 ; alces.algD.merra.sm[2,1] <- alces.rfxbrt.sm.merra ; alces.algD.merra.sm[3,1] <- alces.rfxglm.sm.merra ; alces.algD.merra.sm[4,1] <- alces.rfxglmnet.sm.merra ; alces.algD.merra.sm[5,1] <- alces.rfxmars.sm.merra ; alces.algD.merra.sm[6,1] <- alces.rfxmaxent.sm.merra ; alces.algD.merra.sm[7,1] <- alces.rfxmaxlike.sm.merra ; alces.algD.merra.sm[8,1] <- alces.rfxsvm.sm.merra; alces.algD.merra.sm[2,2] <- 1 ; alces.algD.merra.sm[3,2] <- alces.brtxglm.sm.merra ; alces.algD.merra.sm[4,2] <- alces.brtxglmnet.sm.merra ; alces.algD.merra.sm[5,2] <- alces.brtxmars.sm.merra ; alces.algD.merra.sm[6,2] <- alces.brtxmaxent.sm.merra ; alces.algD.merra.sm[7,2] <- alces.brtxmaxlike.sm.merra ; alces.algD.merra.sm[8,2] <- alces.brtxsvm.sm.merra ; alces.algD.merra.sm[3,3] <- 1 ; alces.algD.merra.sm[4,3] <- alces.glmxglmnet.sm.merra ; alces.algD.merra.sm[5,3] <- alces.glmxmars.sm.merra ; alces.algD.merra.sm[6,3] <- alces.glmxmaxent.sm.merra ; alces.algD.merra.sm[7,3] <- alces.glmnetxmaxlike.sm.merra ; alces.algD.merra.sm[8,3] <- alces.glmnetxsvm.sm.merra ; alces.algD.merra.sm[4,4] <- 1 ; alces.algD.merra.sm[5,4] <- alces.glmnetxmars.sm.merra ; alces.algD.merra.sm[6,4] <- alces.glmnetxmaxent.sm.merra ; alces.algD.merra.sm[7,4] <- alces.glmnetxmaxlike.sm.merra ; alces.algD.merra.sm[8,4] <- alces.glmnetxsvm.sm.merra ; alces.algD.merra.sm[5,5] <- 1 ; alces.algD.merra.sm[6,5] <- alces.marsxmaxent.sm.merra ; alces.algD.merra.sm[7,5] <- alces.marsxmaxlike.sm.merra ; alces.algD.merra.sm[8,5] <- alces.marsxsvm.sm.merra ; alces.algD.merra.sm[6,6] <- 1 ; alces.algD.merra.sm[7,6] <- alces.maxentxmaxlike.sm.merra ; alces.algD.merra.sm[8,6] <- alces.maxentxsvm.sm.merra ; alces.algD.merra.sm[7,7] <- 1 ; alces.algD.merra.sm[8,7] <- alces.maxlikexsvm.sm.merra ; alces.algD.merra.sm[8,8] <- 1 ;

# merra - mm
alces.rfxbrt.mm.merra <- nicheOverlap(alces.merra.mm.raster.rf,alces.merra.mm.raster.brt,stat = 'D')
alces.rfxglm.mm.merra <- nicheOverlap(alces.merra.mm.raster.rf,alces.envirem.mm.raster.glm,stat = 'D')
alces.rfxglmnet.mm.merra <- nicheOverlap(alces.merra.mm.raster.rf,alces.merra.mm.raster.glmnet,stat = 'D')
alces.rfxmars.mm.merra <- nicheOverlap(alces.merra.mm.raster.rf,alces.envirem.mm.raster.mars,stat = 'D')
alces.rfxmaxent.mm.merra <- nicheOverlap(alces.merra.mm.raster.rf,alces.merra.mm.raster.maxent,stat = 'D')
alces.rfxmaxlike.mm.merra <- nicheOverlap(alces.merra.mm.raster.rf,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.rfxsvm.mm.merra <- nicheOverlap(alces.merra.mm.raster.rf,alces.merra.mm.raster.svm,stat = 'D')

alces.brtxglm.mm.merra <- nicheOverlap(alces.merra.mm.raster.brt,alces.envirem.mm.raster.glm,stat = 'D')
alces.brtxglmnet.mm.merra <- nicheOverlap(alces.merra.mm.raster.brt,alces.merra.mm.raster.glmnet,stat = 'D')
alces.brtxmars.mm.merra <- nicheOverlap(alces.merra.mm.raster.brt,alces.envirem.mm.raster.mars,stat = 'D')
alces.brtxmaxent.mm.merra <- nicheOverlap(alces.merra.mm.raster.brt,alces.merra.mm.raster.maxent,stat = 'D')
alces.brtxmaxlike.mm.merra <- nicheOverlap(alces.merra.mm.raster.brt,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.brtxsvm.mm.merra <- nicheOverlap(alces.merra.mm.raster.brt,alces.merra.mm.raster.svm,stat = 'D')

alces.glmxglmnet.mm.merra <- nicheOverlap(alces.merra.mm.raster.glm,alces.merra.mm.raster.glmnet,stat = 'D')
alces.glmxmars.mm.merra <- nicheOverlap(alces.merra.mm.raster.glm,alces.envirem.mm.raster.mars,stat = 'D')
alces.glmxmaxent.mm.merra <- nicheOverlap(alces.merra.mm.raster.glm,alces.merra.mm.raster.maxent,stat = 'D')
alces.glmxmaxlike.mm.merra <- nicheOverlap(alces.merra.mm.raster.glm,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.glmxsvm.mm.merra <- nicheOverlap(alces.merra.mm.raster.glm,alces.merra.mm.raster.svm,stat = 'D')

alces.glmnetxmars.mm.merra <- nicheOverlap(alces.merra.mm.raster.glmnet,alces.envirem.mm.raster.mars,stat = 'D')
alces.glmnetxmaxent.mm.merra <- nicheOverlap(alces.merra.mm.raster.glmnet,alces.merra.mm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.mm.merra <- nicheOverlap(alces.merra.mm.raster.glmnet,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.mm.merra <- nicheOverlap(alces.merra.mm.raster.glmnet,alces.merra.mm.raster.svm,stat = 'D')

alces.marsxmaxent.mm.merra <- nicheOverlap(alces.merra.mm.raster.mars,alces.merra.mm.raster.maxent,stat = 'D')
alces.marsxmaxlike.mm.merra <- nicheOverlap(alces.merra.mm.raster.mars,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.marsxsvm.mm.merra <- nicheOverlap(alces.merra.mm.raster.mars,alces.merra.mm.raster.svm,stat = 'D')

alces.maxentxmaxlike.mm.merra <- nicheOverlap(alces.merra.mm.raster.maxent,alces.envirem.mm.raster.maxlike,stat = 'D')
alces.maxentxsvm.mm.merra <- nicheOverlap(alces.merra.mm.raster.maxent,alces.merra.mm.raster.svm,stat = 'D')

alces.maxlikexsvm.mm.merra <- nicheOverlap(alces.merra.mm.raster.maxent,alces.merra.mm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.mm.merra,alces.rfxglm.mm.merra,alces.rfxglmnet.mm.merra,alces.rfxmars.mm.merra,alces.rfxmaxent.mm.merra,alces.rfxmaxlike.mm.merra,alces.rfxsvm.mm.merra)

alces.algD.merra.mm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.merra.mm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.merra.mm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.merra.mm[1,1] <- 1 ; alces.algD.merra.mm[2,1] <- alces.rfxbrt.mm.merra ; alces.algD.merra.mm[3,1] <- alces.rfxglm.mm.merra ; alces.algD.merra.mm[4,1] <- alces.rfxglmnet.mm.merra ; alces.algD.merra.mm[5,1] <- alces.rfxmars.mm.merra ; alces.algD.merra.mm[6,1] <- alces.rfxmaxent.mm.merra ; alces.algD.merra.mm[7,1] <- alces.rfxmaxlike.mm.merra ; alces.algD.merra.mm[8,1] <- alces.rfxsvm.mm.merra; alces.algD.merra.mm[2,2] <- 1 ; alces.algD.merra.mm[3,2] <- alces.brtxglm.mm.merra ; alces.algD.merra.mm[4,2] <- alces.brtxglmnet.mm.merra ; alces.algD.merra.mm[5,2] <- alces.brtxmars.mm.merra ; alces.algD.merra.mm[6,2] <- alces.brtxmaxent.mm.merra ; alces.algD.merra.mm[7,2] <- alces.brtxmaxlike.mm.merra ; alces.algD.merra.mm[8,2] <- alces.brtxsvm.mm.merra ; alces.algD.merra.mm[3,3] <- 1 ; alces.algD.merra.mm[4,3] <- alces.glmxglmnet.mm.merra ; alces.algD.merra.mm[5,3] <- alces.glmxmars.mm.merra ; alces.algD.merra.mm[6,3] <- alces.glmxmaxent.mm.merra ; alces.algD.merra.mm[7,3] <- alces.glmnetxmaxlike.mm.merra ; alces.algD.merra.mm[8,3] <- alces.glmnetxsvm.mm.merra ; alces.algD.merra.mm[4,4] <- 1 ; alces.algD.merra.mm[5,4] <- alces.glmnetxmars.mm.merra ; alces.algD.merra.mm[6,4] <- alces.glmnetxmaxent.mm.merra ; alces.algD.merra.mm[7,4] <- alces.glmnetxmaxlike.mm.merra ; alces.algD.merra.mm[8,4] <- alces.glmnetxsvm.mm.merra ; alces.algD.merra.mm[5,5] <- 1 ; alces.algD.merra.mm[6,5] <- alces.marsxmaxent.mm.merra ; alces.algD.merra.mm[7,5] <- alces.marsxmaxlike.mm.merra ; alces.algD.merra.mm[8,5] <- alces.marsxsvm.mm.merra ; alces.algD.merra.mm[6,6] <- 1 ; alces.algD.merra.mm[7,6] <- alces.maxentxmaxlike.mm.merra ; alces.algD.merra.mm[8,6] <- alces.maxentxsvm.mm.merra ; alces.algD.merra.mm[7,7] <- 1 ; alces.algD.merra.mm[8,7] <- alces.maxlikexsvm.mm.merra ; alces.algD.merra.mm[8,8] <- 1 ;


# merra - bm
alces.rfxbrt.bm.merra <- nicheOverlap(alces.merra.bm.raster.rf,alces.merra.bm.raster.brt,stat = 'D')
alces.rfxglm.bm.merra <- nicheOverlap(alces.merra.bm.raster.rf,alces.envirem.bm.raster.glm,stat = 'D')
alces.rfxglmnet.bm.merra <- nicheOverlap(alces.merra.bm.raster.rf,alces.merra.bm.raster.glmnet,stat = 'D')
alces.rfxmars.bm.merra <- nicheOverlap(alces.merra.bm.raster.rf,alces.envirem.bm.raster.mars,stat = 'D')
alces.rfxmaxent.bm.merra <- nicheOverlap(alces.merra.bm.raster.rf,alces.merra.bm.raster.maxent,stat = 'D')
alces.rfxmaxlike.bm.merra <- nicheOverlap(alces.merra.bm.raster.rf,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.rfxsvm.bm.merra <- nicheOverlap(alces.merra.bm.raster.rf,alces.merra.bm.raster.svm,stat = 'D')

alces.brtxglm.bm.merra <- nicheOverlap(alces.merra.bm.raster.brt,alces.envirem.bm.raster.glm,stat = 'D')
alces.brtxglmnet.bm.merra <- nicheOverlap(alces.merra.bm.raster.brt,alces.merra.bm.raster.glmnet,stat = 'D')
alces.brtxmars.bm.merra <- nicheOverlap(alces.merra.bm.raster.brt,alces.envirem.bm.raster.mars,stat = 'D')
alces.brtxmaxent.bm.merra <- nicheOverlap(alces.merra.bm.raster.brt,alces.merra.bm.raster.maxent,stat = 'D')
alces.brtxmaxlike.bm.merra <- nicheOverlap(alces.merra.bm.raster.brt,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.brtxsvm.bm.merra <- nicheOverlap(alces.merra.bm.raster.brt,alces.merra.bm.raster.svm,stat = 'D')

alces.glmxglmnet.bm.merra <- nicheOverlap(alces.merra.bm.raster.glm,alces.merra.bm.raster.glmnet,stat = 'D')
alces.glmxmars.bm.merra <- nicheOverlap(alces.merra.bm.raster.glm,alces.envirem.bm.raster.mars,stat = 'D')
alces.glmxmaxent.bm.merra <- nicheOverlap(alces.merra.bm.raster.glm,alces.merra.bm.raster.maxent,stat = 'D')
alces.glmxmaxlike.bm.merra <- nicheOverlap(alces.merra.bm.raster.glm,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.glmxsvm.bm.merra <- nicheOverlap(alces.merra.bm.raster.glm,alces.merra.bm.raster.svm,stat = 'D')

alces.glmnetxmars.bm.merra <- nicheOverlap(alces.merra.bm.raster.glmnet,alces.envirem.bm.raster.mars,stat = 'D')
alces.glmnetxmaxent.bm.merra <- nicheOverlap(alces.merra.bm.raster.glmnet,alces.merra.bm.raster.maxent,stat = 'D')
alces.glmnetxmaxlike.bm.merra <- nicheOverlap(alces.merra.bm.raster.glmnet,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.glmnetxsvm.bm.merra <- nicheOverlap(alces.merra.bm.raster.glmnet,alces.merra.bm.raster.svm,stat = 'D')

alces.marsxmaxent.bm.merra <- nicheOverlap(alces.merra.bm.raster.mars,alces.merra.bm.raster.maxent,stat = 'D')
alces.marsxmaxlike.bm.merra <- nicheOverlap(alces.merra.bm.raster.mars,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.marsxsvm.bm.merra <- nicheOverlap(alces.merra.bm.raster.mars,alces.merra.bm.raster.svm,stat = 'D')

alces.maxentxmaxlike.bm.merra <- nicheOverlap(alces.merra.bm.raster.maxent,alces.envirem.bm.raster.maxlike,stat = 'D')
alces.maxentxsvm.bm.merra <- nicheOverlap(alces.merra.bm.raster.maxent,alces.merra.bm.raster.svm,stat = 'D')

alces.maxlikexsvm.bm.merra <- nicheOverlap(alces.merra.bm.raster.maxent,alces.merra.bm.raster.svm,stat = 'D')

alces.algD.rf <- c(alces.rfxbrt.bm.merra,alces.rfxglm.bm.merra,alces.rfxglmnet.bm.merra,alces.rfxmars.bm.merra,alces.rfxmaxent.bm.merra,alces.rfxmaxlike.bm.merra,alces.rfxsvm.bm.merra)

alces.algD.merra.bm <- matrix(ncol=8,nrow = 8)
rownames(alces.algD.merra.bm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
colnames(alces.algD.merra.bm) <- c('RF','BRT','GLM','GLMNet','MARS','MaxEnt','MaxLike','SVM')
alces.algD.merra.bm[1,1] <- 1 ; alces.algD.merra.bm[2,1] <- alces.rfxbrt.bm.merra ; alces.algD.merra.bm[3,1] <- alces.rfxglm.bm.merra ; alces.algD.merra.bm[4,1] <- alces.rfxglmnet.bm.merra ; alces.algD.merra.bm[5,1] <- alces.rfxmars.bm.merra ; alces.algD.merra.bm[6,1] <- alces.rfxmaxent.bm.merra ; alces.algD.merra.bm[7,1] <- alces.rfxmaxlike.bm.merra ; alces.algD.merra.bm[8,1] <- alces.rfxsvm.bm.merra; alces.algD.merra.bm[2,2] <- 1 ; alces.algD.merra.bm[3,2] <- alces.brtxglm.bm.merra ; alces.algD.merra.bm[4,2] <- alces.brtxglmnet.bm.merra ; alces.algD.merra.bm[5,2] <- alces.brtxmars.bm.merra ; alces.algD.merra.bm[6,2] <- alces.brtxmaxent.bm.merra ; alces.algD.merra.bm[7,2] <- alces.brtxmaxlike.bm.merra ; alces.algD.merra.bm[8,2] <- alces.brtxsvm.bm.merra ; alces.algD.merra.bm[3,3] <- 1 ; alces.algD.merra.bm[4,3] <- alces.glmxglmnet.bm.merra ; alces.algD.merra.bm[5,3] <- alces.glmxmars.bm.merra ; alces.algD.merra.bm[6,3] <- alces.glmxmaxent.bm.merra ; alces.algD.merra.bm[7,3] <- alces.glmnetxmaxlike.bm.merra ; alces.algD.merra.bm[8,3] <- alces.glmnetxsvm.bm.merra ; alces.algD.merra.bm[4,4] <- 1 ; alces.algD.merra.bm[5,4] <- alces.glmnetxmars.bm.merra ; alces.algD.merra.bm[6,4] <- alces.glmnetxmaxent.bm.merra ; alces.algD.merra.bm[7,4] <- alces.glmnetxmaxlike.bm.merra ; alces.algD.merra.bm[8,4] <- alces.glmnetxsvm.bm.merra ; alces.algD.merra.bm[5,5] <- 1 ; alces.algD.merra.bm[6,5] <- alces.marsxmaxent.bm.merra ; alces.algD.merra.bm[7,5] <- alces.marsxmaxlike.bm.merra ; alces.algD.merra.bm[8,5] <- alces.marsxsvm.bm.merra ; alces.algD.merra.bm[6,6] <- 1 ; alces.algD.merra.bm[7,6] <- alces.maxentxmaxlike.bm.merra ; alces.algD.merra.bm[8,6] <- alces.maxentxsvm.bm.merra ; alces.algD.merra.bm[7,7] <- 1 ; alces.algD.merra.bm[8,7] <- alces.maxlikexsvm.bm.merra ; alces.algD.merra.bm[8,8] <- 1 ;
### Further tinkering to use corrplot package ####
alces.algD.bio.sm[upper.tri(alces.algD.bio.sm)] <- t(alces.algD.bio.sm)[upper.tri(alces.algD.bio.sm)]
alces.algD.bio.mm[upper.tri(alces.algD.bio.mm)] <- t(alces.algD.bio.mm)[upper.tri(alces.algD.bio.mm)]
alces.algD.bio.bm[upper.tri(alces.algD.bio.bm)] <- t(alces.algD.bio.bm)[upper.tri(alces.algD.bio.bm)]
alces.algD.envirem.sm[upper.tri(alces.algD.envirem.sm)] <- t(alces.algD.envirem.sm)[upper.tri(alces.algD.envirem.sm)]
alces.algD.envirem.mm[upper.tri(alces.algD.envirem.mm)] <- t(alces.algD.envirem.mm)[upper.tri(alces.algD.envirem.mm)]
alces.algD.envirem.bm[upper.tri(alces.algD.envirem.bm)] <- t(alces.algD.envirem.bm)[upper.tri(alces.algD.envirem.bm)]
alces.algD.merra.sm[upper.tri(alces.algD.merra.sm)] <- t(alces.algD.merra.sm)[upper.tri(alces.algD.merra.sm)]
alces.algD.merra.mm[upper.tri(alces.algD.merra.mm)] <- t(alces.algD.merra.mm)[upper.tri(alces.algD.merra.mm)]
alces.algD.merra.bm[upper.tri(alces.algD.merra.bm)] <- t(alces.algD.merra.bm)[upper.tri(alces.algD.merra.bm)]

### Plotting ####
setwd('..') ; dir.create('D-Stat-Algs') ; setwd('D-Stat-Algs')


corrplot.mixed(alces.algD.bio.sm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - BioClim - SM",mar=c(0,0,1,0))
dev.print(png,'alces.bio.sm.png',width=720,height=720)
corrplot.mixed(alces.algD.bio.mm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - BioClim - MM",mar=c(0,0,1,0))
dev.print(png,'alces.bio.mm.png',width=720,height=720)
corrplot.mixed(alces.algD.bio.bm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - BioClim - BM",mar=c(0,0,1,0))
dev.print(png,'alces.bio.bm.png',width=720,height=720)
corrplot.mixed(alces.algD.envirem.sm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - EnviREM - SM",mar=c(0,0,1,0))
dev.print(png,'alces.envirem.sm.png',width=720,height=720)
corrplot.mixed(alces.algD.envirem.mm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - EnviREM - MM",mar=c(0,0,1,0))
dev.print(png,'alces.envirem.mm.png',width=720,height=720)
corrplot.mixed(alces.algD.envirem.bm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - EnviREM - BM",mar=c(0,0,1,0))
dev.print(png,'alces.envirem.bm.png',width=720,height=720)
corrplot.mixed(alces.algD.merra.sm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - MERRAClim - SM",mar=c(0,0,1,0))
dev.print(png,'alces.merra.sm.png',width=720,height=720)
corrplot.mixed(alces.algD.merra.mm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - MERRAClim - MM",mar=c(0,0,1,0))
dev.print(png,'alces.merra.mm.png',width=720,height=720)
corrplot.mixed(alces.algD.merra.bm, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - MERRAClim - BM",mar=c(0,0,1,0))
dev.print(png,'alces.merra.bm.png',width=720,height=720)

### Averaging alg D-values for all species and Ms, to compare among 
dev.new()

alces.algD.avg <- (alces.algD.bio.bm+alces.algD.bio.mm+alces.algD.bio.sm+alces.algD.envirem.bm+alces.algD.envirem.mm+alces.algD.envirem.sm+alces.algD.merra.bm+alces.algD.merra.mm+alces.algD.merra.sm)/9

corrplot.mixed(alces.algD.avg, order = 'original',col.lim=c(0,1),is.corr = F,title="H.alces - Average",mar=c(0,0,1,0))
dev.print(png,'alces.average.png',width=720,height=720)

# # I also did this averaging for all species at once, with the lines below:
# algD.avg <- (alces.algD.avg+armiger.algD.avg+batesii.algD.avg+boterorum.algD.avg+cervinus.algD.avg+cheiracanthus.algD.avg+elaphus.algD.avg+longicornis.algD.avg+vesanicus.algD.avg)/9
# 
# corrplot.mixed(alces.algD.avg, order = 'original',col.lim=c(0,1),is.corr = F,title="D-Stat Average",mar=c(0,0,1,0))
# dev.print(png,'D.stat.Algs-average.png',width=720,height=720)

setwd('..')
save.image('script.RData')
