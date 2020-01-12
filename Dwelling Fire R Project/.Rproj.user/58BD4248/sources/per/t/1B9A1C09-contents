library(sf)
library(tmap)
library(tmaptools)
library(maptools)
library(plyr)
library(tidyverse)
library(rgdal)
library(spatstat)
library(rgeos)
library(raster)
library(adehabitatHR)
library(spdep)
library(classInt)
library(gstat)
library(spgwr)
library(PBSmapping)
tmap_mode("plot")

# road fire data & base map
dwellingFire <- read_csv("2018 dwelling fire.csv")
wardLondon <- st_read("statistical-gis-boundaries-london/ESRI/London_Ward.shp")
wardLondonOGR <- readOGR("statistical-gis-boundaries-london/ESRI/London_Ward.shp")
# set code for CRS
BNG = "+init=epsg:27700"

# reproject
dwellingFiresf <- st_as_sf(dwellingFire, coords = c("Easting_rounded", "Northing_rounded"), 
                   crs = 27700)
wardLondon <- st_transform(wardLondon, 27700)
wardLondonOGR <- spTransform(wardLondonOGR,BNG)

# define a function for spatial join
spatialJoin <- function(data1, data2) {
  # join OSM and London boroughs
  joined <- st_join(data1, data2, join = st_within)
  
  countno <- as.data.frame(plyr::count(joined$GSS_CODE))
  
  counted <-left_join(data2, countno, by=c("GSS_CODE"="x"))
  
  return(counted)
}

# assign dwelling fire incidents to each ward and rename the column
wardFire <- spatialJoin(dwellingFiresf, wardLondon)
names(wardFire)[names(wardFire) == 'freq'] = "Dwelling Fires"

# point map of dwelling fires
tm_shape(wardLondonOGR) +
  tm_polygons(col = NA, alpha = 0.5) +
  tm_shape(dwellingFiresf) +
  tm_dots(col = "blue")+ 
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Dwelling Fires in London", legend.position = c("left", "bottom"))

# classification plot of dwelling fires in each ward
tm_shape(wardFire) + 
  tm_polygons("Dwelling Fires", 
              style="jenks",
              palette=get_brewer_pal("Reds", n = 5, contrast = c(0, 0.6)),
              midpoint=NA,
              title="Dwelling Fires in Each Ward",
              alpha = 1) + 
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(legend.position = c("left", "bottom"))





## KDE analysis
# the first way
firePoint<-ggplot(dwellingFire, aes(x=Easting_rounded,y=Northing_rounded))+geom_point()+coord_equal()
firePoint
firePoint+stat_density2d(aes(fill = ..level..), geom="polygon")

# the second way
window <- as.owin(wardLondonOGR)
plot(window)
dwellingFire.ppp <- ppp(x=dwellingFire$Easting_rounded,y=dwellingFire$Northing_rounded,window=window)
plot(dwellingFire.ppp,pch=16,cex=0.4, main="Dwelling Fires in London")
plot(density(dwellingFire.ppp, sigma = 1000),main="KDE Plot of Dwelling Fires in London", palette="Blues")

# the third way
# transform the object to SpatialDataFrame
fireOGR <- as(dwellingFiresf, 'Spatial')
# calculate KDE result
kde.output <- kernelUD(fireOGR, h="href", grid = 1000)
kde <- raster(kde.output)
# reproject
projection(kde) <- CRS("+init=EPSG:27700")
bounding_box <- wardLondonOGR@bbox
# create masked map with the scale of London
masked_kde <- mask(kde, wardLondonOGR)

# KDE plot
tm_shape(masked_kde, bbox = bounding_box) + tm_raster("ud", style = "quantile", n = 100, legend.show = FALSE, palette = "OrRd") +
  tm_shape(wardLondonOGR) + tm_borders(alpha=0.3, col = "white") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(#title = "KDE Plot of Dwelling Fires in London",
            #title.position = c(0.1,0.95), 
            legend.position = c("left", "bottom"), frame = FALSE)





## Moran's I
# transform the object to SpatialDataFrame
wardFireOGR <- as(wardFire, 'Spatial')
# fill missing values with 0
wardFireOGR@data$Dwelling.Fires <- wardFireOGR@data$Dwelling.Fires %>% replace_na(0)

#calculate the centroids of all Wards in London
coordsW <- coordinates(wardFireOGR)
#create a neighbours list
FWard_nb <- poly2nb(wardFireOGR, queen=T)
#create a spatial weights object from these weights
Fward.lw <- nb2listw(FWard_nb, style="C")

# global Moran test
I_FWard_Global <- moran.test(wardFireOGR@data$Dwelling.Fires, Fward.lw)
I_FWard_Global
# local Moran's I
I_FWard_Local <- localmoran(wardFireOGR@data$Dwelling.Fires, Fward.lw)

# store local Moran's I & Z-value
wardFireOGR@data$FireI <- I_FWard_Local[,1]
wardFireOGR@data$FireIz <- I_FWard_Local[,4]

# local Moran's I plot
fireIMap <- tm_shape(wardFireOGR) +
  tm_polygons("FireI",
              style="fixed",
              breaks=c(-Inf,-1,-0.5,0,0.5,1,Inf),
              palette=get_brewer_pal("-RdBu", contrast = c(0, 0.6)),
              alpha = .7,
              midpoint=NA,
              legend.format = list(digits=2),
              title="Local Moran's I of Dwelling Fires") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Local Moran's I of Dwelling Fires",
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"), frame = TRUE)

fireIzMap <- tm_shape(wardFireOGR) +
  tm_polygons("FireIz",
              style="fixed",
              breaks=c(-1000,-2.58,-1.96,1.96,2.58,1000),
              palette=get_brewer_pal("-PuOr", contrast = c(0, 0.6)),
              alpha = .7,
              midpoint=NA,
              legend.format = list(digits=2),
              title="Z-value of Local Moran's I \nof Dwelling Fires") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Sig Map of Local Moran's I of Dwelling Fires",
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"), frame = TRUE)
tmap_arrange(fireIMap, fireIzMap)

# significant local Moran's I map of dwelling fires
tm_shape(wardFireOGR) +
  tm_polygons(alpha = 0)+
  tm_shape(wardFireOGR[(wardFireOGR@data$FireIz >= 1.96)|(wardFireOGR@data$FireIz <= -1.96), ]) +
  tm_polygons("FireI",
              style="fixed",
              breaks=c(-Inf,-1,-0.5,0,0.5,1,Inf),
              palette=get_brewer_pal("-RdBu", n = 5, contrast = c(0, 0.6)),
              alpha = 1,
              midpoint=NA,
              legend.format = list(digits=2),
              title="Local Moran's I \nof Dwelling Fires") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(#title = "Local Moran's I of Dwelling Fires in Each Ward of London",
            #title.position = c(0.1,0.95),
            legend.position = c("left", "bottom"), frame = FALSE)





## LISA analysis of dwelling fires
# standarise dwelling fire variable
wardFireOGR@data$sWardFire <- scale(wardFireOGR@data$Dwelling.Fires)
# create a lagged variable
wardFireOGR@data$lag_sWardFire <- lag.listw(Fward.lw, wardFireOGR@data$sWardFire)

# scatter plot for general understanding
plot(x = wardFireOGR@data$sWardFire, y = wardFireOGR@data$lag_sWardFire, main = " Moran Scatterplot")
abline(h = 0, v = 0)
abline(lm(wardFireOGR@data$lag_sWardFire ~ wardFireOGR@data$sWardFire), lty = 3, lwd = 4, col = "red")

# identify the moran plot quadrant for each observation
wardFireOGR@data$quad_sig <- NA
wardFireOGR@data[(wardFireOGR@data$sWardFire >= 0 & wardFireOGR@data$lag_sWardFire >= 0) & (I_FWard_Local[, 5] <= 0.05), "quad_sig"] <- 1
wardFireOGR@data[(wardFireOGR@data$sWardFire <= 0 & wardFireOGR@data$lag_sWardFire <= 0) & (I_FWard_Local[, 5] <= 0.05), "quad_sig"] <- 2
wardFireOGR@data[(wardFireOGR@data$sWardFire >= 0 & wardFireOGR@data$lag_sWardFire <= 0) & (I_FWard_Local[, 5] <= 0.05), "quad_sig"] <- 3
wardFireOGR@data[(wardFireOGR@data$sWardFire <= 0 & wardFireOGR@data$lag_sWardFire >= 0) & (I_FWard_Local[, 5] <= 0.05), "quad_sig"] <- 4
wardFireOGR@data[(I_FWard_Local[, 5] > 0.05), "quad_sig"] <- 5  # All insignificant observations

# LISA plot
breaks <- seq(1, 5, 1)
labels <- c("high-High", "low-Low", "High-Low", "Low-High", "Not Significant")
np <- findInterval(wardFireOGR@data$quad_sig, breaks)
colors <- c("red", "blue", "lightpink", "skyblue2", "white")
plot(wardFireOGR, col = colors[np])
mtext("LISA Map of Dwelling Fires in Each Ward of London", cex = 1.5, side = 3, line = 1)
legend("bottomleft", legend = labels, fill = colors, bty = "n")





## IMD analysis
# read in 2019 London IMD data
imd2019 = read_csv("2019 Ward Average IMD Score.csv")
# join it to ward map
IMDmap<-merge(wardLondonOGR, 
                      imd2019, 
                      by.x="GSS_CODE", 
                      by.y="Ward Code",
                      no.dups = TRUE)

# IMD Moran's I Analysis

# calculate centroid of each ward
coordsIMD <- coordinates(IMDmap)
#create a neighbours list
IMD_nb <- poly2nb(IMDmap, queen=T)
#create a spatial weights object from these weights
IMD.lw <- nb2listw(IMD_nb, style="C")

# global Moran test
I_IMD_Global <- moran.test(IMDmap@data$`2019 Average IMD Score`, IMD.lw)
I_IMD_Global

# local Moran's I test
I_IMD_Local <- localmoran(IMDmap@data$`2019 Average IMD Score`, IMD.lw)
# store local Moran's I & Z-value
IMDmap@data$imdI <- I_IMD_Local[,1]
IMDmap@data$imdIz <- I_IMD_Local[,4]

imdIMap <- tm_shape(IMDmap) +
  tm_polygons("imdI",
              style="fixed",
              breaks=c(-Inf,-1,-0.5,0,0.5,1,Inf),
              palette=get_brewer_pal("-RdBu", contrast = c(0, 0.6)),
              alpha = .7,
              midpoint=NA,
              legend.format = list(digits=2),
              title="Local Moran's I of IMD Scores") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Local Moran's I of IMD Scores",
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"), frame = TRUE)

imdIzMap <- tm_shape(IMDmap) +
  tm_polygons("imdIz",
              style="fixed",
              breaks=c(-1000,-2.58,-1.96,1.96,2.58,1000),
              palette=get_brewer_pal("-PuOr", contrast = c(0, 0.6)),
              alpha = .7,
              midpoint=NA,
              legend.format = list(digits=2),
              title="Z-value of Local Moran's I \nof IMD Scores") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Sig Map of Local Moran's I of IMD Scores",
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"), frame = TRUE)
tmap_arrange(imdIMap,imdIzMap)


# significant local Moran plot
tm_shape(IMDmap) +
  tm_polygons(alpha = 0)+
  tm_shape(IMDmap[(IMDmap@data$imdIz >= 1.96)|(IMDmap@data$imdIz <= -1.96), ]) +
  tm_polygons("imdI",
              style="fixed",
              breaks=c(-Inf,-1,-0.5,0,0.5,1,Inf),
              palette=get_brewer_pal("-RdBu", n = 5, contrast = c(0, 0.6)),
              midpoint=NA,
              legend.format = list(digits=2),
              title="Local Moran's I \nof IMD Score") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(#title = "Local Moran's I of IMD Score in Each Ward of London",
            #title.position = c(0.1,0.95),
            legend.position = c("left", "bottom"), frame = FALSE)


# Overlay Moran Plot of Fires & IMD Score
tm_shape(wardFireOGR) +
  tm_polygons(alpha = 0)+
  tm_shape(wardFireOGR[(wardFireOGR@data$FireIz >= 1.96)|(wardFireOGR@data$FireIz <= -1.96), ]) +
  tm_polygons("FireI",
              style="fixed",
              breaks=c(-Inf,-1,-0.5,0,0.5,1,Inf),
              palette=get_brewer_pal("-RdBu", contrast = c(0, 0.6)),
              midpoint=NA,
              alpha=1,
              legend.show = FALSE) +
  tm_shape(IMDmap[(IMDmap@data$imdIz >= 1.96)|(IMDmap@data$imdIz <= -1.96), ]) +
  tm_polygons("imdI",
              style="fixed",
              breaks=c(-Inf,-1,-0.5,0,0.5,1,Inf),
              palette=get_brewer_pal("-RdBu", contrast = c(0, 0.6)),
              midpoint=NA,
              alpha = 0.5,
              legend.format = list(digits=2),
              title="Local Moran's I of \nDwelling Fires & IMD Scores") +
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(legend.position = c("left", "bottom"), frame = FALSE)





## Spatially Weighted Regression (GWR)
# read in dataset of explicative variables
factorsGWR = read_csv("Fire GWR Factors.csv")
# join it to ward fire dataset
wardFireOGR<-merge(wardFireOGR, 
              factorsGWR, 
              by.x="GSS_CODE", 
              by.y="Ward Code",
              no.dups = TRUE)

# join IMD data into ward fire dataset
wardFireOGR<-merge(wardFireOGR, 
                   IMDmap@data[c("2019 Average IMD Score","GSS_CODE")],
                   by.x="GSS_CODE", 
                   by.y="GSS_CODE",
                   no.dups = TRUE)

# drop all rows with missing values
wardFireOGR <- wardFireOGR[!is.na(wardFireOGR@data$`Ward Name`) ,]

# create Pearson correlation matrix
variables <- c('Dwelling.Fires','Population Density (km2)','2019 Average IMD Score',
               'Median_Annual_Household_Income','Employment rate (16-64) - 2011',
               'Percentage of Buildings','Percentage of Greenspace')
corrMatrix <- cor(wardFireOGR@data[variables],use='complete.obs',method = 'pearson')

## run GWR analysis
coordsW <- coordinates(wardFireOGR)
GWRbandwidth <- gwr.sel(`Dwelling.Fires` ~ `Population Density (km2)` + `2019 Average IMD Score` + `Percentage of Buildings` + `Employment rate (16-64) - 2011`, data = wardFireOGR@data, coords=coordsW,adapt=T)
gwr.model = gwr(`Dwelling.Fires` ~ `Population Density (km2)` + `2019 Average IMD Score` + `Percentage of Buildings` + `Employment rate (16-64) - 2011`, data = wardFireOGR@data, coords=coordsW, adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE)
gwr.model
# store GWR results
results<-as.data.frame(gwr.model$SDF)

#run the significance test
popDensity_sigTest = abs(results$X.Population.Density..km2..) -1.96 * results$X.Population.Density..km2.._se
building_sigTest = abs(results$X.Percentage.of.Buildings.) -1.96 * results$X.Percentage.of.Buildings._se
employment_sigTest = abs(results$X.Employment.rate..16.64....2011.) -1.96 * results$X.Employment.rate..16.64....2011._se
IMD_sigTest = abs(results$X.2019.Average.IMD.Score.) -1.96 * results$X.2019.Average.IMD.Score._se

#store coefficient results
wardFireOGR@data$'coef.Population Density' <- results$X.Population.Density..km2..
wardFireOGR@data$'coef.Percentage of Buildings' <- results$X.Percentage.of.Buildings.
wardFireOGR@data$'coef.Employment Rate' <- results$X.Employment.rate..16.64....2011.
wardFireOGR@data$'coef.IMD Scores' <- results$X.2019.Average.IMD.Score.
#store significance results
wardFireOGR@data$'sig of Population Density'<-popDensity_sigTest
wardFireOGR@data$'sig of Building Percentage'<-building_sigTest
wardFireOGR@data$'sig of Employment Rate'<-employment_sigTest
wardFireOGR@data$'sig of IMD Scores'<-IMD_sigTest


# coefficient plots through Wards
coefPopMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'coef.Population Density', palette = "-Spectral", alpha = 0.5)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Coefficient Map of Population Density", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))

coefIMDMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'coef.IMD Scores', palette = "-Spectral", alpha = 0.5)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Coefficient Map of IMD Scores", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))

coefEmploymentMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'coef.Employment Rate', palette = "-Spectral", alpha = 0.5)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Coefficient Map of Employment Rate", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))

coefBuildingMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'coef.Percentage of Buildings', palette = "-Spectral", alpha = 0.5)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Coefficient Map of Building Percentage", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))
tmap_arrange(coefPopMap, coefIMDMap, coefEmploymentMap, coefBuildingMap)


# significance plots
sigPopMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'sig of Population Density', palette = "-RdGy")+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Significance Map of Population Density", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))

sigIMDMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'sig of IMD Scores', palette = "-RdGy")+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Significance Map of Ward IMD Scores", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))

sigEmploymentMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'sig of Employment Rate', palette = "-RdGy")+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Significance Map of Employment Rate", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))

sigBuildingMap <- tm_shape(wardFireOGR) +
  tm_polygons(col = 'sig of Building Percentage', palette = "-RdGy")+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(title = "Significance Map of Building Percentage", 
            title.position = c(0.01,0.95),
            legend.position = c("left", "bottom"))
tmap_arrange(sigPopMap, sigIMDMap, sigEmploymentMap, sigBuildingMap)


# significant coefficient maps of 4 factors
tm_shape(wardFireOGR) +
  tm_polygons(alpha = 0)+
  tm_shape(wardFireOGR[(wardFireOGR@data$`sig of Population Density` >= 0), ]) +
  tm_polygons(col = 'coef.Population Density', palette = "-RdBu", alpha = 0.7)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(legend.position = c("left", "bottom"), frame = FALSE)

tm_shape(wardFireOGR) +
  tm_polygons(alpha = 0)+
  tm_shape(wardFireOGR[(wardFireOGR@data$`sig of IMD Scores` >= 0), ]) +
  tm_polygons(col = 'coef.IMD Scores', palette = "Oranges", alpha = 0.7)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(legend.position = c("left", "bottom"), frame = FALSE)

tm_shape(wardFireOGR) +
  tm_polygons(alpha = 0)+
  tm_shape(wardFireOGR[(wardFireOGR@data$`sig of Employment Rate` >= 0), ]) +
  tm_polygons(col = 'coef.Employment Rate', palette = get_brewer_pal("-PuOr", n = 4, contrast = c(0, 0.4)), alpha = 0.7)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(legend.position = c("left", "bottom"), frame = FALSE)

tm_shape(wardFireOGR) +
  tm_polygons(alpha = 0)+
  tm_shape(wardFireOGR[(wardFireOGR@data$`sig of Building Percentage`>= 0), ]) +
  tm_polygons(col = 'coef.Percentage of Buildings', palette = "PRGn", alpha = 0.7)+
  tm_compass(position = c("right", "top"),type = "arrow") + 
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(legend.position = c("left", "bottom"), frame = FALSE)

