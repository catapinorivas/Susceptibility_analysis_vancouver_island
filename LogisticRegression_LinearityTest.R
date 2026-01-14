
# Package Installation ----------------------------------------------------


install.packages("caTools")    # For Logistic regression
install.packages("ROCR")
install.packages("sp")
install.packages("raster")
install.packages("corrplot")
install.packages("pROC")
install.packages("car")
install.packages("caret")
install.packages("rms")
install.packages("ResourceSelection")
install.packages("PredictABEL")
install.packages("rasterVis")
install.packages("glmtoolbox")
install.packages('pscl')
install.packages("sf")
install.packages("glmnet")
install.packages("fastDummies")
install.packages("Gifi")
install.packages("MLmetrics")
install.packages("remotes")
# remotes::install_github("rvalavi/myspatial")
# remotes::install_github("hongooi73/glmnetUtils")
install.packages("tibble")
install.packages("lessR")
install.packages("svglite")
install.packages("ggplot2")


## Install required modules
library("dplyr")
library(caTools)
library(pscl)
library(ROCR) 
library(sp)
library(raster)
library(rgdal)
library(pROC)
library(car)
library(caret)
library(rms)
library(ResourceSelection)
library(PredictABEL)
library(rasterVis)
library(glmtoolbox)
library(sf)
library(glmnet)
library(fastDummies)
library(Gifi)
library(ModelMetrics)
library(remotes)
# library(myspatial)
# library(glmnetUtils)
library(tibble)
library(lessR)
library(svglite)
library(ggplot2)

# Thematic variable import (12.5 m) ------------------------------------------

# Thematic variable import (12.5 m) ------------------------------------------
elevation <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\elevation.tif")
ndvi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\NDVI.tif")
plCurv <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\plancurv.tif")
prCurv <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\prcurve.tif")
roadsDist <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\dist_roads.tif")
slope <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\slope.tif")
spi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\SPI.tif")
tpi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\TPIs\\TPI_clasif6cell.tif")
twi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\TWI.tif")
waterDist <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\dist_rivers.tif")
vegage <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\veg_age.tif')
faultDist <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\faults.tif')
map <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\map.tif")


#psa_0p3 <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression Jan 2025\\PSA_0p3.tif')
#psa_1p0 <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression Jan 2025\\PSA_1p0.tif')
#psa_3p0 <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression Jan 2025\\PSA_3p0.tif')
#relief_11 <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression Jan 2025\\RELIEF_11.tif')
#relief_21 <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression Jan 2025\\RELIEF_21.tif')
#relief_31 <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression Jan 2025\\RELIEF_31.tif')
#mat <- raster("C:\\Users\\jpc19\\Documents\\3_MASKED\\3_MASKED\\mat_7")

# Import categorical data
geology <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\geology.tif")
landCover <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\landcover.tif")
aspect <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\aspect.tif")
aspectCriteria <- matrix(c(-5, -1, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, -1, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 380, 9
                           , 1, 2, 3, 4, 5, 6, 7, 8, 1), ncol = 3)
aspect <- reclassify(aspect, aspectCriteria)

# Import and Split Landslide Database -------------------------------------

# Open polygons from shapefile

ls_Polygons <- sf::st_read("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\Inventory_sources\\sources.shp")
lsPolygons_sf <- ls_Polygons$geometry

# Remove z component of sf
allPolygons <- sf::st_zm(lsPolygons_sf)

# Convert sf to spatial object
allPolygons_spatial <- sf::as_Spatial(allPolygons, cast = TRUE)

# Rasterize polygons, with all polygons appearing as 1
allLS <- rasterize(allPolygons_spatial, elevation, field = 1)

# Determine number of landslide pixels
lsCount <- length(lsPolygons_sf)

# Extract downsampled non-landslide points
ls_full <- reclassify(allLS, cbind(NA, 0)) 
ls_full <- raster::mask(x = ls_full, mask = elevation)
allNLS <- reclassify(ls_full, cbind(1, NA))

# Downsample points at 5:1 ratio
NLS_DS <- sampleRandom(allNLS, lsCount * 5, asRaster = TRUE)

# Generate predictor variable Raster Stacks --------------------------------------------------

# Create an initial stack of all pixels in study area
# mat and tri removed due to multicollinearity

predictStackNum <- stack(aspect, vegage, ndvi, roadsDist, waterDist, faultDist, map, tpi, prCurv, plCurv, slope, twi, geology, landCover, elevation, pga, spi)
names(predictStackNum) <- c("aspect", "vegage", "ndvi", "roadsDist", "waterDist", "faultDist", "map", "tpi", "prCurv", "plCurv", "slope", "twi", "geology", "landCover", "elevation", "pga", "spi")

# Average values over each polygon

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

MyMedian <- function(x){
  outVar <- median(x, na.rm = T)
}

StackLS.df <- data.frame(matrix(nrow = length(allPolygons_spatial), ncol = 0))
for (jj in 1:nlayers(predictStackNum)){
  tmpRas <- raster::extract(x = predictStackNum[[jj]], y = allPolygons_spatial)
  tmpRas.list <- lapply(tmpRas, MyMedian)
  tmpRas.df <- t(data.frame(tmpRas.list))
  StackLS.df <- cbind(StackLS.df, tmpRas.df)
}

colnames(StackLS.df) <- names(predictStackNum)

# Seperate raster stacks into LS train/test and NLS
stackNLS <- raster::mask(predictStackNum, NLS_DS)


# Convert seperated stacks into data frames
stackNLS.df <- na.omit(as.data.frame(stackNLS))

data.df <- rbind(StackLS.df, stackNLS.df)



# Create landslides vector ---------------------------

landslides.vec <- c(rep(1, length(StackLS.df$elevation)), rep(0, length(stackNLS.df$elevation)))

lsFull_masked <- raster::mask(ls_full, predictStackNum$elevation)
predictLS.df <- as.data.frame(landslides.vec)


# Assess linearity in continuous predictors -------------------------------

# Remove non-continuous variables


tidwell.df <- data.df

# Transform variables into all positive values

for (ii in 1:ncol(tidwell.df)){
  if (min(tidwell.df[,ii]) <= 0){
    tidwell.df[,ii] <- tidwell.df[,ii] + abs(min(tidwell.df[,ii])) + 1
  }

}


tidwell.df <- cbind(predictLS.df, tidwell.df)
colnames(tidwell.df)[1] <- "landslides"




# Use Box-Tidwell test to assess linearity between logit of outcome and predictors
# Check p-value of interaction (predictor_variable:log(predictor_variable))
# p-values >= 0.05 mean no evidence of non-linear relationship, so can stay continuous
# Otherwise, check transformations
elev1 <- lessR::Logit(landslides ~ elevation + elevation:log(elevation), data = tidwell.df, brief = TRUE) #transform?
ndvi1 <- lessR::Logit(landslides ~ ndvi + ndvi:log(ndvi), data = tidwell.df, brief = TRUE) #transform? 
prCurv1 <- lessR::Logit(landslides ~ prCurv + prCurv:log(prCurv), data = tidwell.df, brief = TRUE) #no transform, but transformed anyway
plCurv1 <- lessR::Logit(landslides ~ plCurv + plCurv:log(plCurv), data = tidwell.df, brief = TRUE) #it seems like it needs to be transformed
roadsDist1 <- lessR::Logit(landslides ~ roadsDist + roadsDist:log(roadsDist), data = tidwell.df, brief = TRUE) #it seems like it needs to be transformed
slope1 <- lessR::Logit(landslides ~ slope + slope:log(slope), data = tidwell.df, brief = TRUE) #it needs transformation
spi1 <- lessR::Logit(landslides ~ spi + spi:log(spi), data = tidwell.df, brief = TRUE) #not enough evidence to suggest transformation
tpi1 <- lessR::Logit(landslides ~ tpi + tpi:log(tpi), data = tidwell.df, brief = TRUE) #transformed
map1 <- lessR::Logit(landslides ~ map + map:log(map), data = tidwell.df, brief = TRUE) #transformed
twi1 <- lessR::Logit(landslides ~ twi + twi:log(twi), data = tidwell.df, brief = TRUE) #it needs to be transformed
waterDist1 <- lessR::Logit(landslides ~ waterDist + waterDist:log(waterDist), data = tidwell.df, brief = TRUE) #no transformation needed
pga1 <- lessR::Logit(landslides ~ pga + pga:log(pga), data = tidwell.df, brief = TRUE) #no transformation
#pgv1 <- lessR::Logit(landslides ~ pgv + pgv:log(pgv), data = tidwell.df, brief = TRUE) #no transformation
faultDist1 <- lessR::Logit(landslides ~ faultDist + faultDist:log(faultDist), data = tidwell.df, brief = TRUE) #transformed
#psa_0p31 <- lessR::Logit(landslides ~ psa_0p3 + psa_0p3:log(psa_0p3), data = tidwell.df, brief = TRUE) #no transformation
#psa_1p01 <- lessR::Logit(landslides ~ psa_1p0 + psa_1p0:log(psa_1p0), data = tidwell.df, brief = TRUE) #no transformation
#psa_3p01 <- lessR::Logit(landslides ~ psa_3p0 + psa_3p0:log(psa_3p0), data = tidwell.df, brief = TRUE) #no transformation
#relief_111 <- lessR::Logit(landslides ~ relief_11 + relief_11:log(relief_11), data = tidwell.df, brief = TRUE) #no transformation
#relief_211 <- lessR::Logit(landslides ~ relief_21 + relief_21:log(relief_21), data = tidwell.df, brief = TRUE) #no transformation
#relief_311 <- lessR::Logit(landslides ~ relief_31 + relief_31:log(relief_31), data = tidwell.df, brief = TRUE) #no transformation



### SUMMARY
# elevation check transformations
# ndvi check transformations 
# prCurv check transformations
# plCurv check transformations
# roadsDist ok (discretize due to lack of meaning though)
# slope check transformations
# twi check transformations


# Plot relationships between continuous vars and log odds -------------------------------------------

# Elevation
elevProbs.vals <- predict(elev1, type = "response")
elevProbs.log <- log(elevProbs.vals / (1 - elevProbs.vals))
elevLR <- lm(tidwell.df$elevation ~ elevProbs.log)
ggplot(data = data.frame(x = tidwell.df$elevation, y = elevProbs.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\elevation.svg")
summary(elevLR)


# ndvi
ndviProbs.vals <- predict(ndvi1, type = "response")
ndviProbs.log <- log(ndviProbs.vals / (1 - ndviProbs.vals))
ndviLR <- lm(tidwell.df$ndvi ~ ndviProbs.log)
plot(tidwell.df$ndvi, ndviProbs.log, cex = 0.25, col = "black")
ggplot(data = data.frame(x = tidwell.df$ndvi, y = ndviProbs.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\ndvi.svg")
summary(ndviLR)

# prCurv
prCurvProbs.vals <- predict(prCurv1, type = "response")
prCurvProbs.log <- log(prCurvProbs.vals / (1 - prCurvProbs.vals))
prCurvLR <- lm(tidwell.df$prCurv ~ prCurvProbs.log)
plot(tidwell.df$prCurv, prCurvProbs.log, cex = 0.25, col = "black")
ggplot(data = data.frame(x = tidwell.df$prCurv, y = prCurvProbs.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\prCurv.svg")
summary(prCurvLR)

# plCurv
plCurvProbs.vals <- predict(plCurv1, type = "response")
plCurvProbs.log <- log(plCurvProbs.vals / (1 - plCurvProbs.vals))
plCurvLR <- lm(tidwell.df$plCurv ~ plCurvProbs.log)
plot(tidwell.df$plCurv, plCurvProbs.log, cex = 0.25, col = "black")
ggplot(data = data.frame(x = tidwell.df$plCurv, y = plCurvProbs.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\plCurv.svg")
summary(plCurvLR)


# roadsDist
roadsDistProbs.vals <- predict(roadsDist1, type = "response")
roadsDistProbs.log <- log(roadsDistProbs.vals / (1 - roadsDistProbs.vals))
roadsDistLR <- lm(tidwell.df$roadsDist ~ roadsDistProbs.log)
ggplot(data = data.frame(x = tidwell.df$roadsDist, y = roadsDistProbs.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\roadsDist.svg")
summary(roadsDistLR)

# slope
slopeProbs.vals <- predict(slope1, type = "response")
slopeProbs.log <- log(slopeProbs.vals / (1 - slopeProbs.vals))
slopeLR <- lm(tidwell.df$slope ~ slopeProbs.log)
ggplot(data = data.frame(x = tidwell.df$slope, y = slopeProbs.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\slope.svg")
# abline(v = 12, col = "red")
# abline(v = 25, col = "red")
# abline(v = 42, col = "red")
summary(slopeLR)


# twi
twiProbs.vals <- predict(twi1, type = "response")
twiProbs.log <- log(twiProbs.vals / (1 - twiProbs.vals))
twiLR <- lm(tidwell.df$twi ~ twiProbs.log)
ggplot(data = data.frame(x = tidwell.df$twi, y = twiProbs.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\twi.svg")
summary(twiLR)


# Perform var transformations and re-plot ---------------------------------------------
squaredTidwell.df <- tidwell.df
squaredTidwell.df[,2:ncol(tidwell.df)] <- squaredTidwell.df[,2:ncol(tidwell.df)] ^ 2
sqrtTidwell.df <- tidwell.df
sqrtTidwell.df[,2:ncol(tidwell.df)] <- sqrt(tidwell.df[,2:ncol(tidwell.df)]) ### Data must be all positive, or apply shift
invTidwell.df <- tidwell.df
invTidwell.df[,2:ncol(tidwell.df)] <- 1 / tidwell.df[,2:ncol(tidwell.df)] ### Data must be all positive, or apply shift


# ndvi

ndviSqrt <- lessR::Logit(landslides ~ ndvi + ndvi:log(ndvi), data = sqrtTidwell.df, brief = TRUE)

ndviSqrt.vals <- predict(ndviSqrt, type = "response")
ndviSqrt.log <- log(ndviSqrt.vals / (1 - ndviSqrt.vals))
ndviSqrt.LR <- lm(sqrtTidwell.df$ndvi ~ ndviSqrt.log)
ggplot(data = data.frame(x = sqrtTidwell.df$ndvi, y = ndviSqrt.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                         
ggsave("LinearityTest_2\\ndvi_sqrt.svg")
summary(ndviSqrt.LR)

ndviInv <- lessR::Logit(landslides ~ ndvi + ndvi:log(ndvi), data = invTidwell.df, brief = TRUE)

ndviInv.vals <- predict(ndviInv, type = "response")
ndviInv.log <- log(ndviInv.vals / (1 - ndviInv.vals))
ndviInv.LR <- lm(invTidwell.df$ndvi ~ ndviInv.log)
ggplot(data = data.frame(x = invTidwell.df$ndvi, y = ndviInv.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\ndvi_inv.svg")
summary(ndviInv.LR)

ndviSquared <- lessR::Logit(landslides ~ ndvi + ndvi:log(ndvi), data = squaredTidwell.df, brief = TRUE)

ndviSquared.vals <- predict(ndviSquared, type = "response")
ndviSquared.log <- log(ndviSquared.vals / (1 - ndviSquared.vals))
ndviSquared.LR <- lm(squaredTidwell.df$ndvi ~ ndviSquared.log)
ggplot(data = data.frame(x = squaredTidwell.df$ndvi, y = ndviSquared.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")
ggsave("LinearityTest_2\\ndvi_squared.svg")
summary(ndviSquared.LR)


# roadsDist - I DECIDED TO DISCRETIZE

roadsDistSqrt <- lessR::Logit(landslides ~ roadsDist + roadsDist:log(roadsDist), data = sqrtTidwell.df, brief = TRUE)

roadsDistSqrt.vals <- predict(roadsDistSqrt, type = "response")
roadsDistSqrt.log <- log(roadsDistSqrt.vals / (1 - roadsDistSqrt.vals))
roadsDistSqrt.LR <- lm(sqrtTidwell.df$roadsDist ~ roadsDistSqrt.log)
ggplot(data = data.frame(x = sqrtTidwell.df$roadsDist, y = roadsDistSqrt.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                             
ggsave("LinearityTest_2\\roadsDist_sqrt.svg")
summary(roadsDistSqrt.LR)

roadsDistInv <- lessR::Logit(landslides ~ roadsDist + roadsDist:log(roadsDist), data = invTidwell.df, brief = TRUE)

roadsDistInv.vals <- predict(roadsDistInv, type = "response")
roadsDistInv.log <- log(roadsDistInv.vals / (1 - roadsDistInv.vals))
roadsDistInv.LR <- lm(invTidwell.df$roadsDist ~ roadsDistInv.log)
ggplot(data = data.frame(x = invTidwell.df$roadsDist, y = roadsDistInv.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                       
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\roadsDist_inv.svg")
summary(roadsDistInv.LR)

roadsDistSquared <- lessR::Logit(landslides ~ roadsDist + roadsDist:log(roadsDist), data = squaredTidwell.df, brief = TRUE)

roadsDistSquared.vals <- predict(roadsDistSquared, type = "response")
roadsDistSquared.log <- log(roadsDistSquared.vals / (1 - roadsDistSquared.vals))
roadsDistSquared.LR <- lm(squaredTidwell.df$roadsDist ~ roadsDistSquared.log)
ggplot(data = data.frame(x = squaredTidwell.df$roadsDist, y = roadsDistSquared.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                             
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\roadsDist_squared.svg")
summary(roadsDistSquared.LR)

# Road dist squared shows linear relationship


# slope

slopeSqrt <- lessR::Logit(landslides ~ slope + slope:log(slope), data = sqrtTidwell.df, brief = TRUE)

slopeSqrt.vals <- predict(slopeSqrt, type = "response")
slopeSqrt.log <- log(slopeSqrt.vals / (1 - slopeSqrt.vals))
slopeSqrt.LR <- lm(sqrtTidwell.df$slope ~ slopeSqrt.log)
ggplot(data = data.frame(x = sqrtTidwell.df$slope, y = slopeSqrt.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                        
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\slopet_sqrt.svg")
summary(slopeSqrt.LR)

slopeInv <- lessR::Logit(landslides ~ slope + slope:log(slope), data = invTidwell.df, brief = TRUE)

slopeInv.vals <- predict(slopeInv, type = "response")
slopeInv.log <- log(slopeInv.vals / (1 - slopeInv.vals))
slopeInv.LR <- lm(invTidwell.df$slope ~ slopeInv.log)
ggplot(data = data.frame(x = invTidwell.df$slope, y = slopeInv.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                                
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\slope_inv.svg")
summary(slopeInv.LR)

slopeSquared <- lessR::Logit(landslides ~ slope + slope:log(slope), data = squaredTidwell.df, brief = TRUE)

slopeSquared.vals <- predict(slopeSquared, type = "response")
slopeSquared.log <- log(slopeSquared.vals / (1 - slopeSquared.vals))
slopeSquared.LR <- lm(squaredTidwell.df$slope ~ slopeSquared.log)
ggplot(data = data.frame(x = squaredTidwell.df$slope, y = slopeSquared.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                              
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\slope_squared.svg")
summary(slopeSquared.LR)




# twi

twiSqrt <- lessR::Logit(landslides ~ twi + twi:log(twi), data = sqrtTidwell.df, brief = TRUE)

twiSqrt.vals <- predict(twiSqrt, type = "response")
twiSqrt.log <- log(twiSqrt.vals / (1 - twiSqrt.vals))
twiSqrt.LR <- lm(sqrtTidwell.df$twi ~ twiSqrt.log)
ggplot(data = data.frame(x = sqrtTidwell.df$twi, y = twiSqrt.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                      
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\twi_sqrt.svg")
summary(twiSqrt.LR)

twiInv <- lessR::Logit(landslides ~ twi + twi:log(twi), data = invTidwell.df, brief = TRUE)

twiInv.vals <- predict(twiInv, type = "response")
twiInv.log <- log(twiInv.vals / (1 - twiInv.vals))
twiInv.LR <- lm(invTidwell.df$twi ~ twiInv.log)
ggplot(data = data.frame(x = invTidwell.df$twi, y = twiInv.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                                       
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\twi_inv.svg")
summary(twiInv.LR)

twiSquared <- lessR::Logit(landslides ~ twi + twi:log(twi), data = squaredTidwell.df, brief = TRUE)

twiSquared.vals <- predict(twiSquared, type = "response")
twiSquared.log <- log(twiSquared.vals / (1 - twiSquared.vals))
twiSquared.LR <- lm(squaredTidwell.df$twi ~ twiSquared.log)
ggplot(data = data.frame(x = squaredTidwell.df$twi, y = twiSquared.log), aes(x = x, y = y)) + geom_point(size = 1, shape = 1) + 
  geom_smooth(method='lm', formula= y~x, linetype = "dashed", color = "red")                                                                         
ggsave("C:\\Users\\GT\\OneDrive - Simon Fraser University (1sfu)\\Documents\\RProjects\\20230124_Meager\\LinearityTest_2\\twi_squared.svg")
summary(twiSquared.LR)




# Resources ----------------------------------------------------------------

# https://janhove.github.io/analysis/2015/10/16/nonlinear-relationships
# https://r4ds.github.io/bookclub-islr/addendum---logistic-regression-assumptions.html
# https://www.youtube.com/watch?v=O7gRceyeyT8&t=1123s
# https://rforhr.com/logistic.html

# https://datascienceplus.com/fitting-polynomial-regression-r/
# https://openforecast.org/sba/types-of-variables-transformations.html