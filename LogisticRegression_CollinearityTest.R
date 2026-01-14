#PEARSON CORRELATION TEST
# Package Installation ----------------------------------------------------


install.packages("caTools")    # For Logistic regression
install.packages("ROCR")
install.packages("sp")
install.packages("raster")
# install.packages("arcgisbinding", repos="https://r.esri.com", type="win.binary")
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
install.packages("Hmisc")

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
# library(rvalavi/myspatial)
library(remotes)
# library(myspatial)
# library(glmnetUtils)
library(tibble)
library(Hmisc)



## Import continuous data

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

# Generate  Continuous variable Raster Stacks --------------------------------------------------

# Create an initial stack of all pixels in study area


predictStack <- stack(elevation, aspect, vegage, ndvi, roadsDist, waterDist, faultDist, map, tpi, prCurv, plCurv, slope, twi, geology, landCover, spi)
names(predictStack) <- c("elevation", "aspect", "vegage", "ndvi", "roadsDist", "waterDist", "faultDist", "map", "tpi", "prCurv", "plCurv", "slope", "twi", "geology", "landCover", "spi")

# Average values over each polygon

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

MyMedian <- function(x){
  outVar <- median(x, na.rm = T)
}

StackLS.df <- data.frame(matrix(nrow = length(allPolygons_spatial), ncol = 0))
for (jj in 1:nlayers(predictStack)){
  # Train
  tmpRas <- raster::extract(x = predictStack[[jj]], y = allPolygons_spatial)
  tmpRas.list <- lapply(tmpRas, MyMedian)
  tmpRas.df <- t(data.frame(tmpRas.list))
  StackLS.df <- cbind(StackLS.df, tmpRas.df)
}

colnames(StackLS.df) <- names(predictStack)

# Seperate raster stacks into LS train/test and NLS
stackNLS <- raster::mask(predictStack, NLS_DS)


# Convert seperated stacks into data frames
stackNLS.df <- na.omit(as.data.frame(stackNLS))

data.df <- rbind(StackLS.df, stackNLS.df)
data.mat <- data.matrix(data.df)




# Correlation using GVIF -------------------------------------------------------



# add the response variable
data.df$response <- c(rep(1, nrow(StackLS.df)), rep(0, nrow(stackNLS.df)))


y <- c(rep(1, nrow(StackLS.df)), rep(0, nrow(stackNLS.df)))
data.df$response <- y  # Add the binary column to dataframe


model <- glm(response ~ ., data = data.df, family = binomial)

vif(model)






# Categorical data generation and correlation using Spearman's correlation-------------------------------------------------------

predictStackCat <- stack(aspect, landCover, geology)
names(predictStackCat) <- c("aspect", "landCover", "geology")

predictStackCat.df <- na.omit(as.data.frame(predictStackCat))

dataCat.mat <- data.matrix(predictStackCat.df)

spearmanCorr <- Hmisc::rcorr(dataCat.mat, type = "spearman")
corrplot::corrplot(spearmanCorr$r, type = "lower", method = "color")


# View continuous variable correlation using Pearson's correlation-----------------------------------------------

# Check numerical variables using pearson correlation
pearsonCorr <- Hmisc::rcorr(data.mat, type = "pearson")

# Plot correlation values
corrplot::corrplot(pearsonCorr$r, type = "lower", method = "color")

# Plot p-values associated with pearson correlation (replace NA with 0)
pVals <- pearsonCorr$P
pVals[is.na(pVals)] <- 0
corrplot::corrplot(pVals, type = "lower", method = "color")

## Find high correlation indices
# https://www.sciencedirect.com/science/article/pii/S0048969720307415 (source for 0.7 threshold)

corrplot::corrplot((pearsonCorr$r > 0.7) | (pearsonCorr$r < -0.7), type = "lower", method = "color")


# Check if highly correlated varaibles have p-values less than 0.05
pVals_logical <- pearsonCorr$P < 0.05
pVals_logical[is.na(pVals_logical)] <- 0
corrplot::corrplot(pVals_logical, type = "lower", method = "color")

round(pearsonCorr$r, 2)
which(abs(pearsonCorr$r) > 0.7 & pearsonCorr$r != 1, arr.ind = TRUE)


corrplot::corrplot(pearsonCorr$r,
                   type = "lower",
                   method = "color",
                   addCoef.col = "black",
                   number.cex = 0.7,
                   tl.col = "red",
                   tl.srt = 45)  


