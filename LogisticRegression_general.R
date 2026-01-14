# General suscept analysis

#Package Installation ----------------------------------------------------

install.packages("sp")
install.packages("raster")
install.packages("caret")
install.packages("sf")
install.packages("glmnet")
install.packages("tibble")
install.packages("ranger")
install.packages(c("mlr3", "mlr3pipelines", "mlr3learners"))
install.packages("mlr3viz")
install.packages("mlr3mbo")
install.packages("mlr3spatial")
install.packages("mlr3pipelines")
install.packages('future')
install.packages('precrec')
install.packages("dplyr")
install.packages('ggplot2')
install.packages('lattice')
install.packages('paradox')
install.packages("terra")
install.packages("remotes")
install.packages("parallelly")
install.packages("mlr3spatial")
install.packages("mlr3")  # si no lo habías instalado ya
library(mlr3)


remotes::install_github("jpconnel/mlr3spatial", force = TRUE)
# remotes::install_github("mlr-org/mlr3spatial@prob")



## Install required modules
library(sp)
library(raster)
library(caret)
library(sf)
library(glmnet)
library(tibble)
library(ranger)
library(mlr3learners)
library(mlr3)
library(mlr3viz)
library(mlr3mbo)
library(mlr3tuning)
library(mlr3spatial)
library(mlr3pipelines)
library(future)
library(precrec)
library(dplyr)
library(terra)
library(remotes)

# Thematic variable import (12.5 m) ------------------------------------------

plCurv <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\plancurv.tif")

slope <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\slope.tif")
spi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\SPI.tif")
tpi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\TPI_clasif6cell.tif")

#pga <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\PGAA.tif")
#pgv <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\PGV.tif")



# Import categorical data
geology <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\geology.tif")
aspect <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\aspect.tif")




# Re-classify non-linear numeric vars -----------

# https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/pretty

facts <- list("geology", "tpi")

aspectCriteria <- matrix(c(-5, -1, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, -1, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 380, 9
                           , 1, 2, 3, 4, 5, 6, 7, 8, 1), ncol = 3)
aspect <- reclassify(aspect, aspectCriteria)
facts <- append(facts, "aspect")



#TPI criteria applying Weiss 2001
#tpi_class <- tpi20cell # create empty raster
#tpi_class[] <- NA  # applying conditions, starting with NA
#tpi_class[tpi20cell <= -1] <- 1 # 1. valley (TPI <= -1)
#tpi_class[tpi20cell > -1 & tpi20cell <= -0.5] <- 2 # 2. lower slope (-1 < TPI <= -0.5)
#tpi_class[tpi20cell > -0.5 & tpi20cell <= 0.5 & slope <= 5] <- 3 # 3. flat slope (-0.5 < TPI <= 0.5 and slope <= 5)
#tpi_class[tpi20cell > -0.5 & tpi20cell <= 0.5 & slope > 5] <- 4 # 4. middle slope (-0.5 < TPI <= 0.5 and slope > 5)
#tpi_class[tpi20cell > 0.5 & tpi20cell <= 1] <- 5 # 5. upper slope (0.5 < TPI <= 1)
#tpi_class[tpi20cell > 1] <- 6 # 6. ridge (TPI > 1)
#writeRaster(tpi_class, "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Logistic regression April 2025\\TPIs\\TPI_clasif20cell.tif", overwrite = TRUE)

plCurvCriteria <- matrix(c(-300, -0.05, 0.05, -0.05, 0.05, 300, 2, 1, 3), ncol = 3)
plCurv <- raster::reclassify(plCurv, plCurvCriteria)
facts <- append(facts, "plCurv")

slopeCriteria <- matrix(
  c(-Inf, 5, 1,
    5, 10, 2,
    10, 15, 3,
    15, 20, 4,
    20, 25, 5,
    25, 30, 6,
    30, 45, 7,
    45, Inf, 8),
  ncol = 3, byrow = TRUE
)
slope <- raster::reclassify(slope, slopeCriteria)

# # twi (inverse)
#twi <- 1 / (twi + 1)

# Import Landslide Database -------------------------------------

# Open polygons from shapefile

ls_Polygons <- sf::st_read("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Inventory and sources\\sources_may2025utm.shp")
lsPolygons_sf <- ls_Polygons$geometry


lsPolygons_sf <- ls_Polygons$geometry
lsCentroids <- sf::st_point_on_surface(lsPolygons_sf)
allPolygons <- sf::st_zm(lsPolygons_sf)
allPoly_spatial <- sf::as_Spatial(allPolygons, cast = TRUE)

# Initialize stack and calculate landslide properties -------------

# Create an initial stack of all pixels in study area

print("Initializing stack")
# JASON ADDED SPI
predictStack <- stack(aspect, tpi, plCurv, slope, geology, spi, pga)
names(predictStack) <- c("aspect", "tpi", "plCurv", "slope", "geology", "spi", "pga")


# Find representative values for each landslide polygon


stackLS.df <- data.frame(matrix(nrow = length(allPoly_spatial), ncol = 0))
rownames(stackLS.df) <- names(allPoly_spatial)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

MyMedian <- function(x){
  outVar <- median(x, na.rm = T)
}  #check the Median


for (jj in 1:nlayers(predictStack)){
  if (names(predictStack)[jj] %in% facts){
    
    tmpRas <- raster::extract(x = predictStack[[jj]], y = allPoly_spatial)
    tmpRas.list <- lapply(tmpRas, Mode)
    tmpRas.df <- t(data.frame(tmpRas.list))
    stackLS.df <- cbind(stackLS.df, tmpRas.df)
    
    
  }
  else {
    
    tmpRas <- raster::extract(x = predictStack[[jj]], y = allPoly_spatial)
    tmpRas.list <- lapply(tmpRas, MyMedian)
    tmpRas.df <- t(data.frame(tmpRas.list))
    stackLS.df <- cbind(stackLS.df, tmpRas.df)
    
  }
}

colnames(stackLS.df) <- names(predictStack)

# Determine number of landslides
lsCount <- length(lsPolygons_sf)

rasterLS <- rasterize(allPoly_spatial, aspect, field = 1)
ls_full <- reclassify(rasterLS, cbind(NA, 0)) 
ls_full <- raster::mask(x = ls_full, mask = aspect)
allNLS <- reclassify(ls_full, cbind(1, NA))

# Perform statistical analysis -------------





# Initiate resampling loop
learnersList <- list()
predictionsList <- list()
dataList <- list()
nRepeats <- 100  #change this, according the number of repeats I want
outerResamp <- rsmp("subsampling", ratio = 0.7, repeats = nRepeats)


for (ii in 1:nRepeats){
  
  NLS_DS <- sampleRandom(allNLS, lsCount * 5, asRaster = TRUE)
  
  stackNLS <- raster::mask(predictStack, NLS_DS) # Changed from NLS_DS
  stackNLS.df <- na.omit(as.data.frame(stackNLS))
  
  # weightsLayer <- reclassify(ls_full, cbind(1,5))
  # weightsLayer <- reclassify(weightsLayer, cbind(0,1))
  # predictStack_v2 <- stack(ls_full, predictStack, weightsLayer)
  # names(predictStack_v2)[1] <- "landslides"
  # names(predictStack_v2)[nlayers(predictStack_v2)] <- "weights"
  
  # predict.df <- na.omit(as.data.frame(predictStack_v2))
  # predict.df$landslides <- as.factor(predict.df$landslides)
  
  #stackLS.df <- na.omit(as.data.frame()) # Adding this line to match the structure
  NLS_vec <- 1:nrow(stackNLS.df)
  
  data.df <- rbind(stackLS.df, stackNLS.df)
  
  # Create landslide and weights vectors ---------------------------
  
  outcome.vec <- c(rep(1, length(stackLS.df$aspect)), rep(0, length(stackNLS.df$aspect)))
  
  # Weight of 5 applied
  wts.vec <- c(rep(5, length(stackLS.df$aspect)), rep(1, length(stackNLS.df$aspect)))
  
  data.df <- cbind(outcome.vec, data.df, wts.vec)
  names(data.df)[1] <- 'landslides'
  names(data.df)[length(names(data.df))] <- "weights"
  data.df$landslides <- as.factor(data.df$landslides)
  
  
  # Convert Categorical variables to factors ---------------------
  
  data.df$aspect <- as.factor(data.df$aspect)
  #data.df$relief <- as.factor(data.df$relief)
  #data.df$waterDist <- as.factor(data.df$waterDist)
  #data.df$faultDist <- as.factor(data.df$faultDist)
  #data.df$twi <- as.factor(data.df$twi)
  #data.df$map <- as.factor(data.df$map)
  data.df$tpi <- as.factor(data.df$tpi)
  #data.df$prCurv <- as.factor(data.df$prCurv)
  data.df$plCurv <- as.factor(data.df$plCurv)
  data.df$slope <- as.factor(data.df$slope)
  data.df$geology <- as.factor(data.df$geology)
  #data.df$elevation <- as.factor(data.df$elevation)
  # data.df$pga <- as.factor(data.df$pga)
  # data.df$pgv <- as.factor(data.df$pgv)
  # 
  
  
  
  # data.df$trimlines <- as.factor(data.df$trimlines)
  # data.df$roadsDist <- as.factor(data.df$roadsDist)
  
  #----------------------------------------------------------------------
  
  task_lr <- TaskClassif$new('task_lr', backend = data.df, target = 'landslides')
  task_lr$set_col_roles(cols = 'weights', roles = 'weights_learner')
  task_lr$positive <- '1'
  
  dataList[[ii]] <- data.df
  outerResamp$instantiate(task_lr) # SWITCH
  
  
  graph_lr <- po("encode", method = 'treatment') %>>% po("scale") %>>% lrn("classif.cv_glmnet",
                                                                           predict_type = 'prob', 
                                                                           type.measure = 'auc', 
                                                                           predict_sets = c("train", "test"))
  
  # graph_lr = mlr3pipelines::po("scale", param_vals = list(method = "treatment"))
  
  
  
  
  graph_lr$param_set$values <- mlr3misc::insert_named(graph_lr$param_set$values, 
                                                      list(classif.cv_glmnet.standardize.response = FALSE))
  # graph_lr$keep_results = TRUE
  
  graphLearner_lr <- GraphLearner$new(graph_lr)
  graphLearner_lr$predict_sets <- c("test", "train")
  
  learnerClone <- graphLearner_lr$clone(deep = TRUE)
  learnersList[[ii]] <- learnerClone$train(task_lr, row_ids = outerResamp$train_set(ii))
  predictionsList[[ii]] <- list(test = learnerClone$predict(task_lr, row_ids = outerResamp$test_set(ii)), 
                                train = learnerClone$predict(task_lr, row_ids = outerResamp$train_set(ii)))
  trainPred <- learnersList[[ii]]$predict(task_lr, row_ids = outerResamp$train_set(ii))
  testPred <- learnersList[[ii]]$predict(task_lr, row_ids = outerResamp$test_set(ii))
  print(ii)
}


# End resampling loop
rdata <- as_result_data(task_lr, learnersList, predictions = predictionsList, 
                        resampling = outerResamp, iterations = as.integer(1:nRepeats))
lr_resample <- ResampleResult$new(rdata)
saveRDS(lr_resample, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_resample50") # SWITCH
saveRDS(dataList, file= "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_resample_DataList50")


# Extract model training performance metrics ------------------------------


# Training data fit (ENSURE THIS IS DIFFERENT FROM TESTING)  SKIPPABLE
trainingResults <- lr_resample$score(measures = c(msr('classif.auc', id = "auc.train", predict_sets = "train"), 
                                                  msr('classif.bbrier', id = "bbrier.train", predict_sets = "train"), 
                                                  msr('classif.dor', id = "dor.train", predict_sets = "train"), 
                                                  msr('classif.tn', id = "tn.train", predict_sets = "train"), 
                                                  msr('classif.tp', id = "tp.train", predict_sets = "train"), 
                                                  msr('classif.fn', id = "fn.train", predict_sets = "train"), 
                                                  msr('classif.fp', id = "fp.train", predict_sets = "train"), 
                                                  msr('classif.ce', id = "ce.train", predict_sets = "train"), 
                                                  msr('classif.fbeta', id = "fbeta.train", predict_sets = "train"), 
                                                  msr('classif.logloss', id = "logloss.train", predict_sets = "train"), 
                                                  msr('classif.bacc', id = "bacc.train", predict_sets = "train")))
head(trainingResults)
saveRDS(trainingResults, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_trainingMetrics_cedar_dalex50")


# Extract model testing performance data ----------------------------------

# Test Results
testingResults <- lr_resample$score(measures = c(msr('classif.auc', id = "auc.test", predict_sets = "test"), 
                                                 msr('classif.bbrier', id = "bbrier.test", predict_sets = "test"), 
                                                 msr('classif.dor', id = "dor.test", predict_sets = "test"), 
                                                 msr('classif.tn', id = "tn.test", predict_sets = "test"), 
                                                 msr('classif.tp', id = "tp.test", predict_sets = "test"), 
                                                 msr('classif.fn', id = "fn.test", predict_sets = "test"), 
                                                 msr('classif.fp', id = "fp.test", predict_sets = "test"), 
                                                 msr('classif.ce', id = "ce.test", predict_sets = "test"), 
                                                 msr('classif.fbeta', id = "fbeta.test", predict_sets = "test"), 
                                                 msr('classif.logloss', id = "logloss.test", predict_sets = "test"), 
                                                 msr('classif.bacc', id = "bacc.test", predict_sets = "test"), 
                                                 msr('classif.acc', id = "acc.test", predict_sets = "test")))
head(testingResults)
saveRDS(testingResults, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_testingMetrics_cedar_dalex50") # SWITCH


# Extract model coefficients ----------------------------------------------



lambda_1se <- vector()
betaCoeffs <- data.frame(matrix(nrow = nrow(lr_resample$learners[[1]]$model$classif.cv_glmnet$model$glmnet.fit$beta), ncol = 0))
rownames(betaCoeffs) <- lr_resample$learners[[1]]$model$encode$outtasklayout$id
a0vals <- vector()

for (ii in 1:length(lr_resample$learners)){
  indivLrn <- lr_resample$learners[[ii]]
  lambda_1se <- append(lambda_1se, indivLrn$model$classif.cv_glmnet$model$lambda.1se)  
  lambdaIndex <- na.omit(match(lambda_1se[ii], indivLrn$model$classif.cv_glmnet$model$lambda))
  a0vals <- append(a0vals, indivLrn$model$classif.cv_glmnet$model$glmnet.fit$a0[lambdaIndex])
  
  indivBeta <- indivLrn$model$classif.cv_glmnet$model$glmnet.fit$beta[,lambdaIndex]
  if (ii == 1){
    betaCoeffs <- data.frame(names(indivBeta), indivBeta)
    colnames(betaCoeffs) <- c("names", "1")
    
  }
  else{
    betaAddOn <- data.frame(names(indivBeta), indivBeta)
    names(betaAddOn) <- c("names", ii)
    betaCoeffs <- full_join(betaCoeffs, betaAddOn, by="names")
  }
}
saveRDS(betaCoeffs, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_resample_betas50") # SWITCH
saveRDS(a0vals, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_resample_fs_a050") # SWITCH
saveRDS(lambda_1se, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_resample_lambda1se50") # SWITCH



#___________________________ROC CURVE______________________________


install.packages("pROC")
install.packages("ggplot2")
install.packages("dplyr")
library(pROC)
library(ggplot2)
library(dplyr)

roc_list <- list()

for (i in seq_along(predictionsList)) {
  pred <- predictionsList[[i]]$test
  
  if (!"1" %in% colnames(pred$prob)) next
  
  scores <- as.numeric(pred$prob[, "1"])
  labels <- as.numeric(as.character(pred$truth))
  
  if (length(scores) == length(labels) && all(c(0,1) %in% labels)) {
    roc_obj <- roc(response = labels, predictor = scores, quiet = TRUE)
    roc_list[[i]] <- roc_obj
  } else {
    cat("️  Iteration", i, "wthout both clases\n")
  }
}

fpr_grid <- seq(0, 1, length.out = 100)
tpr_matrix <- matrix(NA, nrow = length(roc_list), ncol = length(fpr_grid))

for (i in seq_along(roc_list)) {
  roc_obj <- roc_list[[i]]
  interp_tpr <- approx(1 - roc_obj$specificities, roc_obj$sensitivities, xout = fpr_grid)$y
  tpr_matrix[i, ] <- interp_tpr
}


mean_tpr <- apply(tpr_matrix, 2, mean, na.rm = TRUE)
sd_tpr <- apply(tpr_matrix, 2, sd, na.rm = TRUE)


roc_df <- data.frame(
  FPR = fpr_grid,
  TPR_mean = mean_tpr,
  TPR_upper = mean_tpr + sd_tpr,
  TPR_lower = mean_tpr - sd_tpr
)

#Plot
ggplot(roc_df, aes(x = FPR, y = TPR_mean)) +
  geom_line(color = "blue", size = 1.2) +
  geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper), fill = "blue", alpha = 0.2) +
  labs(title = "ROC average curve",
       x = "False positive rate (FPR)",
       y = "True positive rate (TPR)") +
  theme_minimal()


#AUC VALUE
auc_vals <- sapply(roc_list, auc)
mean_auc <- mean(auc_vals)
cat("Average AUC:", round(mean_auc, 3), "\n")



#betas___________________________________
#Graph 1 : vertical

# Load the beta coefficients from the RDS file.
betaCoeffs <- readRDS("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\lr_resample_betas50")

# Check the first coefficients to ensure they were loaded correctly
head(betaCoeffs)

# Calculate the average of the beta coefficients for each variable
betaCoeffs_mean <- rowMeans(betaCoeffs[,-1], na.rm = TRUE)  # Ignore the "names" column and average the remaining ones

# Create a dataframe with the variable names and their average coefficients
betaCoeffs_df <- data.frame(variable = betaCoeffs$names, coef_mean = betaCoeffs_mean)

# Display the variables sorted by the absolute value of their average coefficients
betaCoeffs_df <- betaCoeffs_df[order(abs(betaCoeffs_df$coef_mean), decreasing = TRUE), ]
print(betaCoeffs_df)


# Load ggplot2
library(ggplot2)

ggplot(betaCoeffs_df, aes(x = reorder(variable, coef_mean), y = coef_mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(coef_mean, 2)), 
            hjust = ifelse(betaCoeffs_df$coef_mean >= 0, -0.1, 1.1),  # cambia posición si es positivo o negativo
            size = 3.5) +
  coord_flip() +
  labs(title = "Variable Importance in the Model",
       x = "Variables",
       y = "Average Coefficients (\u03B2)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # centrar título
  expand_limits(y = c(min(betaCoeffs_df$coef_mean) * 1.2, max(betaCoeffs_df$coef_mean) * 1.2))  # espacio para el texto




#GRAPH 2: absolute values and red,blue labels, vertical

betaCoeffs_mean <- rowMeans(betaCoeffs[,-1], na.rm = TRUE)  # Ignore name column 
betaCoeffs_df <- data.frame(variable = betaCoeffs$names, coef_mean = betaCoeffs_mean)

# Transform to absolute values
betaCoeffs_df$coef_abs <- abs(betaCoeffs_df$coef_mean)

# Labels: "Increase" if positive, "Decrease" if negative
betaCoeffs_df$Effect <- ifelse(betaCoeffs_df$coef_mean > 0, "Increases landslide classification", "Decreases landslide classification")

# Order variables by magnitude of coefficients (absolute value)
betaCoeffs_df <- betaCoeffs_df %>% arrange(desc(coef_abs))

# Plot
ggplot(betaCoeffs_df, aes(x = reorder(variable, coef_abs), y = coef_abs, fill = Effect)) +
  geom_bar(stat = "identity") +  # Bars
  coord_flip() +  # Rotate plot (horizontal)
  scale_fill_manual(values = c("Increases landslide classification" = "#d73027", "Decreases landslide classification" = "#4575b4")) + 
  labs(title = "Variable Contribution to Landslide Susceptibility",
       x = "Variables",
       y = "Absolute Beta Coefficients (|β|)",
       fill = NULL) +  
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top")



#Graph 3, horizontal, red/blue labels


betaCoeffs_mean <- rowMeans(betaCoeffs[,-1], na.rm = TRUE)  
betaCoeffs_df <- data.frame(variable = betaCoeffs$names, coef_mean = betaCoeffs_mean)
betaCoeffs_df$coef_abs <- abs(betaCoeffs_df$coef_mean)

betaCoeffs_df$Effect <- ifelse(betaCoeffs_df$coef_mean > 0, "Increases landslide classification", "Decreases landslide classification")

betaCoeffs_df <- betaCoeffs_df %>% arrange(desc(coef_abs))

ggplot(betaCoeffs_df, aes(x = variable, y = coef_abs, fill = Effect)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("Increases landslide classification" = "#d73027", "Decreases landslide classification" = "#4575b4")) + 
  labs(title = "Variable Contribution to Landslide Susceptibility",
       x = "Variables",
       y = "Absolute Beta Coefficients (|β|)",
       fill = NULL) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Rotate labels axis x
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top")





# Generate susceptibility maps ----------------------------------------------





library(raster)
spatialPredStack <- stack()

for (jj in seq_along(lr_resample$learners)) {
  tmpLearner <- lr_resample$learners[[jj]]
  tmpPred <- predict_spatial(newdata = predictStack, learner = tmpLearner, format = 'raster')
  writeRaster(tmpPred, filename = sprintf("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 2 filtered\\Predictions_results_test%02d.tif", jj), overwrite = TRUE)
  print(jj)
}


#__________AVERAGE RASTER__________


# tiffs route
ruta <- "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Review_Logistic regression May 2025\\Results - Model 4 PGAA\\"
# List of the TIF files in the directory
file_list <- list.files(path = ruta, pattern = "\\.tif$", full.names = TRUE)


if(length(file_list) > 0){
  # Load all rasters
  stacked_rasters <- stack(file_list)
  
  # average value of rasters
  mean_raster <- calc(stacked_rasters, fun = mean)
  
  # save the raster
  output_file <- file.path(ruta, "mean_genersuscept_raster.tif")
  writeRaster(mean_raster, filename = output_file, format="GTiff", overwrite=TRUE)
  
  print(paste("The average raster has been saved in:", output_file))
} else {
  print("No TIFs files in the directory")
}



