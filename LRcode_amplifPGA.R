# Amplified PGA 

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
ndvi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\NDVI.tif")
slope <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\slope.tif")
spi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\SPI.tif")
tpi <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\TPI_clasif6cell.tif") #first run with 6 cells
pga <- raster('C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\PGAA.tif')
# Import categorical data
geology <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\geology.tif")
aspect <- raster("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\aspect.tif")


cat("NDVI:\n");    print(summary(values(ndvi)));    print(unique(values(ndvi)))
cat("Slope:\n");   print(summary(values(slope)));   print(unique(values(slope)))
cat("SPI:\n");     print(summary(values(spi)));     print(unique(values(spi)))
cat("TPI:\n");     print(summary(values(tpi)));     print(unique(values(tpi)))
cat("Geology:\n"); print(summary(values(geology))); print(unique(values(geology)))
cat("Aspect:\n");  print(summary(values(aspect)));  print(unique(values(aspect)))
cat("PGA:\n");     print(summary(values(pga)));     print(unique(values(pga)))

library(raster)

# Apila los raster
rasters <- stack(ndvi, slope, spi, tpi, geology, aspect, pga)

# Crea una máscara: NA donde haya NA en cualquier capa
mask_all <- calc(rasters, fun=function(x) if(any(is.na(x))) NA else 1)

# Aplica la máscara a cada capa
for(i in 1:nlayers(rasters)) {
  rasters[[i]] <- mask(rasters[[i]], mask_all)
}

# Split layers
ndvi    <- rasters[[1]]
slope   <- rasters[[2]]
spi     <- rasters[[3]]
tpi     <- rasters[[4]]
geology <- rasters[[5]]
aspect  <- rasters[[6]]
pga  <- rasters[[7]]

# Re-classify non-linear numeric vars -----------

# https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/pretty

facts <- list("geology", "tpi")

aspectCriteria <- matrix(c(-5, -1, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, -1, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 380, 9
                           , 1, 2, 3, 4, 5, 6, 7, 8, 1), ncol = 3)
aspect <- reclassify(aspect, aspectCriteria)
facts <- append(facts, "aspect")

ndviCriteria <- matrix(c(-1, 0, 0.2, 0.5, 0.7, 0, 0.2, 0.5, 0.7, 1, 1, 2, 3, 4, 5), ncol = 3)
ndvi <- raster::reclassify(ndvi, ndviCriteria)
facts <- append(facts, "ndvi")

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

unique(values(slope))

# Apila los rasters reclasificados
rasters <- stack(ndvi, slope, spi, tpi, geology, aspect, pga)

# Suma los NA por celda
na_count <- calc(rasters, function(x) sum(is.na(x)))

# Verifica que todas las celdas con NA tengan NA en todas las capas
table(values(na_count))

# Import Landslide Database -------------------------------------

# Open polygons from shapefile

ls_Polygons <- sf::st_read("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\Inventory_sources\\sources.shp")

lsPolygons_sf <- ls_Polygons$geometry
lsCentroids <- sf::st_point_on_surface(lsPolygons_sf)
allPolygons <- sf::st_zm(lsPolygons_sf)
allPoly_spatial <- as(allPolygons, "Spatial")

# Initialize stack and calculate landslide properties -------------

# Create an initial stack of all pixels in study area

print("Initializing stack")

predictStack <- stack(aspect, ndvi, tpi, slope, geology, spi, pga)
names(predictStack) <- c("aspect", "ndvi", "tpi", "slope", "geology", "spi", "pga")
print(names(predictStack))

# Find representative values for each landslide polygon

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

MyMedian <- function(x){
  median(x, na.rm = TRUE)
}


print("Extracting values for each polygon")
ls_columns <- list()
n_poly <- length(allPoly_spatial)
for (jj in 1:nlayers(predictStack)){
  layer_name <- names(predictStack)[jj]
  tmpRas <- raster::extract(x = predictStack[[jj]], y = allPoly_spatial)
  # This ensures every polygon gets a value, even if it's NA
  if (layer_name %in% facts){
    tmpRas.list <- lapply(tmpRas, function(x) if (length(x) == 0) NA else Mode(x))
  } else {
    tmpRas.list <- lapply(tmpRas, function(x) if (length(x) == 0) NA else MyMedian(x))
  }
  ls_columns[[layer_name]] <- unlist(tmpRas.list)
}

stackLS.df <- data.frame(
  aspect = ls_columns$aspect,
  ndvi = ls_columns$ndvi,
  tpi = ls_columns$tpi,
  slope = ls_columns$slope,
  geology = ls_columns$geology,
  spi = ls_columns$spi,
  pga = ls_columns$pga
)

print("Dataframe structure")
str(stackLS.df)
print("Dataframe summary:")
summary(stackLS.df)
str(stackLS.df)
summary(stackLS.df)


# Optional: check that all columns are the same length
for (name in names(ls_columns)) {
  cat(name, "length:", length(ls_columns[[name]]), "\n")
}
stackLS.df <- as.data.frame(ls_columns)


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
nRepeats <- 50   #change this, according the number of repeats I want
outerResamp <- rsmp("subsampling", ratio = 0.7, repeats = nRepeats)

aspect_levels <- as.character(1:9)      # reclass values
ndvi_levels   <- as.character(1:5)
tpi_levels    <- as.character(1:6)
slope_levels  <- as.character(1:8)
geology_levels<- as.character(1:6)
pga_levels<- as.character(1:7)


print(aspect_levels)
print(ndvi_levels)
print(tpi_levels)
print(slope_levels)
print(geology_levels)
print(pga_levels)

for (ii in 1:nRepeats){
  
  NLS_DS <- sampleRandom(allNLS, lsCount * 5, asRaster = TRUE)
  
  stackNLS <- raster::mask(predictStack, NLS_DS)
  stackNLS.df <- na.omit(as.data.frame(stackNLS))
  NLS_vec <- 1:nrow(stackNLS.df)
  
  print(colnames(stackLS.df))
  print(colnames(stackNLS.df))
  
  data.df <- rbind(stackLS.df, stackNLS.df)
  
  # Create landslide and weights vectors
  outcome.vec <- c(rep(1, length(stackLS.df$aspect)), rep(0, length(stackNLS.df$aspect)))
  wts.vec <- c(rep(5, length(stackLS.df$aspect)), rep(1, length(stackNLS.df$aspect)))
  data.df <- cbind(outcome.vec, data.df, wts.vec)
  names(data.df)[1] <- 'landslides'
  names(data.df)[length(names(data.df))] <- "weights"
  data.df$landslides <- as.factor(data.df$landslides)
  
  # --- AQUÍ: Forzar niveles de los factores ---
  data.df$aspect <- factor(data.df$aspect, levels = aspect_levels)
  data.df$ndvi <- factor(data.df$ndvi, levels = ndvi_levels)
  data.df$tpi <- factor(data.df$tpi, levels = tpi_levels)
  data.df$slope <- factor(data.df$slope, levels = slope_levels)
  data.df$geology <- factor(data.df$geology, levels = geology_levels)
  
  cat("\n--- Chequeo de niveles y NA en la iteración", ii, "---\n")
  
  cat("Aspect levels in data.df:\n")
  print(table(data.df$aspect, useNA = "ifany"))
  cat("NDVI levels in data.df:\n")
  print(table(data.df$ndvi, useNA = "ifany"))
  cat("TPI levels in data.df:\n")
  print(table(data.df$tpi, useNA = "ifany"))
  cat("Slope levels in data.df:\n")
  print(table(data.df$slope, useNA = "ifany"))
  cat("Geology levels in data.df:\n")
  print(table(data.df$geology, useNA = "ifany"))
  
  # Chequeo de columnas con solo NA
  cols_to_check <- setdiff(names(data.df), c("landslides", "weights"))
  na_cols <- sapply(data.df[, cols_to_check], function(x) all(is.na(x)))
  if (any(na_cols)) {
    cat("Columnas con solo NA en esta iteración:\n")
    print(names(na_cols[na_cols]))
  }
  
  # Chequeo de dimensiones y nombres de columnas
  cat("Columnas en data.df:\n")
  print(colnames(data.df))
  cat("Tabla de landslides (presencia/ausencia):\n")
  print(table(data.df$landslides))
  cat("Dimensiones de data.df:\n")
  print(dim(data.df))
  
  
  
  # Eliminar columnas con solo NA (excepto landslides y weights)
  cols_to_check <- setdiff(names(data.df), c("landslides", "weights"))
  na_cols <- sapply(data.df[, cols_to_check], function(x) all(is.na(x)))
  data.df <- data.df[, c("landslides", cols_to_check[!na_cols], "weights")]
  
  # --------------------------------------------
  
  # El resto de tu código sigue igual
  task_lr <- TaskClassif$new('task_lr', backend = data.df, target = 'landslides')
  task_lr$set_col_roles(cols = 'weights', roles = 'weight')
  task_lr$positive <- '1'
  
  dataList[[ii]] <- data.df
  outerResamp$instantiate(task_lr)
  
  graph_lr <- po("encode", method = 'treatment') %>>% po("scale") %>>% lrn("classif.cv_glmnet",
                                                                           predict_type = 'prob', 
                                                                           type.measure = 'auc', 
                                                                           predict_sets = c("train", "test"))
  graph_lr$param_set$values <- mlr3misc::insert_named(graph_lr$param_set$values, 
                                                      list(classif.cv_glmnet.standardize.response = FALSE))
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
saveRDS(lr_resample, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions 2\\lr_resample3") # SWITCH
saveRDS(dataList, file= "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions 2\\lr_resample_DataList3")


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
saveRDS(trainingResults, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions 2\\lr_trainingMetrics_cedar_dalex3")


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
saveRDS(testingResults, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions 2\\lr_testingMetrics_cedar_dalex3") # SWITCH



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
saveRDS(betaCoeffs, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions 2\\lr_resample_betas3") # SWITCH
saveRDS(a0vals, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions 2\\lr_resample_fs_a03") # SWITCH
saveRDS(lambda_1se, file = "C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions 2\\lr_resample_lambda1se3") # SWITCH

betaCoeffs <- readRDS("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions\\lr_resample_betas")
# Plot coefficient distributions, etc.

write.csv(betaCoeffs, file = "C:/Users/Lenovo/OneDrive - Simon Fraser University (1sfu)/Documents/PhD work/Amplif PGA Analysis/Predictions/lr_resample_betas.csv", row.names = FALSE)
#write.csv(data.frame(a0 = a0vals), file = "C:/Users/Lenovo/OneDrive - Simon Fraser University (1sfu)/Documents/PhD work/Amplif PGA Analysis/Predictions/lr_resample_fs_a0.csv", row.names = FALSE)
#write.csv(data.frame(lambda_1se = lambda_1se), file = "C:/Users/Lenovo/OneDrive - Simon Fraser University (1sfu)/Documents/PhD work/Amplif PGA Analysis/Predictions/lr_resample_lambda1se.csv", row.names = FALSE)

# Load results
testingResults <- readRDS("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions\\lr_testingMetrics_cedar_dalex")
boxplot(testingResults$auc.test, main = "AUC (Test) across 100 iterations")

#___________________________ROC CURVE________ Need to check this

lr_resample <- readRDS("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\Amplif PGA Analysis\\Predictions\\lr_resample")


library(pROC)
library(ggplot2)

# Create a data frame for ggplot
roc_df <- data.frame(
  FPR = fpr_grid,
  mean_TPR = mean_tpr,
  lower_TPR = pmax(mean_tpr - sd_tpr, 0),
  upper_TPR = pmin(mean_tpr + sd_tpr, 1)
)

# Calculate mean AUC
mean_auc <- mean(sapply(rocs, auc))
mean_auc_text <- paste0("Mean AUC = ", round(mean_auc, 3))

# Plot with ggplot2
ggplot(roc_df, aes(x = FPR, y = mean_TPR)) +
  geom_line(size = 1.2, color = "blue") +
  geom_ribbon(aes(ymin = lower_TPR, ymax = upper_TPR), fill = "lightblue", alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Mean ROC Curve Across Iterations",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  annotate("text", x = 0.65, y = 0.1, label = mean_auc_text, size = 5, color = "black") +
  theme_minimal(base_size = 15) +
  coord_equal()




# Generate susceptibility maps ----------------------------------------------

aspect_levels  <- as.character(1:9)
ndvi_levels    <- as.character(1:5)
tpi_levels     <- as.character(1:6)
slope_levels   <- as.character(1:8)
geology_levels <- as.character(1:6)



library(raster)
spatialPredStack <- stack()

for (jj in seq_along(lr_resample$learners)) {
  tmpLearner <- lr_resample$learners[[jj]]
  tmpPred <- predict_spatial(newdata = predictStack, learner = tmpLearner, format = 'raster')
  writeRaster(tmpPred, filename = sprintf("C:/Users/Lenovo/OneDrive - Simon Fraser University (1sfu)/Documents/PhD work/Amplif PGA Analysis/Predictions 2/Predictions_results_test%03d.tif", jj), overwrite = TRUE)
  print(jj)
}

#betas___________________________________
#Graph 1 : vertical

# Load the beta coefficients from the RDS file.
betaCoeffs <- readRDS("C:\\Users\\Lenovo\\OneDrive - Simon Fraser University (1sfu)\\Documents\\PhD work\\General Suscept Analysis\\Predictions\\lr_resample_betas")

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



#__________AVERAGE RASTER__________


# tiffs route
ruta <- "C:/Users/Lenovo/OneDrive - Simon Fraser University (1sfu)/Documents/PhD work/Amplif PGA Analysis/Predictions 2"

# List of the TIF files in the directory
file_list <- list.files(path = ruta, pattern = "\\.tif$", full.names = TRUE)


if(length(file_list) > 0){
  # Load all rasters
  stacked_rasters <- stack(file_list)
  
  # average value of rasters
  mean_raster <- calc(stacked_rasters, fun = mean)
  
  # save the raster
  output_file <- file.path(ruta, "mean_PGAamplif_raster.tif")
  writeRaster(mean_raster, filename = output_file, format="GTiff", overwrite=TRUE)
  
  print(paste("The average raster has been saved in:", output_file))
} else {
  print("No TIFs files in the directory")
}




#_________________End of the process ___________