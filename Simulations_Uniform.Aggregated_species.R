##########################################
#               SIMULATIONS              #
#        UNIFORM/AGGREGATED SPECIES      #
##########################################


#Charge the packages required
library(terra)
library(dismo)
library(spatstat)
library(gstat)


# CREATE PREDICTIVE LAYERS ------------------------------------------------

 
set.seed(1234)

# We start creating a clime variable where we can control some parameters of
# these variables (mean, range of variation and spatial autocorrelation) 
#Create the function
generateVarb <- function(nrow, ncol, Vmean, Vsd, Vrange) {
  
  grd <- expand.grid(x=1:ncol, y=1:nrow)
  
  # Variable simulation data
  model_varb <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=Vmean, 
                      model=vgm(psil=Vsd**2, range=Vrange, model="Sph"),
                      nmax=40)
  
  variab <- predict(model_varb, newdata=grd, nsim=1)
  colnames(variab) <- c("x", "y", "value")
  
  
  variab_df <- data.frame(variab)
  variab_df
}

clim <- generateVarb(nrow=1000, ncol=1000, Vmean=15, Vsd=6, Vrange=450)

clim_rst <- rast(clim, type="xyz")
plot(clim_rst)

rm(clim)

# We create other layer, vegetation. With low spatial autocorrelation
vegt <- generateVarb(nrow=1000, ncol=1000, Vmean=8, Vsd=3, Vrange=10)

vegt_rst <- rast(vegt, type="xyz")
plot(vegt_rst)

rm(vegt)

#create a raster with the extension desired (we have yet the coordinates well)
r_base.t <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1, ymin=0, ymax=1)

#put this extension to the predictors created above
clim_rst@ptr[["extent"]] <- r_base.t@ptr[["extent"]]
vegt_rst@ptr[["extent"]] <- r_base.t@ptr[["extent"]]

#change the resolution of the layers.
#we employ these scales after in the Maxent predictions
#clime variables
clim_rst.100 <- aggregate(clim_rst, fact = 10, fun = mean)  #from 1000x1000 units to 100x100, joining by 10 units
clim_rst.10 <- aggregate(clim_rst, fact = 100, fun = mean)  #from 1000x1000 units to 10x10, joining by 100 units

plot(clim_rst)
plot(clim_rst.100)
plot(clim_rst.10)

#vegetation variables
vegt_rst.100 <- aggregate(vegt_rst, fact = 10, fun = mean, na.rm = TRUE)  #from 1000x1000 units to 100x100, joining by 10 units
vegt_rst.10 <- aggregate(vegt_rst, fact = 100, fun = mean, na.rm = TRUE)  #from 1000x1000 units to 10x10, joining by 100 units

plot(vegt_rst)
plot(vegt_rst.100)
plot(vegt_rst.10)

#standardize all the layers to allow comparison in regression and save on the computer
#vegetation
#scale the variable
rm(vegt_std.i)
rm(vegt_std_list)
vegt_std_list <- list()
for (i in 1:length(ls(pattern="vegt_rst*"))) {
  vegt_std.i <<- scale(lapply(ls(pattern="vegt_rst*"), get)[[i]])
  vegt_std_list <- c(vegt_std_list, vegt_std.i)
}
rm(i);rm(vegt_std.i)
#save the layer
for (i in 1:length(vegt_std_list)) {
  NUM <- as.character(nrow(vegt_std_list[[i]]))
  route_rst<-paste("mypath/VEGT/vegt",NUM,".asc", sep="")
  writeRaster(vegt_std_list[[i]], route_rst, NAflag = -9999, overwrite= TRUE)
}
rm(i);rm(NUM);rm(route_rst);rm(vegt_std_list)

#clima
#scale the variable
rm(clim_std.i)
rm(clim_std_list)
clim_std_list <- list()
for (i in 1:length(ls(pattern="clim_rst*"))) {
  clim_std.i <<- scale(lapply(ls(pattern="clim_rst*"), get)[[i]])
  clim_std_list <- c(clim_std_list, clim_std.i)
}
rm(i);rm(clim_std.i)
#save the layer
for (i in 1:length(clim_std_list)) {
  NUM <- as.character(nrow(clim_std_list[[i]]))
  route_rst<-paste("mypath/CLIM/clim",NUM,".asc", sep="")
  writeRaster(clim_std_list[[i]], route_rst, NAflag = -9999, overwrite= TRUE)
}
rm(i);rm(NUM);rm(route_rst);rm(clim_std_list)


# CREATE THE EXPECTED ABUNDANCE LAYERS ------------------------------------


#Create an error layer, with normal distribution, for the regression
r_error <- VEGT
r_error[] <- rnorm(ncell(r_error), mean = 0, sd = 0.1)
plot(r_error)


## UNIFORM PATTERN SPECIES -------------------------------------------------


#create an EXPECTED ABUNDANCE VALUE weighted by the CLIM and VEGT variables
ADQ <- 0.5*CLIM + 1*VEGT + r_error
plot(ADQ) #we can visualize

# transform the normal distribution of EXPECTED ABUNDANCE to the logistic distribution. 
# To obtain a range of positive values and considering a maximum of animals 
# present per cell  
ADQ_lgt <- logistic(ADQ, 1, 7, 3) 
par(mfrow = c(1,1))
plot(ADQ[], ADQ_lgt[]) #and visualize the relationship between the transformed and
#original values #TO SEE THE CURVE OF THE LOGISTIC FUNCTION

#The logistic distribution do not have 0, we plantain round values to obtain 0
#round the values to obtain zeros, or make a threshold of the bottom part 
ADQ_lgt_R <- round(ADQ_lgt, 2)


#Save the layers created
setwd("mypath")
route_rst<-paste("especialista/ADQ.asc", sep="")
writeRaster(ADQ, route_rst, NAflag = -9999, overwrite= TRUE)

route_rst<-paste("especialista/error_nodist.asc", sep="")
writeRaster(r_error, route_rst, NAflag = -9999, overwrite= TRUE)

route_rst<-paste("especialista/ADQ_lgt_R_nodist.asc", sep="")
writeRaster(ADQ_lgt_R, route_rst, NAflag = -9999, overwrite= TRUE)


## AGGREGATED PATTERN SPECIES ----------------------------------------------


#create an EXPECTED ABUNDANCE VALUE weighted by the CLIM and VEGT variables
ADQ <- 1*CLIM + 0.5*VEGT + r_error
plot(ADQ) #we can visualize

# transform the normal distribution of EXPECTED ABUNDANCE to the logistic distribution. 
# To obtain a range of positive values and considering a maximum of animals 
# present per cell  
ADQ_lgt <- logistic(ADQ, 1, 7, 3) 
par(mfrow = c(1,1))
plot(ADQ[], ADQ_lgt[]) #and visualize the relationship between the transformed and
#original values #TO SEE THE CURVE OF THE LOGISTIC FUNCTION

#The logistic distribution do not have 0, we plantain two forms to obtain 0
#round the values to obtain zeros, or make a threshold of the bottom part 
ADQ_lgt_R <- round(ADQ_lgt, 2)

#Save the layers created
route_rst<-paste("generalista/ADQ.asc", sep="")
writeRaster(ADQ, route_rst, NAflag = -9999, overwrite= TRUE)

route_rst<-paste("generalista/error_nodist.asc", sep="")
writeRaster(r_error, route_rst, NAflag = -9999, overwrite= TRUE)

route_rst<-paste("generalista/ADQ_lgt_R_nodist.asc", sep="")
writeRaster(ADQ_lgt_R, route_rst, NAflag = -9999, overwrite= TRUE)


# DISTRIBUTE ANIMALS ON SPACE -----------------------------------------------


#(DO IT FOR UNIFORM AND AGGREGATED SPATIAL PATTTERN SPECIES)
ADQ_lgt_R <- rast("mypath/especialista/ADQ_lgt_R_nodist.asc")
ADQ_lgt_R <- rast("mypath/generalista/ADQ_lgt_R_nodist.asc")


####Now, what we have to do, is obtaining simulated animals from this pattern of 
#heterogeneous LAMBDA 
r_base.1000 <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1000, ymin=0, ymax=1000, )
#put this extension to the rasters saved
ADQ_1000 <- ADQ_lgt_R
ADQ_1000@ptr[["extent"]] <- r_base.1000@ptr[["extent"]]

generatePoints <- function(x, min, max) {
  w <- owin(c(min,max), c(min,max))
  im.obj <- as.im(apply(matrix(x, max, max, byrow=TRUE), 2, rev), W=w)
  occurs <- rpoispp(im.obj)
  occurs
}

occurs <- generatePoints(ADQ_1000, 0, 1000)
par(mfrow = c(1,1))
plot(ADQ_lgt_R, asp = 1)
points(occurs, cex = 0.3, pch=19, col='red')


#Obtain a raster layer with the value of the true abundance of animals per cell 
occus <- vect(data.frame(x = occurs$x, y = occurs$y), 
              geom=c("x", "y")) #pass the occurrences to dataframe and spatvector
occurs_r <- rasterize(occus, ADQ_lgt_R, fun = sum)
plot(occurs_r, asp = 1)
#zoom()
plot(ADQ_lgt_R)
plot(occurs_r)
plot(ADQ_lgt_R, occurs_r)



# OBTAIN ABUNDANCE VALUES -------------------------------------------------


#create a gride of the extent that we want, and sample 400 of them
sam_grid <- rGridSampling(n= 400, xrange = c(0,1000), yrange = c(0,1000), 5)

#Plot the units sampled
#plot(ADQ_dist)
plot(ADQ_1000)
points(sam_grid, col='red', pch=16)

#Create a grid of all the units available to sample
rGrid <- function(xrange, yrange, cellsize) {
  s <- cellsize/2
  grd <- expand.grid(x=seq(xrange[1]+s, xrange[2]-s, s),
                     y=seq(yrange[1]+s, yrange[2]-s, s))
}

total_grid <- rGrid(xrange = c(0,1000), yrange = c(0,1000), 5)

ampl <- ext(700, 900, 500, 700)
zoom(ADQ_1000, ampl)
points(total_grid, col='blue', pch=16, cex=0.01)
points(sam_grid, col='red', pch=16)

#Extract the value of abundance of each sampled unit from the abundance layer
ABN_samp <- terra::extract(occurs_r, sam_grid, fun=sum, df=TRUE)
names(ABN_samp) <- c("SampleID", "Abundance")
ABN_samp[is.na(ABN_samp)] <- 0 #change NAs values per 0



# OBTAIN DETECTION NON-DETECTION VALUES -----------------------------------


#transform the abundance data into DETECTION/NON-DETECTION data, considering an individual
#imperfect detection, and considering if any animal have been detected per sample unit,
#considers the sample with or without detection
abun <- ABN_samp
ipd <- 0.5

sample_prove <- i_occuSample(abun, ipd)

#see the number of sampling units without and with detections
table(sample_prove$occurrences)
table(sample_prove$abundance)

#put a dataframe with the coordinates of samples, the abundances and the occurrences
sampling_units <- data.frame(cbind(x = sam_grid$x,
                                   y = sam_grid$y,
                                   sampleID = sample_prove$sampleID,
                                   abundance = sample_prove$abundance, 
                                   occurrences = sample_prove$occurrences))

#visualize the abundances
par(xpd=TRUE)
plot(ADQ_1000)
points(occurs, col = "red", pch = 16, cex = 0.05)
points(sampling_units, col="black", bg=sampling_units$abundance, pch=21)
legend("topright", inset=c(-0.2, 0), legend=c(sort(unique(sampling_units$abundance))), 
       col = "black", pt.bg = c(sort(unique(sampling_units$abundance))), bty = "n", pch = 21)

#visualize the detections no-detections
par(xpd=TRUE)
plot(ADQ_1000)
points(occurs, col = "red", pch = 16, cex = 0.05)
points(sampling_units, col="black", bg=sampling_units$occurrences, pch=21)
legend("topright", inset=c(-0.2, 0), legend=c(sort(unique(sampling_units$occurrences))), 
       col = "black", pt.bg = c(range(sampling_units$occurrences)), bty = "n", pch = 21)



# SAVE ALL DATA REQUIRED --------------------------------------------------


setwd("mypath")
#Save the data frame with all the information
write.table(sampling_units, 'generalista/sampling_units.csv', row.names = FALSE, sep = ",")
#And save as vector layer format
v.sampling_units <- vect(sampling_units, geom=c("x","y"))
terra::writeVector(v.sampling_units, 'generalista/sampling_units.shp', overwrite=TRUE)
#Save the data frame with the locations with detections only [This is what we are going to employ on Maxent]
final_data <- sampling_units[sampling_units$occurrences == 1,]
final_data <- data.frame(final_data[,c(1,2)])
h <- data.frame(ID = cbind(c(rep.int("Virtual_Species", nrow(final_data)))))
final_data <- list(ID=h, x=final_data$x, y=final_data$y)
write.table(final_data, "generalista/occ_species.csv", row.names = FALSE, sep = ",")

#We save the original abundance layer (all the values, not only sampled units)
terra::writeRaster(occurs_r, "generalista/true_abundance.asc", overwrite=TRUE)



# POSITIVE VALUES FOR MAXENT ----------------------------------------------


#Charge VEGETATION layers
setwd("mypath")
lista_rasters<-list()
for (j in 1:length(list.files("VEGT/", pattern = "^vegt.*\\.asc$"))) {
  lista_rasters[[j]]<-rast(paste("VEGT/", list.files("VEGT/", pattern = "^vegt.*\\.asc$")[j], sep=""))
}
rm(j)

lista_rasters

lista_post <- list()
for (i in 1:length(lista_rasters)) {
  lista_post[[i]] <- lista_rasters[[i]] + abs(lista_rasters[[3]]@ptr[["range_min"]]) +1
}
rm(i)

VEGT_10<-lista_post[[1]]
VEGT_100<-lista_post[[2]]
VEGT_1000<-lista_post[[3]]

#Save the layers
savepath <- file.path("mypath")
if (!dir.exists(savepath)) {
  dir.create(savepath, recursive=TRUE)
  print(paste("Creating directory", savepath))
} else {
  print(paste("Directory", savepath, "already exists. Overwriting."))
}

for(i in 1:length(ls(pattern = "VEGT_"))) {
  terra::writeRaster(get(ls(pattern = "VEGT_")[i]), 
                     paste0("mypath/VEGT", 
                            gsub("VEGT_", "", ls(pattern = "VEGT_")[i]),
                            ".asc"), NAflag = -9999, overwrite = TRUE)
}
rm(i)

#Charge CLIME layers
setwd("mypath")
lista_rasters<-list()
for (j in 1:length(list.files("CLIM/", pattern = "^clim.*\\.asc$"))) {
  lista_rasters[[j]]<-rast(paste("CLIM/", list.files("CLIM/", pattern = "^clim.*\\.asc$")[j], sep=""))
}
rm(j)

lista_rasters

lista_post <- list()
for (i in 1:length(lista_rasters)) {
  lista_post[[i]] <- lista_rasters[[i]] + abs(lista_rasters[[3]]@ptr[["range_min"]]) +1
}
rm(i)

CLIM_10<-lista_post[[1]]
CLIM_100<-lista_post[[2]]
CLIM_1000<-lista_post[[3]]

#Save the layers
savepath <- file.path(base, "mypath")
if (!dir.exists(savepath)) {
  dir.create(savepath, recursive=TRUE)
  print(paste("Creating directory", savepath))
} else {
  print(paste("Directory", savepath, "already exists. Overwriting."))
}

for(i in 1:length(ls(pattern = "CLIM_"))) {
  terra::writeRaster(get(ls(pattern = "CLIM_")[i]), 
                     paste0("mypath/CLIM", 
                            gsub("CLIM_", "", ls(pattern = "CLIM_")[i]),
                            ".asc"), NAflag = -9999, overwrite = TRUE)
}
rm(i)



# MAXENT SETTINGS ---------------------------------------------------------


#All data employed were on this repository
#mypath

#We employ MaxEnt software
#The parameters that have been taken into account were: 
#In SAMPLES load csv file of virtual species: occ_especialista y occ_generalista. One each time
#In ENVIRONMENTAL LAYERS load predictor layers considered for the model (CLIME and VEGETATION)
#The FEATURES that we have considered were: LINEAR, QUADRATIC and PRODUCT. 
#Options activated: CREATE RESPONSE CURVES; MAKE PICTURES OF PREDICTIONS; DO JACKKNIFE TO MEASURE VARIABLES IMPORTANCE
#Output format: LOGISTIC
#Output file type: ASC
#SETTINGS:
#(Basic)
#Random test percentage = 30
#Regularization multiplier = 1
#Max number of background points = 10000
#Replicates = 10
#Replicate run type = Subsample

#(Advanced)
#Bias file = None


# CREATE DIFFERENCE 10-100 DIMENSIONS -------------------------------------

#Join all the layers of the same resolution in one multilayer raster
COV.1000 <- c(CLIM_1000, VEGT_1000)
COV.100 <- c(CLIM_100, VEGT_100)
COV.10 <- c(CLIM_10, VEGT_10)
names(COV.1000) <- c("CLIM_1000", "VEGT_1000")
names(COV.100) <- c("CLIM_100", "VEGT_100")
names(COV.10) <- c("CLIM_10", "VEGT_10")

#Create a vector layer with all the centroid of spatial units
locs <- as.points(COV.1000$CLIM_1000)

#Extract for each location, the predictor values at 1000x1000, 100x100 and 10x10 dimensions
COV.ext.1000 <- terra::extract(COV.1000, locs)
COV.ext.100 <- terra::extract(COV.100, locs)
COV.ext.10 <- terra::extract(COV.10, locs)

library(plyr)
COV_TOT1 <- join(COV.ext.1000, COV.ext.100, by="ID")
COV_TOT2 <- join(COV.ext.100, COV.ext.10, by="ID")

#Calculate the distance (absolute value) between each COVARIATE at 100x100 and  
#1000x1000 dimensions for each location
#This is the DIFFERENCE between 1000x1000 and 100x100 dimensions
for (i in 2:3) {
  COV_TOT1[i+4] <- abs(COV_TOT1[i] - COV_TOT1[i+2])
  name <- names(COV_TOT1[i])
}
for (i in 2:3) {
  COV_TOT1[i+4] <- COV_TOT1[i] - COV_TOT1[i+2]
  name <- names(COV_TOT1[i])
}
names(COV_TOT1)
names(COV_TOT1) <- c("ID", "CLIM_1000", "VEGT_1000", "CLIM_100", "VEGT_100", "DIF_C.1000-100", "DIFF_V.1000-100")

#Calculate the distance between each COVARIATE at 10x10 and  
#100x100 dimensions for each location
#This is the DIFFERENCE between 100x100 and 10x10 dimensions
for (i in 2:3) {
  COV_TOT2[i+4] <- abs(COV_TOT2[i] - COV_TOT2[i+2])
  name <- names(COV_TOT2[i])
}
for (i in 2:3) {
  COV_TOT2[i+4] <- COV_TOT2[i] - COV_TOT2[i+2]
  name <- names(COV_TOT2[i])
}
names(COV_TOT2)
names(COV_TOT2) <- c("ID", "CLIM_100", "VEGT_100", "CLIM_10", "VEGT_10", "DIF_C.100-10", "DIFF_V.100-10")

#Join all the data together
COV_TOT <- join(COV_TOT1[,c(1,2,3,6,7)], COV_TOT2[,c(1,6,7)], by = "ID")

#Return again to vector format
v.COV_TOT <- cbind(locs, COV_TOT[,-1])
v.COV_TOT <- v.COV_TOT[,-1]

#Create a list for the rasterized layers
raster_layers <- list()

#Rasterize each column of the vector layer
for (i in names(v.COV_TOT)) {
  raster_layer <- rasterize(v.COV_TOT, COV.1000, field = i)
  raster_layers[[i]] <- raster_layer
}

#Combine all the layers in a multilayer spatraster
COV_raster <- rast(raster_layers)

#Visualize
plot(COV_raster)
names(COV_raster) <- c("CLIM_1000", "VEGT_1000", "DIF_C.1000-100", 
                       "DIF_V.1000-100", "DIF_C.100-10", "DIF_V.100-10")


# PREPARE A DATASET WITH ALL THE INFO -------------------------------------

#Charge the layers with the prediction at lower resolution for both species
Abund <- c(unif.species0, aggr.species0)
names(Abund) <- c("Unif.species", "Aggr.species")

#Charge Maxent predictions
setwd("mypath")
unif.pred <- rast("especialista/10/Virtual_Species_avg.asc")
aggr.pred <- rast("generalista/10/Virtual_Species_avg.asc")

unif.pred@cpp[["extent"]] <- COV_raster@cpp[["extent"]]
aggr.pred@cpp[["extent"]] <- COV_raster@cpp[["extent"]]
Abund@cpp[["extent"]] <- COV_raster@cpp[["extent"]]

#change to the resolution of the rest of layers
unif.pred_1000 <- terra::disagg(unif.pred, fact=100, method="near")
aggr.pred_1000 <- terra::disagg(aggr.pred, fact=100, method="near")

#And add the CRS to be equal in all the cases
crs(unif.pred_1000) <- crs(COV_raster)
crs(aggr.pred_1000) <- crs(COV_raster)
crs(Abund) <- crs(COV_raster)

#Join all the info: Predictions at lower resolution, true abundance, covariates
#of difference between resolutions (intermediate-scale covariates)
TOT.DATA <- c(Abund, unif.pred_1000, aggr.pred_1000, COV_raster)

names(TOT.DATA) <- c("Abund.unif", "Abund.aggr", "Pred.unif", "Pred.aggr", 
                     "CLIM_1000", "VEGT_1000", "DIF_C.1000-100", 
                     "DIF_V.1000-100", "DIF_C.100-10", "DIF_V.100-10")



# SUITABILITY-ABUNDANCE ANALYSIS ------------------------------------------

#prepare variables
Total.data.df <- as.data.frame(TOT.DATA)
data1 <- data.frame(scale(Total.data.df[,3:10]))
data2 <- data.frame(Total.data.df[,1:2])
model.data <- cbind(data2, data1)

names(Total.data.df) <- names(model.data)

#(DO IT FOR UNIFORM AND AGGREGATED SPATIAL PATTTERN SPECIES)
setwd("mypath")

#Create a data frame to save all the results
results <- data.frame(
  Replica = integer(),
  ModelType = character(),
  ModelLevel = character(),
  DevianceExplained = numeric(),
  stringsAsFactors = FALSE
)

#Loop of 10 replicates per species and model level
set.seed(123) # Set seed to permit reproducible results
for (i in 1:10) {
  # Random sample of 400 of model data rows
  sampled_data <- model.data[sample(nrow(model.data), size = 400), ]
  
  # UNIFORM SPECIES
  fitU.MA <- glm(Abund.unif ~ Pred.unif, data = sampled_data, family = poisson())
  fitU.MAME <- glm(Abund.unif ~ Pred.unif + DIF_C.100.10 + DIF_V.100.10,
                   data = sampled_data, family = poisson())
  fitU.MAMEMI <- glm(Abund.unif ~ Pred.unif + DIF_C.100.10 + DIF_V.100.10 +
                       VEGT_1000,
                     data = sampled_data, family = poisson())
  
  # Explained diviance calculation
  DMA <- ((summary(fitU.MA)$null.deviance - summary(fitU.MA)$deviance) / summary(fitU.MA)$null.deviance) * 100
  DMAME <- ((summary(fitU.MAME)$null.deviance - summary(fitU.MAME)$deviance) / summary(fitU.MAME)$null.deviance) * 100
  DMAMEMI <- ((summary(fitU.MAMEMI)$null.deviance - summary(fitU.MAMEMI)$deviance) / summary(fitU.MAMEMI)$null.deviance) * 100
  
  # Save results of fitU
  results <- rbind(results, data.frame(
    Replica = i,
    ModelType = "fitU",
    ModelLevel = c("MA", "MAME", "MAMEMI"),
    DevianceExplained = c(DMA, DMAME, DMAMEMI)
  ))
  
  # AGGREGATED SPECIES
  fitG.MA <- glm(Abund.aggr ~ Pred.aggr, data = sampled_data, family = poisson())
  fitG.MAME <- glm(Abund.aggr ~ Pred.aggr + DIF_C.100.10 + DIF_V.100.10,
                   data = sampled_data, family = poisson())
  fitG.MAMEMI <- glm(Abund.aggr ~ Pred.aggr + DIF_C.100.10 + DIF_V.100.10 +
                       VEGT_1000,
                     data = sampled_data, family = poisson())
  
  # Explained deviance calculation
  DMA2 <- ((summary(fitG.MA)$null.deviance - summary(fitG.MA)$deviance) / summary(fitG.MA)$null.deviance) * 100
  DMAME2 <- ((summary(fitG.MAME)$null.deviance - summary(fitG.MAME)$deviance) / summary(fitG.MAME)$null.deviance) * 100
  DMAMEMI2 <- ((summary(fitG.MAMEMI)$null.deviance - summary(fitG.MAMEMI)$deviance) / summary(fitG.MAMEMI)$null.deviance) * 100
  
  # Save results of fitG
  results <- rbind(results, data.frame(
    Replica = i,
    ModelType = "fitG",
    ModelLevel = c("MA", "MAME", "MAMEMI"),
    DevianceExplained = c(DMA2, DMAME2, DMAMEMI2)
  ))
  print(paste("loop model", i))
}

# Visualize results
print(results)

#Check and prepare data for analysis
str(results)

results$Type.f <- as.factor(results$ModelType)
results$Res.f <- as.factor(results$ModelLevel)

#Quantify and check interaction
#Complete model (with interaction)
lm_full <- lm(DevianceExplained ~ Res.f * Type.f, data = results)
summary(lm_full)
dev_full <- summary(lm_full)$r.squared

#Reduced model (without interaction)
lm_reduced <- lm(DevianceExplained ~ Res.f + Type.f, data = results)
dev_reduced <- summary(lm_reduced)$r.squared

interaction_effect <- dev_full - dev_reduced

interaction_effect * 100


# Now visualize the interaction and independdent effects
library(ggplot2)

# Crear el gráfico con barras de error
ggplot(results, aes(x = ModelLevel, y = DevianceExplained, color = ModelType, group = ModelType)) +
  # Líneas conectando los puntos
  stat_summary(fun = mean, geom = "line", size = 1) +
  # Puntos de la media
  stat_summary(fun = mean, geom = "point", size = 3) +
  # Barras de error
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +
  labs(
    title = "Interacción entre ModelLevel y ModelType",
    x = "Model Level",
    y = "Deviance Explained (%)"
  ) +
  scale_color_manual(values = c("cyan2", "coral")) + # Personalizar colores de las líneas
  theme_minimal()



# QUANTILE REGRESSIONS ---------------------------------------------------


model.data$Abund.aggr01 <- model.data$Abund.aggr/max(model.data$Abund.aggr)
model.data$Abund.unif01 <- model.data$Abund.unif/max(model.data$Abund.unif)

model.data$Abund.aggrlog <- log(model.data$Abund.aggr01 + 1)
model.data$Abund.uniflog <- log(model.data$Abund.unif01 + 1)


# Random sample of 400 model data rows
set.seed(3333)
sampled_data <- model.data[sample(nrow(model.data), size = 400), ]


library(ggplot2)
library(qgam)

#Aggregate columns of covariates that determine suitability for each model level
sampled_data$Pred_agg_1 <- sampled_data$Pred.aggr
sampled_data$Pred_agg_2 <- sampled_data$Pred.aggr + sampled_data$DIF_C.100.10 + sampled_data$DIF_V.100.10
sampled_data$Pred_agg_3 <- sampled_data$Pred.aggr + sampled_data$DIF_C.100.10 + sampled_data$DIF_V.100.10 + sampled_data$VEGT_1000


#Define percentiles
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)

#BROAD scale
#Adjust non-linear models per percentil with qgam
quantile_models <- lapply(tau, function(q) {
  qgam(Abund.aggr01 ~ s(Pred_agg_1, bs = "cs"), qu = q, data = sampled_data)
})

#Extract predictions of the models for plotting
predictions <- do.call(rbind, lapply(seq_along(tau), function(i) {
  data.frame(
    Pred.aggr = sampled_data$Pred_agg_1,
    Quantile = tau[i],
    Prediction = predict(quantile_models[[i]], newdata = sampled_data)
  )
}))

#Create the plot with the non-linear quantile regressions
ggplot(sampled_data, aes(Pred_agg_1, Abund.aggr01)) +
  geom_point(alpha = 0.5) +
  geom_line(data = predictions, aes(x = Pred.aggr, y = Prediction, group = Quantile, color = as.factor(Quantile)), size = 1) +
  scale_color_viridis_d(name = "Quantile") +
  ggtitle("Regresión cuantilica suitability_abundancia (no lineal)") +
  theme_minimal()

#BROAD+INTERMEDIATE scale
#Adjust non-linear models per percentil with qgam
quantile_models <- lapply(tau, function(q) {
  qgam(Abund.aggr01 ~ s(Pred_agg_2, bs = "cs"), qu = q, data = sampled_data)
})

#Extract predictions of the models for plotting
predictions <- do.call(rbind, lapply(seq_along(tau), function(i) {
  data.frame(
    Pred.aggr = sampled_data$Pred_agg_2,
    Quantile = tau[i],
    Prediction = predict(quantile_models[[i]], newdata = sampled_data)
  )
}))

#Create the plot with the non-linear quantile regressions
ggplot(sampled_data, aes(Pred_agg_2, Abund.aggr01)) +
  geom_point(alpha = 0.5) +
  geom_line(data = predictions, aes(x = Pred.aggr, y = Prediction, group = Quantile, color = as.factor(Quantile)), size = 1) +
  scale_color_viridis_d(name = "Quantile") +
  ggtitle("Regresión cuantilica suitability_abundancia (no lineal)") +
  theme_minimal()

#BROAD+INTERMEDIATE+LOCAL scale
#Adjust non-linear models per percentil with qgam
quantile_models <- lapply(tau, function(q) {
  qgam(Abund.aggr01 ~ s(Pred_agg_3, bs = "cs"), qu = q, data = sampled_data)
})

#Extract predictions of the models for plotting
predictions <- do.call(rbind, lapply(seq_along(tau), function(i) {
  data.frame(
    Pred.aggr = sampled_data$Pred_agg_3,
    Quantile = tau[i],
    Prediction = predict(quantile_models[[i]], newdata = sampled_data)
  )
}))

#Create the plot with the non-linear quantile regressions
ggplot(sampled_data, aes(Pred_agg_3, Abund.aggr01)) +
  geom_point(alpha = 0.5) +
  geom_line(data = predictions, aes(x = Pred.aggr, y = Prediction, group = Quantile, color = as.factor(Quantile)), size = 1) +
  scale_color_viridis_d(name = "Quantile") +
  ggtitle("Regresión cuantilica suitability_abundancia (no lineal)") +
  theme_minimal()

#UNIFORM
#Aggregate columns of covariates that determine suitability for each model level
sampled_data$Pred_uni_1 <- sampled_data$Pred.unif
sampled_data$Pred_uni_2 <- sampled_data$Pred.unif + sampled_data$DIF_C.100.10 + sampled_data$DIF_V.100.10
sampled_data$Pred_uni_3 <- sampled_data$Pred.unif + sampled_data$DIF_C.100.10 + sampled_data$DIF_V.100.10 + sampled_data$VEGT_1000


#Define percentiles
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)

#BROAD scale
#Adjust non-linear models per percentil with qgam
quantile_models <- lapply(tau, function(q) {
  qgam(Abund.unif01 ~ s(Pred_uni_1, bs = "cs"), qu = q, data = sampled_data)
})

#Extract predictions of the models for plotting
predictions <- do.call(rbind, lapply(seq_along(tau), function(i) {
  data.frame(
    Pred.aggr = sampled_data$Pred_uni_1,
    Quantile = tau[i],
    Prediction = predict(quantile_models[[i]], newdata = sampled_data)
  )
}))

#Create the plot with the non-linear quantile regressions
ggplot(sampled_data, aes(Pred_uni_1, Abund.unif01)) +
  geom_point(alpha = 0.5) +
  geom_line(data = predictions, aes(x = Pred.aggr, y = Prediction, group = Quantile, color = as.factor(Quantile)), size = 1) +
  scale_color_viridis_d(name = "Quantile") +
  ggtitle("Regresión cuantilica suitability_abundancia (no lineal)") +
  theme_minimal()

#BROAD+INTERMEDIATE scale
#Adjust non-linear models per percentil with qgam
quantile_models <- lapply(tau, function(q) {
  qgam(Abund.unif01 ~ s(Pred_uni_2, bs = "cs"), qu = q, data = sampled_data)
})

#Extract predictions of the models for plotting
predictions <- do.call(rbind, lapply(seq_along(tau), function(i) {
  data.frame(
    Pred.aggr = sampled_data$Pred_uni_2,
    Quantile = tau[i],
    Prediction = predict(quantile_models[[i]], newdata = sampled_data)
  )
}))

#Create the plot with the non-linear quantile regressions
ggplot(sampled_data, aes(Pred_uni_2, Abund.unif01)) +
  geom_point(alpha = 0.5) +
  geom_line(data = predictions, aes(x = Pred.aggr, y = Prediction, group = Quantile, color = as.factor(Quantile)), size = 1) +
  scale_color_viridis_d(name = "Quantile") +
  ggtitle("Regresión cuantilica suitability_abundancia (no lineal)") +
  theme_minimal()

#BROAD+INTERMEDIATE+LOCAL scale
#Adjust non-linear models per percentil with qgam
quantile_models <- lapply(tau, function(q) {
  qgam(Abund.unif01 ~ s(Pred_uni_3, bs = "cs"), qu = q, data = sampled_data)
})

#Extract predictions of the models for plotting
predictions <- do.call(rbind, lapply(seq_along(tau), function(i) {
  data.frame(
    Pred.aggr = sampled_data$Pred_uni_3,
    Quantile = tau[i],
    Prediction = predict(quantile_models[[i]], newdata = sampled_data)
  )
}))

#Create the plot with the non-linear quantile regressions
ggplot(sampled_data, aes(Pred_uni_3, Abund.unif01)) +
  geom_point(alpha = 0.5) +
  geom_line(data = predictions, aes(x = Pred.aggr, y = Prediction, group = Quantile, color = as.factor(Quantile)), size = 1) +
  scale_color_viridis_d(name = "Quantile") +
  ggtitle("Regresión cuantilica suitability_abundancia (no lineal)") +
  theme_minimal()



# END ---------------------------------------------------------------------


