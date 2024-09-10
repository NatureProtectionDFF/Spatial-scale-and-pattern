##########################################
#               SIMULATIONS              #
#        UNIFORM/AGGREGATED SPECIES      #
#        DAVID FERRER  19/11/2023        #
#       VERSION 1000 X 1000 SQUARES      #
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
  route_rst<-paste("capas_generalista.especialista/VEGT/vegt",NUM,".asc", sep="")
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
  route_rst<-paste("capas_generalista.especialista/CLIM/clim",NUM,".asc", sep="")
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
setwd("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones")
route_rst<-paste("capas_generalista.especialista/especialista/ADQ.asc", sep="")
writeRaster(ADQ, route_rst, NAflag = -9999, overwrite= TRUE)

route_rst<-paste("capas_generalista.especialista/especialista/error_nodist.asc", sep="")
writeRaster(r_error, route_rst, NAflag = -9999, overwrite= TRUE)

route_rst<-paste("capas_generalista.especialista/especialista/ADQ_lgt_R_nodist.asc", sep="")
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
ADQ_lgt_R <- rast("capas_generalista.especialista/especialista/ADQ_lgt_R_nodist.asc")
ADQ_lgt_R <- rast("capas_generalista.especialista/generalista/ADQ_lgt_R_nodist.asc")


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


setwd("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones/capas_generalista.especialista")
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
setwd("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones")
lista_rasters<-list()
for (j in 1:length(list.files("capas_generalista.especialista/VEGT/", pattern = "^vegt.*\\.asc$"))) {
  lista_rasters[[j]]<-rast(paste("capas_generalista.especialista/VEGT/", list.files("capas_generalista.especialista/VEGT/", pattern = "^vegt.*\\.asc$")[j], sep=""))
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
base <- "D:/david.ferrer"
savepath <- file.path(base, "SIMULACIONES/predictores/especialista.generalista")
if (!dir.exists(savepath)) {
  dir.create(savepath, recursive=TRUE)
  print(paste("Creating directory", savepath))
} else {
  print(paste("Directory", savepath, "already exists. Overwriting."))
}

for(i in 1:length(ls(pattern = "VEGT_"))) {
  terra::writeRaster(get(ls(pattern = "VEGT_")[i]), 
                     paste0("D:/david.ferrer/SIMULACIONES/predictores/especialista.generalista/VEGT", 
                            gsub("VEGT_", "", ls(pattern = "VEGT_")[i]),
                            ".asc"), NAflag = -9999, overwrite = TRUE)
}
rm(i)

#Charge CLIME layers
setwd("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones")
lista_rasters<-list()
for (j in 1:length(list.files("capas_generalista.especialista/CLIM/", pattern = "^clim.*\\.asc$"))) {
  lista_rasters[[j]]<-rast(paste("capas_generalista.especialista/CLIM/", list.files("capas_generalista.especialista/CLIM/", pattern = "^clim.*\\.asc$")[j], sep=""))
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
base <- "D:/david.ferrer"
savepath <- file.path(base, "SIMULACIONES/predictores/especialista.generalista")
if (!dir.exists(savepath)) {
  dir.create(savepath, recursive=TRUE)
  print(paste("Creating directory", savepath))
} else {
  print(paste("Directory", savepath, "already exists. Overwriting."))
}

for(i in 1:length(ls(pattern = "CLIM_"))) {
  terra::writeRaster(get(ls(pattern = "CLIM_")[i]), 
                     paste0("D:/david.ferrer/SIMULACIONES/predictores/especialista.generalista/CLIM", 
                            gsub("CLIM_", "", ls(pattern = "CLIM_")[i]),
                            ".asc"), NAflag = -9999, overwrite = TRUE)
}
rm(i)



# MAXENT SETTINGS ---------------------------------------------------------


#All data employed were on this repository
#D:\david.ferrer\SIMULACIONES

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



# SUITABILITY-ABUNDANCE ANALYSIS ------------------------------------------


#(DO IT FOR UNIFORM AND AGGREGATED SPATIAL PATTTERN SPECIES)
setwd("D:/david.ferrer/SIMULACIONES")

## UNIFORM PATTERN SPECIES -------------------------------------------------

### Resolution 100x100 cell size --------------------------------------------

#Charge the prediction obtained with Maxent
pred <- rast("resultados/especialista.generalista_NODIST/especialista/10/Virtual_Species_avg.asc")
r_base.1000 <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1000, ymin=0, ymax=1000, )
pred@cpp[["extent"]] <- r_base.1000@cpp[["extent"]]

data_oc <- read.csv("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones/capas_generalista.especialista/especialista/sampling_units.csv")
oc_sv <- vect(data_oc, geom=c("x", "y"))

plot(pred)
points(oc_sv, pch=16)

preds <- extract(pred, oc_sv, method="simple")

comp_data <- data.frame(cbind(data_oc, preds))

names(comp_data) <- c("x", "y", "SampleID", "Abundance", "Occurrence",
                      "ID",  "Suitability_predicted_by_Maxent")

comp_data$Abundance <- comp_data$Abundance/max(comp_data$Abundance)

library(ggplot2)
ggplot(data = comp_data, 
       aes(x = Suitability_predicted_by_Maxent, y = Abundance)) + 
  geom_point() +
  ggtitle("Abundance-Suitability relationship")

tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)

ggplot(comp_data, aes(Suitability_predicted_by_Maxent,Abundance)) + 
  geom_point() + 
  geom_quantile(quantiles = tau) +
  geom_smooth(method = 'glm', col="red", size = 2) +
  ggtitle("Regresion cuantilica suitability_abundancia Simulacion 10S")

#Now we do a GLM to obtain the linear adjustment of the relationship.
fit10S <- glm(Abundance ~ Suitability_predicted_by_Maxent, data = comp_data, family = poisson())
summary(fit10S)

#Obtain % explicated Devianze
with(summary(fit10S), 1 - deviance/null.deviance)*100

#Now we do the Quantile regresion
library(quantreg)
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)
fitQ10S <- rq(Abundance ~ Suitability_predicted_by_Maxent, tau, data = comp_data) 
plot(fitQ10S)


### Resolution 10x10 cell size ----------------------------------------------

#Charge the prediction obtained with Maxent
pred <- rast("resultados/especialista.generalista_NODIST/especialista/100/Virtual_Species_avg.asc")
r_base.1000 <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1000, ymin=0, ymax=1000, )
pred@cpp[["extent"]] <- r_base.1000@cpp[["extent"]]

data_oc <- read.csv("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones/capas_generalista.especialista/especialista/sampling_units.csv")
oc_sv <- vect(data_oc, geom=c("x", "y"))

plot(pred)
points(oc_sv, pch=16)

preds <- extract(pred, oc_sv, method="simple")

comp_data <- data.frame(cbind(data_oc, preds))

names(comp_data) <- c("x", "y", "SampleID", "Abundance", "Occurrence",
                      "ID",  "Suitability_predicted_by_Maxent")

comp_data$Abundance <- comp_data$Abundance/max(comp_data$Abundance)

library(ggplot2)
ggplot(data = comp_data, 
       aes(x = Suitability_predicted_by_Maxent, y = Abundance)) + 
  geom_point() +
  ggtitle("Abundance-Suitability relationship")

tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)

ggplot(comp_data, aes(Suitability_predicted_by_Maxent,Abundance)) + 
  geom_point() + 
  geom_quantile(quantiles = tau) +
  geom_smooth(method = 'glm', col="red", size = 2) +
  ggtitle("Regresion cuantilica suitability_abundancia Simulacion 100S")

#Now we do a GLM toobtain the linear adjustment of the relationship.
fit100S <- glm(Abundance ~ Suitability_predicted_by_Maxent, data = comp_data, family = poisson())
summary(fit100S)

#Obtain % explicated Devianze
with(summary(fit100S), 1 - deviance/null.deviance)*100

#Now we do the Quantile regresion
library(quantreg)
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)
fitQ100S <- rq(Abundance ~ Suitability_predicted_by_Maxent, tau, data = comp_data) 
plot(fitQ100S)


### Resolution 1x1 cell size ------------------------------------------------

#Charge the prediction obtained with Maxent
pred <- rast("resultados/especialista.generalista_NODIST/especialista/1000/Virtual_Species_avg.asc")
r_base.1000 <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1000, ymin=0, ymax=1000, )
pred@cpp[["extent"]] <- r_base.1000@cpp[["extent"]]

data_oc <- read.csv("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones/capas_generalista.especialista/especialista/sampling_units.csv")
oc_sv <- vect(data_oc, geom=c("x", "y"))

plot(pred)
points(oc_sv, pch=16)

preds <- extract(pred, oc_sv, method="simple")

comp_data <- data.frame(cbind(data_oc, preds))

names(comp_data) <- c("x", "y", "SampleID", "Abundance", "Occurrence",
                      "ID",  "Suitability_predicted_by_Maxent")

comp_data$Abundance <- comp_data$Abundance/max(comp_data$Abundance)

library(ggplot2)
ggplot(data = comp_data, 
       aes(x = Suitability_predicted_by_Maxent, y = Abundance)) + 
  geom_point() +
  ggtitle("Abundance-Suitability relationship")

tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)

ggplot(comp_data, aes(Suitability_predicted_by_Maxent,Abundance)) + 
  geom_point() + 
  geom_quantile(quantiles = tau) +
  geom_smooth(method = 'glm', col="red", size = 2) +
  ggtitle("Regresion cuantilica suitability_abundancia Simulacion 1000S")

#Now we do a GLM toobtain the linear adjustment of the relationship.
fit1000S <- glm(Abundance ~ Suitability_predicted_by_Maxent, data = comp_data, family = poisson())
summary(fit1000S)

#Obtain % explicated Devianze
with(summary(fit1000S), 1 - deviance/null.deviance)*100

#Now we do the Quantile regresion
library(quantreg)
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)
fitQ1000S <- rq(Abundance ~ Suitability_predicted_by_Maxent, tau, data = comp_data) 
plot(fitQ1000S)


## AGGREGATED PATTERN SPECIES ----------------------------------------------

### Resolution 100x100 cell size --------------------------------------------

#Charge the prediction obtained with Maxent
pred <- rast("resultados/especialista.generalista_NODIST/generalista/10/Virtual_Species_avg.asc")
r_base.1000 <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1000, ymin=0, ymax=1000, )
pred@cpp[["extent"]] <- r_base.1000@cpp[["extent"]]

data_oc <- read.csv("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones/capas_generalista.especialista/generalista/sampling_units.csv")
oc_sv <- vect(data_oc, geom=c("x", "y"))

plot(pred)
points(oc_sv, pch=16)

preds <- extract(pred, oc_sv, method="simple")

comp_data <- data.frame(cbind(data_oc, preds))

names(comp_data) <- c("x", "y", "SampleID", "Abundance", "Occurrence",
                      "ID",  "Suitability_predicted_by_Maxent")

comp_data$Abundance <- comp_data$Abundance/max(comp_data$Abundance)

library(ggplot2)
ggplot(data = comp_data, 
       aes(x = Suitability_predicted_by_Maxent, y = Abundance)) + 
  geom_point() +
  ggtitle("Abundance-Suitability relationship")

tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)


ggplot(comp_data, aes(Suitability_predicted_by_Maxent,Abundance)) + 
  geom_point() + 
  geom_quantile(quantiles = tau) +
  geom_smooth(method = 'glm', col="red", size = 2) +
  ggtitle("Regresion cuantilica suitability_abundancia Simulacion 10G")

#Now we do a GLM toobtain the linear adjustment of the relationship.
fit10G <- glm(Abundance ~ Suitability_predicted_by_Maxent, data = comp_data, family = poisson())
summary(fit10G)

#Obtain % explicated Devianze
with(summary(fit10G), 1 - deviance/null.deviance)*100

#Now we do the Quantile regresion
library(quantreg)
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)
fitQ10G <- rq(Abundance ~ Suitability_predicted_by_Maxent, tau, data = comp_data) 
plot(fitQ10G)

### Resolution 10x10 cell size --------------------------------------------

#Charge the prediction obtained with Maxent
pred <- rast("resultados/especialista.generalista_NODIST/generalista/100/Virtual_Species_avg.asc")
r_base.1000 <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1000, ymin=0, ymax=1000, )
pred@cpp[["extent"]] <- r_base.1000@cpp[["extent"]]

data_oc <- read.csv("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones/capas_generalista.especialista/generalista/sampling_units.csv")
oc_sv <- vect(data_oc, geom=c("x", "y"))

plot(pred)
points(oc_sv, pch=16)

preds <- extract(pred, oc_sv, method="simple")

comp_data <- data.frame(cbind(data_oc, preds))

names(comp_data) <- c("x", "y", "SampleID", "Abundance", "Occurrence",
                      "ID",  "Suitability_predicted_by_Maxent")

comp_data$Abundance <- comp_data$Abundance/max(comp_data$Abundance)

library(ggplot2)
ggplot(data = comp_data, 
       aes(x = Suitability_predicted_by_Maxent, y = Abundance)) + 
  geom_point() +
  ggtitle("Abundance-Suitability relationship")

tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)

ggplot(comp_data, aes(Suitability_predicted_by_Maxent,Abundance)) + 
  geom_point() + 
  geom_quantile(quantiles = tau) +
  geom_smooth(method = 'glm', col="red", size = 2) +
  ggtitle("Regresion cuantilica suitability_abundancia Simulacion 100G")

#Now we do a GLM to obtain the linear adjustment of the relationship.
fit100G <- glm(Abundance ~ Suitability_predicted_by_Maxent, data = comp_data, family = poisson())
summary(fit100G)

#Obtain % explicated Devianze
with(summary(fit100G), 1 - deviance/null.deviance)*100

#Now we do the Quantile regresion
library(quantreg)
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)
fitQ100G <- rq(Abundance ~ Suitability_predicted_by_Maxent, tau, data = comp_data) 
plot(fitQ100G)

### Resolution 1x1 cell size --------------------------------------------

#Charge the prediction obtained with Maxent
pred <- rast("resultados/especialista.generalista_NODIST/generalista/1000/Virtual_Species_avg.asc")
r_base.1000 <- rast(ncol=1000, nrow=1000, xmin=0, xmax=1000, ymin=0, ymax=1000, )
pred@cpp[["extent"]] <- r_base.1000@cpp[["extent"]]

data_oc <- read.csv("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/ESTANCIAS/estancia UCLM/simulaciones/capas_generalista.especialista/generalista/sampling_units.csv")
oc_sv <- vect(data_oc, geom=c("x", "y"))

plot(pred)
points(oc_sv, pch=16)

preds <- extract(pred, oc_sv, method="simple")

comp_data <- data.frame(cbind(data_oc, preds))

names(comp_data) <- c("x", "y", "SampleID", "Abundance", "Occurrence",
                      "ID",  "Suitability_predicted_by_Maxent")

comp_data$Abundance <- comp_data$Abundance/max(comp_data$Abundance)

library(ggplot2)
ggplot(data = comp_data, 
       aes(x = Suitability_predicted_by_Maxent, y = Abundance)) + 
  geom_point() +
  ggtitle("Abundance-Suitability relationship")

tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)

ggplot(comp_data, aes(Suitability_predicted_by_Maxent,Abundance)) + 
  geom_point() + 
  geom_quantile(quantiles = tau) +
  geom_smooth(method = 'glm', col="red", size = 2) +
  ggtitle("Regresion cuantilica suitability_abundancia Simulacion 1000G")

#Now we do a GLM to obtain the linear adjustment of the relationship.
fit1000G <- glm(Abundance ~ Suitability_predicted_by_Maxent, data = comp_data, family = poisson())
summary(fit1000G)

#Obtain % explicated Devianze
with(summary(fit1000G), 1 - deviance/null.deviance)*100

#Now we do the Quantile regresion
library(quantreg)
tau <- c(0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 0.99)
fitQ1000G <- rq(Abundance ~ Suitability_predicted_by_Maxent, tau, data = comp_data) 
plot(fitQ1000G)



# END ---------------------------------------------------------------------


