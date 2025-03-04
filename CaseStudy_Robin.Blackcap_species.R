##################################################
#                  CASE STUDY                    #
#              ROBIN AND BLACKCAP                #
##################################################


#Determine the work directory
setwd("mypath")

library(terra)


########
#MAXENT#
########
#All data employed were on this repository
#mypath

#We employ MaxEnt software
#Version 3.4.4  
#The parameters that have been taken into account were: 
#In SAMPLES load csv file of each species (Erithacus rubecula, Sylvia atricapilla)
#In ENVIRONMENTAL LAYERS load predictor layers considered for the model (sparse vegetation, cropland, shrubland, open foret, closed forest, herbaceous, precipitation, temperature)
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
#Bias file = Human Footprint layer


# GLM models ---------------------------------------------

#Charge the dataset
Complete_data <- readRDS("review_data.Casestudy.rds")

#prepare variables
data <- data.frame(scale(Complete_data[,c(3:17)]))
data$VEGETACION <- Complete_data$VEGETACION
data$ERRUB_500M <- Complete_data$ERRUB_500M
data$SYAT_500M <- Complete_data$SYAT_500M


#Create a data frame to save all the results
results <- data.frame(
  Replica = integer(),
  ModelType = character(),
  ModelLevel = character(),
  DevianceExplained = numeric(),
  stringsAsFactors = FALSE
)

#Loop of 10 replicates per species
set.seed(123) # Set seed to permit reproducible results

for (i in 1:10) {
  # Random sample of 70% of model data rows
  sampled_data <- data[sample(nrow(data), size = 0.7 * nrow(data)), ]
  
  # BLACKCAP MODELS
  fitC.MA <- glm(SYAT_500M ~ MXCAPI, data = sampled_data, family = poisson())
  fitC.MAME <- glm(SYAT_500M ~ MXCAPI + DIFF_ABIERTO + DIFF_ARBUSTO + 
                     DIFF_BOSQUEABIE + DIFF_BOSQUECERR,
                   data = sampled_data, family = poisson())
  fitC.MAMEMI <- glm(SYAT_500M ~ MXCAPI + + DIFF_ABIERTO + DIFF_ARBUSTO + 
                       DIFF_BOSQUEABIE + DIFF_BOSQUECERR + VEGETACION,
                     data = sampled_data, family = poisson())
  
  # Explained deviance calculation
  DMA <- ((summary(fitC.MA)$null.deviance - summary(fitC.MA)$deviance) / summary(fitC.MA)$null.deviance) * 100
  DMAME <- ((summary(fitC.MAME)$null.deviance - summary(fitC.MAME)$deviance) / summary(fitC.MAME)$null.deviance) * 100
  DMAMEMI <- ((summary(fitC.MAMEMI)$null.deviance - summary(fitC.MAMEMI)$deviance) / summary(fitC.MAMEMI)$null.deviance) * 100
  
  # Save results of fitC
  results <- rbind(results, data.frame(
    Replica = i,
    ModelType = "fitC",
    ModelLevel = c("MA", "MAME", "MAMEMI"),
    DevianceExplained = c(DMA, DMAME, DMAMEMI)
  ))
  
  # ROBIN MODELS
  fit1MAMEMI<- glm(ERRUB_500M ~ MXPETI + DIFF_HERBACEO +  
                     DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_BOSQUECERR + 
                     DIFF_BOSQUEABIE + DIFF_ARBUSTO + VEGETACION, data = data, family = poisson())
  
  fitP.MA <- glm(ERRUB_500M ~ MXPETI, data = sampled_data, family = poisson())
  fitP.MAME <- glm(ERRUB_500M ~ MXPETI + DIFF_HERBACEO +  
                     DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_BOSQUECERR + 
                     DIFF_BOSQUEABIE,
                   data = sampled_data, family = poisson())
  fitP.MAMEMI <- glm(ERRUB_500M ~ MXPETI + DIFF_HERBACEO +  
                       DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_BOSQUECERR + 
                       DIFF_BOSQUEABIE + VEGETACION,
                     data = sampled_data, family = poisson())
  
  # Explained deviance calculation
  DMA2 <- ((summary(fitP.MA)$null.deviance - summary(fitP.MA)$deviance) / summary(fitP.MA)$null.deviance) * 100
  DMAME2 <- ((summary(fitP.MAME)$null.deviance - summary(fitP.MAME)$deviance) / summary(fitP.MAME)$null.deviance) * 100
  DMAMEMI2 <- ((summary(fitP.MAMEMI)$null.deviance - summary(fitP.MAMEMI)$deviance) / summary(fitP.MAMEMI)$null.deviance) * 100
  
  # Save results of fitP
  results <- rbind(results, data.frame(
    Replica = i,
    ModelType = "fitP",
    ModelLevel = c("MA", "MAME", "MAMEMI"),
    DevianceExplained = c(DMA2, DMAME2, DMAMEMI2)
  ))
  print(paste("loop model", i))
}

#Visualize results
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

#Now we visualize the interaction and independent effect
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


# Venn Diagrams models ----------------------------------------------------



#GLM MODELS (Run all the levels and combinations, to extract the deviance explained by all the levels: M1, M2 and M3)
#ROBIN#
fit1MA <- glm(ERRUB_500M ~ MXPETI, data = data, family = poisson()) #Create a model with only broad-scale covariates

fit1MI <- glm(ERRUB_500M ~ VEGETACION, data = data, family = poisson()) #Create a model with only local-scale covariates

#Use a stepwise to know which predictors at intermediate-scale level we need to consider
SELECTED<-step(glm(ERRUB_500M~ MXPETI ,poisson(link = "log"), data = data), scope= ~DIFF_ARBUSTO+DIFF_HERBACEO+DIFF_AGRICOLA+DIFF_ABIERTO+DIFF_BOSQUECERR+DIFF_BOSQUEABIE+MXPETI, direction=c("both"))

fit1ME <- glm(ERRUB_500M ~ DIFF_HERBACEO +  #Create a model with only intermediate-scale covariates
                DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_BOSQUECERR + 
                DIFF_BOSQUEABIE + DIFF_ARBUSTO, data = data, family = poisson())

fit1MAME <- glm(ERRUB_500M ~ MXPETI + DIFF_HERBACEO +  #Create a model with only broad- and intermediate-scale covariates
                  DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_BOSQUECERR + 
                  DIFF_BOSQUEABIE + DIFF_ARBUSTO, data = data, family = poisson())

fit1MAMI<- glm(ERRUB_500M ~ MXPETI + VEGETACION, data = data, family = poisson()) #Create a model with only broad- and local-scale covariates

fit1MEMI<- glm(ERRUB_500M ~ DIFF_HERBACEO + DIFF_AGRICOLA + #Create a model with only intermediate- and local-scale covariates
                 DIFF_ABIERTO + DIFF_BOSQUECERR + DIFF_BOSQUEABIE + DIFF_ARBUSTO +
                 VEGETACION, data = data, family = poisson())

fit1MAMEMI<- glm(ERRUB_500M ~ MXPETI + DIFF_HERBACEO +  #Create a model with broad-, intermediate- and local-scale covariates
                   DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_BOSQUECERR + 
                   DIFF_BOSQUEABIE + DIFF_ARBUSTO + VEGETACION, data = data, family = poisson())

#DEVIANCE CALCULATION#
MA <- summary(fit1MA)
ME <- summary(fit1ME)
MI <- summary(fit1MI)
MAMI <- summary(fit1MAMI)
MEMI <- summary(fit1MEMI)
MAME <- summary(fit1MAME)
MAMEMI <- summary(fit1MAMEMI)

#Extract the Deviance explained by each model
DMA <- ((MA$null.deviance - MA$deviance)/MA$null.deviance)*100
DME <- ((ME$null.deviance - ME$deviance)/ME$null.deviance)*100
DMI <- ((MI$null.deviance - MI$deviance)/MI$null.deviance)*100
DMAMI <- ((MAMI$null.deviance - MAMI$deviance)/MAMI$null.deviance)*100
DMEMI <- ((MEMI$null.deviance - MEMI$deviance)/MEMI$null.deviance)*100
DMAME <- ((MAME$null.deviance - MAME$deviance)/MAME$null.deviance)*100
DMAMEMI <- ((MAMEMI$null.deviance - MAMEMI$deviance)/MAMEMI$null.deviance)*100

#Calculate the deviance explained only by each level (Broad-, Intermediate-, Local-scale)
FMA <- DMAMEMI - DMEMI
FME <- DMAMEMI - DMAMI
FMI <- DMAMEMI - DMAME

#Calculate the deviance explained by combination of pair levels
IMAME <- DMAMEMI - FMA - FME - DMI
IMEMI <- DMAMEMI - FME - FMI - DMA
IMAMI <- DMAMEMI - FMA - FMI - DME

#Calculate the deviance explained by combination of all levels
IMAMEMI <- DMAMEMI - FMA - FME - FMI - IMAME - IMAMI - IMEMI

#Change % of deviance explained by each level to % of contribution of each level
DLIST <- c(FMA, FME, FMI, IMAME, IMEMI, IMAMI, IMAMEMI)
CLIST <- c()
for (i in 1:length(DLIST)) {
  value.i <- (DLIST[[i]]/DMAMEMI)*100
  CLIST <- c(CLIST, value.i)
}

#Round the decimal values
Area1 <- round(CLIST[1]+CLIST[4]+CLIST[6]+CLIST[7], 2)
Area2 <- round(CLIST[2]+CLIST[4]+CLIST[5]+CLIST[7], 2)
Area3 <- round(CLIST[3]+CLIST[5]+CLIST[6]+CLIST[7], 2)
Area12 <- round(CLIST[4]+CLIST[7], 2)
Area23 <- round(CLIST[5]+CLIST[7], 2)
Area13 <- round(CLIST[6]+CLIST[7], 2)
Area123 <- round(CLIST[7], 2)

#Create a Venn diagram of these values
library(VennDiagram)

# move to new plotting page
grid.newpage()


# create Venn diagram with three sets
overrideTriple=TRUE
draw.triple.venn(area1=Area1, 
                 area2=Area2, 
                 area3=Area3, 
                 n12=Area12, 
                 n23=Area23, 
                 n13=Area13, 
                 n123=Area123, 
                 category=c("Broad","Intermediate","Local"),
                 rotation = 1, reverse = FALSE, col="Black",
                 fill=c("yellow", "blue","green"), 
                 lty = "blank")

#BLACKCAP#

fit2MA <- glm(SYAT_500M ~ MXCAPI, data = data, family = poisson()) #Create a model with only broad-scale covariates

fit2MI <- glm(SYAT_500M ~ VEGETACION, data = data, family = poisson()) #Create a model with only local-scale covariates

#Use a stepwise to know which predictors at intermediate-scale level we need to consider
SELECTED<-step(glm(SYAT_500M~ MXCAPI ,poisson(link = "log"), data = data), scope= ~DIFF_ARBUSTO+DIFF_HERBACEO+DIFF_AGRICOLA+DIFF_ABIERTO+DIFF_BOSQUECERR+DIFF_BOSQUEABIE+MXCAPI, direction=c("both"))

fit2ME <- glm(SYAT_500M ~ DIFF_HERBACEO +  #Create a model with only intermediate-scale covariates
                DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_ARBUSTO + 
                DIFF_BOSQUEABIE + DIFF_BOSQUECERR, data = data, family = poisson())

fit2MAME <- glm(SYAT_500M ~ MXCAPI + DIFF_HERBACEO +  #Create a model with only broad- and intermediate-scale covariates
                  DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_ARBUSTO + 
                  DIFF_BOSQUEABIE + DIFF_BOSQUECERR, data = data, family = poisson())

fit2MAMI<- glm(SYAT_500M ~ MXCAPI + VEGETACION, data = data, family = poisson()) #Create a model with only broad- and local-scale covariates

fit2MEMI<- glm(SYAT_500M ~ DIFF_HERBACEO + DIFF_AGRICOLA + #Create a model with only intermediate- and local-scale covariates
                 DIFF_ABIERTO + DIFF_ARBUSTO + DIFF_BOSQUEABIE + DIFF_BOSQUECERR +
                 VEGETACION, data = data, family = poisson())

fit2MAMEMI<- glm(SYAT_500M ~ MXCAPI + DIFF_HERBACEO +  #Create a model with broad-, intermediate and local-scale covariates
                   DIFF_AGRICOLA + DIFF_ABIERTO + DIFF_ARBUSTO + 
                   DIFF_BOSQUEABIE + DIFF_BOSQUECERR + VEGETACION, data = data, family = poisson())

#DEVIANCE CALCULATION#
MA <- summary(fit2MA)
ME <- summary(fit2ME)
MI <- summary(fit2MI)
MAMI <- summary(fit2MAMI)
MEMI <- summary(fit2MEMI)
MAME <- summary(fit2MAME)
MAMEMI <- summary(fit2MAMEMI)

#Extract the Deviance explained by each model
DMA <- ((MA$null.deviance - MA$deviance)/MA$null.deviance)*100
DME <- ((ME$null.deviance - ME$deviance)/ME$null.deviance)*100
DMI <- ((MI$null.deviance - MI$deviance)/MI$null.deviance)*100
DMAMI <- ((MAMI$null.deviance - MAMI$deviance)/MAMI$null.deviance)*100
DMEMI <- ((MEMI$null.deviance - MEMI$deviance)/MEMI$null.deviance)*100
DMAME <- ((MAME$null.deviance - MAME$deviance)/MAME$null.deviance)*100
DMAMEMI <- ((MAMEMI$null.deviance - MAMEMI$deviance)/MAMEMI$null.deviance)*100

#Calculate the deviance explained only by each level (Broad-, intermediate-, local-scale)
FMA <- DMAMEMI - DMEMI
FME <- DMAMEMI - DMAMI
FMI <- DMAMEMI - DMAME

#Calculate the deviance explained by combination of pair levels
IMAME <- DMAMEMI - FMA - FME - DMI
IMEMI <- DMAMEMI - FME - FMI - DMA
IMAMI <- DMAMEMI - FMA - FMI - DME

#Calculate the deviance explained by combination of all levels
IMAMEMI <- DMAMEMI - FMA - FME - FMI - IMAME - IMAMI - IMEMI

#Change % of deviance explained by each level to % of contribution of each level
DLIST <- c(FMA, FME, FMI, IMAME, IMEMI, IMAMI, IMAMEMI)
CLIST <- c()
for (i in 1:length(DLIST)) {
  value.i <- (DLIST[[i]]/DMAMEMI)*100
  CLIST <- c(CLIST, value.i)
}

#Round the decimal values
Area1 <- round(CLIST[1]+CLIST[4]+CLIST[6]+CLIST[7], 2)
Area2 <- round(CLIST[2]+CLIST[4]+CLIST[5]+CLIST[7], 2)
Area3 <- round(CLIST[3]+CLIST[5]+CLIST[6]+CLIST[7], 2)
Area12 <- round(CLIST[4]+CLIST[7], 2)
Area23 <- round(CLIST[5]+CLIST[7], 2)
Area13 <- round(CLIST[6]+CLIST[7], 2)
Area123 <- round(CLIST[7], 2)

# move to new plotting page
grid.newpage()


# create Venn diagram with these values
overrideTriple=TRUE
draw.triple.venn(area1=Area1, 
                 area2=Area2, 
                 area3=Area3, 
                 n12=Area12, 
                 n23=Area23, 
                 n13=Area13, 
                 n123=Area123, 
                 category=c("Broad","Intermediate","Local"),
                 rotation = 1, reverse = FALSE, col="Black",
                 fill=c("yellow","blue","green"), 
                 lty = "blank")



