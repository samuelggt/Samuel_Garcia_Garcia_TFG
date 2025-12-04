# Clear workspace
rm(list = ls())
# Save original parameters
op <- par()

# Instalar paquetes necesarios

if(!require("glmnet")) install.packages("glmnet"); library(glmnet)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("pROC")) install.packages("pROC"); library(pROC)
if(!require("caret")) install.packages("caret"); library(caret)
if(!require("themis")) install.packages("themis"); library(themis)
if(!require("ROSE")) install.packages("ROSE"); library(ROSE)

#-------Cargar Datos

datos <- read.csv("Data/STC_data.csv")

#-------Crear variable binaria

datos$BinDx <- NA
datos$BinDx[datos$Dx<3] <- 0
datos$BinDx[datos$Dx==3] <- 1

# Escribirla como factor (enet lo requiere)
  datos$BinDx <- factor(datos$BinDx) 

#-------Dataframes Variables ultrasonido

# Dataframe con todas las variables 
  datos0v3_av <- datos[, c(names(datos)[19:106], "BinDx")]
  datos0v3_av <- na.omit(datos0v3_av[which(datos$BinDx %in% "0" | datos$BinDx %in% "1"),])

# Dataframe con las variables de imagenes transversales
  datos0v3_tv <- datos[, c(names(datos)[19:62], "BinDx")]
  datos0v3_tv <- na.omit(datos0v3_tv[which(datos$BinDx %in% "0" | datos$BinDx %in% "1"),])

# Dataframe con las variables de imagenes longitudinales
  datos0v3_lv <- datos[, c(names(datos)[63:106], "BinDx")]
  datos0v3_lv <- na.omit(datos0v3_lv[which(datos$BinDx %in% "0" | datos$BinDx %in% "1"),])

#-------Definición matrices para cada conjunto de datos

X_av <- model.matrix(BinDx ~. , data=datos0v3_av)[,-1]
X_tv <- model.matrix(BinDx ~. , data=datos0v3_tv)[,-1]
X_lv <- model.matrix(BinDx ~. , data=datos0v3_lv)[,-1]

y_av <- datos0v3_av$BinDx
y_tv <- datos0v3_tv$BinDx
y_lv <- datos0v3_lv$BinDx

# Vector de los métodos de muestreo
  smpMet <- c("none", "down", "up", "smote", "rose")

#-------Función red elástica para todos los métodos de muestreo que devuelve los AUC
enet_loocv <- function(tag, datos, X, y) {
  
  # Definir vector AUCs
    aucs <- numeric(length(smpMet))
    names(aucs) <- smpMet
  
  # Bucle para los métodos de muestreo
    for (mm in 1:length(smpMet)) {
      
      #Mostrar progreso
        cat(paste0("\nProcesando método: ", smpMet[mm], "\n"))
      
      # Definir vector predicciones
        preds <- numeric(nrow(datos))
      
      # Model inputs
      # Sampling = "none" (traincontrol da error, se debe poner por separado)
        if(mm == 1){
          fitControl <- trainControl(method = "repeatedcv",
                                     number = 5,
                                     repeats = 5)
        }
      # Resto de metodos de muestreo
        else{
        fitControl <- trainControl(method = "repeatedcv",
                                   number = 5,
                                   repeats = 5,
                                   sampling = smpMet[mm])
        }
        
        #Bucle loocv
        for (ss in 1:nrow(datos)) {
          
          # Mostrar progreso
            cat(paste0("..", ss))
          
          # Train-test
            trainSet <- datos[-ss,]
            testSet <- datos[ss,]
          
          # Modelos entrenamiento  
          model_tr <- train(BinDx ~ ., data = trainSet,
                              method = "glmnet",
                              metric = "Kappa",
                              trControl = fitControl,
                              family = "binomial")
          
          # Train-test
            xtrain <- X[-ss,]
            ytrain <- y[-ss]
            xtest <- X[ss,]
          
          # Diseño mejor modelo
            enet_best <- glmnet(x = xtrain,
                                y = ytrain,
                                alpha = model_tr$bestTune$alpha,
                                lambda = model_tr$bestTune$lambda,
                                family = binomial(link = "logit"))
          
          # Guardar valor predicciones
            preds[ss] <- predict(enet_best, s = model_tr$bestTune$lambda, newx = xtest, type = "response")
            
        } # Fin de bucle LOOCV
      
      # Obtener y guardar AUC en la coordenada del vector
        aucs[mm] <- roc(datos$BinDx, preds)$auc
        
    } # Fin de bucle métodos de muestreo
  
  # Devolver vector AUCS
    return(aucs)
}

#-------Ejecutar la función y guardar RDS con el vector de aucs

saveRDS(enet_loocv(av, datos0v3_av, X_av, y_av), 
        "Matrices/enet_aucs_av.rds")

saveRDS(enet_loocv(tv, datos0v3_tv, X_tv, y_tv), 
        "Matrices/enet_aucs_tv.rds")

saveRDS(enet_loocv(lv, datos0v3_lv, X_lv, y_lv), 
        "Matrices/enet_aucs_lv.rds")

#-------Cargar RDS

aucs_av <- readRDS("Matrices/enet_aucs_av.rds")
aucs_tv <- readRDS("Matrices/enet_aucs_tv.rds")
aucs_lv <- readRDS("Matrices/enet_aucs_lv.rds")

#-------Mostrar AUCS

# Modelo completo
  aucs_av
  
# Modelo transversal
  aucs_tv
  
# Modelo longitudinal
  aucs_lv

