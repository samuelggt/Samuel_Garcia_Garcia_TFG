# Clear workspace
rm(list = ls())

# Save original parameters
op <- par()

# Instalar paquetes necesarios

if(!require("glmnet")) install.packages("glmnet"); library(glmnet)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("pROC")) install.packages("pROC"); library(pROC)

#-------Cargar Datos

datos <- read.csv("Data/STC_data.csv")

#-------Crear variable binaria

datos$BinDx <- NA
datos$BinDx[datos$Dx<3] <- 0
datos$BinDx[datos$Dx==3] <- 1

#-------Dataframes Variables ultrasonido

  # Todas las variables ultrasonido 
    datos0v3_av <- datos[, c(names(datos)[19:106], "BinDx")]
    datos0v3_av <- na.omit(datos0v3_av[which(datos$BinDx %in% "0" | datos$BinDx %in% "1"),])
    
  # Variables de imagenes transversales
    datos0v3_tv <- datos[, c(names(datos)[19:62], "BinDx")]
    datos0v3_tv <- na.omit(datos0v3_tv[which(datos$BinDx %in% "0" | datos$BinDx %in% "1"),])
    
  # Variables de imagenes longitudinales
    datos0v3_lv <- datos[, c(names(datos)[63:106], "BinDx")]
    datos0v3_lv <- na.omit(datos0v3_lv[which(datos$BinDx %in% "0" | datos$BinDx %in% "1"),])

  
#-------Definición matrices para cada conjunto de datos

# Secuencia lambdas para Lasso
  lambda_seq <- 10^seq(2, -2, by = -.1) 
    
# Todas las variables
  X_av <- model.matrix(BinDx ~. , data=datos0v3_av)[,-1]    
  y_av <- datos0v3_av$BinDx   
  lasso_coef_matrix_av <- matrix(0, nrow = nrow(datos0v3_av), ncol = ncol(X_av) + 1)
  colnames(lasso_coef_matrix_av) <- rownames(coef(cv.glmnet(X_av, y_av, alpha = 1, 
                                                            lambda = lambda_seq)))
  zvalues_matrix_av <- lasso_coef_matrix_av
  colnames(zvalues_matrix_av) <- colnames(lasso_coef_matrix_av)
  preds_av <- vector("numeric",length = nrow(datos0v3_av))

# Variables transversales
  X_tv <- model.matrix(BinDx ~. , data=datos0v3_tv)[,-1]
  y_tv <- datos0v3_tv$BinDx
  lasso_coef_matrix_tv <- matrix(0, nrow = nrow(datos0v3_tv), ncol = ncol(X_tv) + 1)
  colnames(lasso_coef_matrix_tv) <- rownames(coef(cv.glmnet(X_tv, y_tv, alpha = 1, 
                                                            lambda = lambda_seq)))
  zvalues_matrix_tv <- lasso_coef_matrix_tv
  colnames(zvalues_matrix_tv) <- colnames(lasso_coef_matrix_tv)
  preds_tv <- vector("numeric",length = nrow(datos0v3_tv))

# Variables longitudinales
  X_lv <- model.matrix(BinDx ~. , data=datos0v3_lv)[,-1]
  y_lv <- datos0v3_lv$BinDx
  lasso_coef_matrix_lv <- matrix(0, nrow = nrow(datos0v3_lv), ncol = ncol(X_lv) + 1)
  colnames(lasso_coef_matrix_lv) <- rownames(coef(cv.glmnet(X_lv, y_lv, alpha = 1, 
                                                         lambda = lambda_seq)))
  zvalues_matrix_lv <- lasso_coef_matrix_lv
  colnames(zvalues_matrix_lv) <- colnames(lasso_coef_matrix_lv)
  preds_lv <- vector("numeric",length = nrow(datos0v3_lv))


#-------Función Lasso leave-one-out que devuelve las matrices coefs, zvalues y preds

loocv_lasso <- function(datos0v3, X, y,
                            lasso_coef_matrix, zvalues_matrix, preds) {
  
  for(ss in 1:nrow(datos0v3)){
    
    # Mostrar progreso
     cat(paste0("..",ss))
    
    # Train n test
     xtrain <- X[-ss,]
     ytrain <- y[-ss]
     xtest <- X[ss,]
    
    # Diseño del modelo Lasso
     cv_output <- cv.glmnet(xtrain,
                           ytrain,
                           alpha = 1,
                           lambda = lambda_seq, 
                           nfolds = 5,
                           family = binomial(link = "logit"))
     best_lam <- cv_output$lambda.min
     lasso_best <- glmnet(xtrain,
                         ytrain,
                         alpha = 1,
                         lambda = best_lam,
                         family = binomial(link = "logit"))
    
    # Guardar coeficientes
     coefs <- as.numeric(coef(lasso_best, s = best_lam))
     lasso_coef_matrix[ss, ] <- coefs
    
    # Seleccionar las variables cuyos coeficientes son distintos de 0
     cv_vars <- colnames(lasso_coef_matrix)[(which(lasso_coef_matrix[ss,] != 0)[-1])]
    
    # Diseño modelo de regresión logística con las variables seleccionadas
      if (length(cv_vars) != 0){
      
        cv_formula <- as.formula(paste("ytrain ~", paste(cv_vars, collapse=" + ")))
        fit <- glm(cv_formula, data=datos0v3[-ss,], family = binomial(link = "logit"))
      
        #Guardar z valores
        zvalues <- summary(fit)$coefficients[,"z value"]
        cv_vars <- names(zvalues[!is.na(zvalues)])
        zvalues_matrix[ss, cv_vars] <- zvalues
      }
    
    
    # Predicción
     preds[ss] <- predict(lasso_best, s = best_lam, newx = xtest, type = "response")
    
  }
  
  return(list(coefs = lasso_coef_matrix,
              zvalues = zvalues_matrix,
              preds = preds))
}

#-------Ejecutar el lasso y guardar RDS con la lista de matrices

saveRDS(loocv_lasso(datos0v3_av, X_av, y_av,
                    lasso_coef_matrix_av, zvalues_matrix_av, preds_av), 
        "Matrices/ultrasonido_lasso_matrixes_av.rds")

saveRDS(loocv_lasso(datos0v3_tv, X_tv, y_tv,
                    lasso_coef_matrix_tv, zvalues_matrix_tv, preds_tv),
        "Matrices/ultrasonido_lasso_matrixes_tv.rds")

saveRDS(loocv_lasso(datos0v3_lv, X_lv, y_lv,
                    lasso_coef_matrix_lv, zvalues_matrix_lv, preds_lv), 
        "Matrices/ultrasonido_lasso_matrixes_lv.rds")

#-------Cargar RDS

matrices_av <- readRDS("Matrices/ultrasonido_lasso_matrixes_av.rds")
matrices_tv <- readRDS("Matrices/ultrasonido_lasso_matrixes_tv.rds")
matrices_lv <- readRDS("Matrices/ultrasonido_lasso_matrixes_lv.rds")

#-------Función Análisis modelo (Curva ROC, valores z, proporción)

analisis_modelo_lasso <- function(datos, matrices) {
  # Curva ROC
    rocobj <- roc(datos$BinDx, matrices$preds)
    print(ggroc(rocobj) + 
          coord_fixed() + 
          labs(title = "Curva ROC", x = "Especificidad", y = "Sensibilidad"))
    cat("AUC:", rocobj$auc, "\n\n")
  
  
  # Filtrar filas con valores z < 10
    cat("Dimensión matriz zvalues:", dim(matrices$zvalues), "\n") #Comprobar dimension previa
    filas_validas <- apply(abs(matrices$zvalues), 1, function(x) all(x < 10))
    matrices$zvalues <- matrices$zvalues[filas_validas, ]
    cat("Dimensión matriz zvalues tras quitar los >10:", dim(matrices$zvalues), "\n\n") #Dimension filtrada
  
  # Media de valores z por variable
    z_col_means <- apply(matrices$zvalues, 2, function(col) {
    mean(col)
    })
    cat("Medias de valores z distintas de cero:\n")
    print(z_col_means[z_col_means != 0])
  
  # Porcentaje de aparición de coeficientes distintos de 0
    coef100 <- apply(matrices$coefs, 2, function(x) {
    100 * sum(x != 0) / length(x)
    })
    cat("\nVariables con coeficientes no cero en más del 70% de iteraciones:\n")
    print(coef100[coef100 > 70])
  
  return(list(
    rocobj = rocobj,
    z_col_means = z_col_means,
    coef100 = coef100
  ))
}

#-------Ejecutar análisis de los modelos Lasso

# Modelo completo
  analisis_lasso_av <- analisis_modelo_lasso(datos0v3_av, matrices_av) 

# Modelo transversal
  analisis_lasso_tv <-analisis_modelo_lasso(datos0v3_tv, matrices_tv) 

# Modelo longitudinal
  analisis_lasso_lv <- analisis_modelo_lasso(datos0v3_lv, matrices_lv)
  
#-------Selección de variables >70% y |zvalue| > 0.5
  
cv_vars_av <- names(analisis_lasso_av$coef100)[(analisis_lasso_av$coef100 > 70) & 
                                                  (abs(analisis_lasso_av$z_col_means) > 0.5)][-1]
cv_vars_tv <- names(analisis_lasso_tv$coef100)[(analisis_lasso_tv$coef100 > 70) & 
                                                  (abs(analisis_lasso_tv$z_col_means) > 0.5)][-1]
cv_vars_lv <- names(analisis_lasso_lv$coef100)[(analisis_lasso_lv$coef100 > 70) & 
                                                  (abs(analisis_lasso_lv$z_col_means) > 0.5)][-1]

#-------Función modelo regresión logística y análisis (Curva ROC, coeficiente, valores z y p)

analisis_logit <- function(datos, variables_seleccionadas) {
  
  # Crear fórmula para la regresión logística
    cv_formula <- as.formula(paste("BinDx ~", paste(variables_seleccionadas, collapse = " + ")))
  
  # Definir matrices
    n <- nrow(datos)
  
    logit_coef_matrix <- matrix(0, nrow = n, ncol = length(variables_seleccionadas) + 1)
    colnames(logit_coef_matrix) <- c("(Intercept)", variables_seleccionadas)
  
    logit_zvalues_matrix <- logit_coef_matrix
    colnames(logit_zvalues_matrix) <- colnames(logit_coef_matrix)
    logit_pvalues_matrix <- logit_coef_matrix
    colnames(logit_pvalues_matrix) <- colnames(logit_coef_matrix)
    
    logit_pred <- numeric(n)
  
  # LOOCV
    for (ss in 1:n) {
      
      # Mostrar progreso
        cat(paste0("..", ss))
      
      # Train n test
        datostrain <- datos[-ss, ]
        datostest <- datos[ss, ]
      
      # Diseño del modelo
        fit <- glm(cv_formula, data = datostrain, family = binomial(link = "logit"))
      
      # Guardar coeficientes
        coefs <- as.numeric(coef(fit))
        logit_coef_matrix[ss,] <- coefs
      
      # Guardar zvalores y pvalores
        sumfit <- summary(fit)$coefficients
        zvalues <- sumfit[, "z value"]
        pvalues <- sumfit[, "Pr(>|z|)"]
      
        logit_zvalues_matrix[ss, names(zvalues)] <- zvalues
        logit_pvalues_matrix[ss, names(pvalues)] <- pvalues
      
      # Vector de predicciones
        logit_pred[ss] <- predict(fit, newdata = datostest, type = "response")
  }
  
  # Curva ROC y AUC
    rocobj <- roc(datos$BinDx, logit_pred)
    print(ggroc(rocobj) + 
            coord_fixed() + 
            labs(title = "Curva ROC", x = "Especificidad", y = "Sensibilidad"))
    cat("AUC:", rocobj$auc, "\n\n")
  
  # Coeficiente medio
    coefs_mean <- colMeans(logit_coef_matrix)
    cat("Coeficientes:\n")
    print(coefs_mean[coefs_mean != 0])
  
  # Medias de valores z
    z_col_means <- colMeans(logit_zvalues_matrix)
    cat("\nMedias de z:\n")
    print(z_col_means[z_col_means != 0])
  
  # Medias de valores p
    p_col_means <- colMeans(logit_pvalues_matrix)
    cat("\nMedias de p:\n")
    print(p_col_means[p_col_means != 0])
  
  # Invisible return
    invisible(list(
      coef_matrix = logit_coef_matrix,
      zvalues_matrix = logit_zvalues_matrix,
      pvalues_matrix = logit_pvalues_matrix,
      preds = logit_pred,
      auc = rocobj$auc,
      coefs_mean = coefs_mean,
      z_col_means = z_col_means,
      p_col_means = p_col_means
    ))
}

#-------Análisis de los modelos logísticos con las variables seleccionadas

# Modelo completo
  analisis_logit(datos0v3_av, cv_vars_av)
  
# Modelo transversal
  analisis_logit(datos0v3_tv, cv_vars_tv)
  
# Modelo longitudinal
  analisis_logit(datos0v3_lv, cv_vars_lv)
