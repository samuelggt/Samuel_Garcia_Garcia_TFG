# Clear workspace
rm(list = ls())
  
# Save original parameters
op <- par()

# Instalar paquetes necesarios
  
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("pROC")) install.packages("pROC"); library(pROC)

#-------Cargar Datos

datos <- read.csv("Data/STC_data.csv")

#-------Crear variable binaria

datos$BinDx <- NA
datos$BinDx[datos$Dx<3] <- 0
datos$BinDx[datos$Dx==3] <- 1

#-------Dataframe Variables Electrodiagnóstico

datos0v3 <- datos[, c(names(datos)[14:18], "BinDx")]
eldx_vars <- colnames(datos0v3)[1:5]
eldx_formula <- as.formula(paste("BinDx ~", paste(eldx_vars, collapse=" + ")))

#-------Eliminar filas con NA

datos0v3 <- na.omit(datos0v3[which(datos$BinDx %in% "0" | datos$BinDx %in% "1"),])

#-------Definición de matrices

# Matriz de coeficientes
  eldx_coef_matrix <- matrix(0, nrow = nrow(datos0v3), ncol = length(eldx_vars) + 1)
  colnames(eldx_coef_matrix) <- c("(Intercept)", eldx_vars)

# Matriz de valores z
  eldx_zvalues_matrix <- eldx_coef_matrix
  colnames(eldx_zvalues_matrix) <- colnames(eldx_coef_matrix)

# Matriz de valores p
  eldx_pvalues_matrix <- eldx_coef_matrix
  colnames(eldx_pvalues_matrix) <- colnames(eldx_coef_matrix)

# Vector predicciones
  eldx_pred <- vector("numeric",length = nrow(datos0v3))

#-------Validación cruzada leave-one-out

for (ss in 1:nrow(datos0v3)){
  
  # Mostrar progreso
    cat(paste0("..",ss))
  
  # Train n test
    datostrain <- datos0v3[-ss,]
    datostest <- datos0v3[ss,]
  
  # Diseño del modelo
    fit <- glm(eldx_formula, data = datostrain, family = binomial(link = "logit"))
  
  # Guardar coeficientes
    coefs <- as.numeric(coef(fit))
    eldx_coef_matrix[ss, ] <- coefs
  
  # Guardar zvalores y pvalores
    sumfit <- summary(fit)$coefficients
    zvalues <- sumfit[, "z value"]
    pvalues <- sumfit[, "Pr(>|z|)"]
  
    eldx_zvalues_matrix[ss, names(zvalues)] <- zvalues
    eldx_pvalues_matrix[ss, names(pvalues)] <- pvalues
  
  # Vector de predicciones
    eldx_pred[ss] <- predict(fit, newdata = datostest, type = "response")
  
}

#-------Análisis
  
# Curva ROC y AUC
  rocobj <- roc(datos0v3$BinDx, eldx_pred)
  ggroc(rocobj) + coord_fixed() + labs(title = "Curva ROC", x = "Especificidad", y = "Sensibilidad")
  rocobj$auc

# Coeficiente medio
  colMeans(eldx_coef_matrix)

# Media de los valores z para cada variable
  colMeans(eldx_zvalues_matrix)

# Media de los valores p para cada variable
  colMeans(eldx_pvalues_matrix)
