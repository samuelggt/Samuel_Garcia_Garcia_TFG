# Samuel_Garcia_Garcia_TFG
Repositorio que contiene los códigos que se han diseñado para el trabajo de fin de grado de matemáticas de la Universidad de Zaragoza: Selección de características y modelos predictivos para el diagnóstico clínico mediante ultrasonido, Samuel García García, Junio 2025.

El archivo "Modelo Variables Electrodiagnósticas.R" contiene el código para el diseño y análisis de un modelo de regresión logística basado en las cinco variables electrodiagnósticas, utilizando validación cruzada leave-one-out (LOOCV).

El archivo "Modelos Lasso Variables US.R" incluye el código para la construcción de tres modelos Lasso con variables de ultrasonido cuantitativo (transversales, longitudinales y la combinación de ambas), aplicando LOOCV. También se desarrollan modelos de regresión logística clásicos utilizando como predictores las variables seleccionadas por Lasso.

El archivo "Modelos Red Elástica Variables US.R" contiene el código destinado al cálculo del AUC de los modelos ecográficos empleando red elástica con diferentes técnicas de muestreo
