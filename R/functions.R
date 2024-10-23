#' Cálculo de Tamaño Muestral y Potencia para Análisis de Supervivencia (Modelo de Cox)
#' 
#' @param alfa Nivel de significación
#' @param potencia Potencia deseada (opcional si se proporciona n)
#' @param n Tamaño muestral (opcional si se proporciona potencia)
#' @param HR Hazard Ratio a detectar
#' @param IC_superior Límite superior del intervalo de confianza del HR
#' @param IC_inferior Límite inferior del intervalo de confianza del HR (opcional si se proporciona EE)
#' @param EE Error estándar del coeficiente (opcional si se proporcionan los IC)
#' @param n_previo Tamaño de muestra del estudio previo (requerido para método "ee")
#' @param sigma_x Desviación estándar del predictor X (requerido para método "rho")
#' @param psi Probabilidad de que una observación no esté censurada (requerido para método "rho")
#' @param metodo Método de cálculo: "rho" o "ee" (error estándar)
#' @param valores_rho Vector de valores rho para calcular tamaños muestrales (si metodo = "rho")
#' @return Un data frame con resultados según el método elegido
#' @export
SampleHazardRatioMulti <- function(alfa, potencia = NULL, n = NULL, 
                                   HR = NULL,
                                   IC_superior = NULL, IC_inferior = NULL,
                                   EE = NULL, n_previo = NULL,
                                   sigma_x = NULL, psi = NULL,
                                   metodo = "rho",
                                   valores_rho = seq(0, 0.9, by = 0.1)) {
  
  # Verificaciones según el método
  if (metodo == "rho") {
    if (is.null(sigma_x) || is.null(psi)) {
      stop("Para método 'rho' se requiere sigma_x y psi")
    }
    if (is.null(HR)) {
      stop("Se requiere especificar HR")
    }
  } else if (metodo == "ee") {
    if (is.null(HR)) {
      stop("Para método 'ee' se requiere HR")
    }
    if (is.null(EE) && (is.null(IC_superior) || is.null(IC_inferior))) {
      stop("Para método 'ee' se requiere EE o ambos límites del IC")
    }
    if (is.null(n_previo)) {
      stop("Para método 'ee' se requiere n_previo")
    }
  } else {
    stop("El método debe ser 'rho' o 'ee'")
  }
  
  # Si se proporcionaron los IC, calcular el EE
  if (!is.null(IC_superior) && !is.null(IC_inferior)) {
    EE <- (log(IC_superior) - log(HR)) / qnorm(0.975)
    cat("Error Estándar calculado:", EE, "\n")
  }
  
  if (metodo == "rho") {
    beta_a <- log(HR)
    z_alfa <- qnorm(1 - alfa/2)
    
    if (!is.null(potencia)) {
      z_gamma <- qnorm(potencia)
      
      calcular_n <- function(rho) {
        n <- ((z_alfa + z_gamma)^2) / ((beta_a * sigma_x)^2 * psi * (1 - rho^2))
        return(ceiling(n))
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        tamaño_muestral = sapply(valores_rho, calcular_n)
      )
      
      resultados$eventos_esperados <- ceiling(resultados$tamaño_muestral * psi)
    } else {
      calcular_potencia <- function(rho) {
        z_gamma <- sqrt(n * (beta_a * sigma_x)^2 * psi * (1 - rho^2)) - z_alfa
        potencia <- pnorm(z_gamma)
        return(potencia)
      }
      
      resultados <- data.frame(
        rho = valores_rho,
        potencia = sapply(valores_rho, calcular_potencia)
      )
      
      resultados$tamaño_muestral <- n
      resultados$eventos_esperados <- ceiling(n * psi)
    }
  } else {
    # Método basado en error estándar
    z_alfa <- qnorm(1 - alfa/2)
    if (!is.null(potencia)) {
      z_gamma <- qnorm(potencia)
      n <- ceiling((z_alfa + z_gamma)^2 * n_previo * EE^2 / (log(HR))^2)
      resultados <- data.frame(
        tamaño_muestral = n
      )
    } else {
      z_gamma <- sqrt(n * (log(HR))^2 / (n_previo * EE^2)) - z_alfa
      potencia <- pnorm(z_gamma)
      resultados <- data.frame(
        potencia = potencia,
        tamaño_muestral = n
      )
    }
  }
  
  return(resultados)
}

#' Cálculo Logístico para Estudios de Supervivencia
#'
#' @param n_total Tamaño muestral total requerido
#' @param eventos_esperados Número de eventos esperados
#' @param tasa_reclutamiento Tasa de reclutamiento por mes
#' @param tiempo_seguimiento Tiempo de seguimiento en meses
#' @param tasa_perdida_seguimiento Tasa esperada de pérdida en el seguimiento
#' @param tasa_rechazo Tasa de rechazo esperada
#' @return Una lista con los cálculos logísticos del estudio
#' @export
logistica_estudio_supervivencia <- function(n_total, 
                                           eventos_esperados,
                                           tasa_reclutamiento,
                                           tiempo_seguimiento,
                                           tasa_perdida_seguimiento,
                                           tasa_rechazo) {
  
  # Ajuste por pérdidas en el seguimiento
  muestra_con_perdidas <- n_total / (1 - tasa_perdida_seguimiento)
  
  # Ajuste por rechazos
  muestra_para_contactar <- muestra_con_perdidas / (1 - tasa_rechazo)
  
  # Cálculo de tiempos
  tiempo_reclutamiento <- ceiling(muestra_para_contactar / tasa_reclutamiento)
  tiempo_total_estudio <- tiempo_reclutamiento + tiempo_seguimiento
  
  # Resultados
  resultados <- list(
    tamaño_muestral = n_total,
    eventos_esperados = eventos_esperados,
    muestra_con_perdidas = ceiling(muestra_con_perdidas),
    muestra_para_contactar = ceiling(muestra_para_contactar),
    meses_reclutamiento = tiempo_reclutamiento,
    meses_total_estudio = tiempo_total_estudio
  )
  
  # Imprimir resumen
  cat("
Resumen de logística del estudio:
")
  cat("----------------------------------------
")
  cat("Tamaño muestral requerido:", n_total, "
")
  cat("Eventos esperados:", eventos_esperados, "
")
  cat("Muestra considerando pérdidas:", ceiling(muestra_con_perdidas),
      "(", tasa_perdida_seguimiento*100, "% de pérdidas)
")
  cat("Muestra a contactar:", ceiling(muestra_para_contactar),
      "(", tasa_rechazo*100, "% de rechazo)
")
  cat("Meses de reclutamiento:", tiempo_reclutamiento,
      "(", tasa_reclutamiento, "personas por mes)
")
  cat("Meses totales de estudio:", tiempo_total_estudio,
      "(incluyendo", tiempo_seguimiento, "meses de seguimiento)
")
  cat("----------------------------------------
")
  
  return(invisible(resultados))
}
