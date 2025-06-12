#' @title Função Para Estudo De Simulação Com Splines
#' 
#' @param fun_real função verdadeira.
#' @param SNR razão sinal-ruído desejada para simulação.
#' @param rep número de replicações.

# library(doFuture)
# library(splines)

simu_spline <- function(fun_real, SNR, rep=500) {
  # criando objetos
  MSE <- numeric(rep)
  n <- length(fun_real)
  x <- seq(1, n)/n
  
  result <- foreach(i = 1:rep, .combine='rbind',
                    .options.future=list(seed=TRUE)) %dofuture% {
    # gerando amostra
    sd_ruido <- sd(fun_real)/SNR              # calculando sd(ruido)
    ruido <- rnorm(n, mean=0, sd=sd_ruido)    # erro
    fun_ruido <- fun_real + ruido             # adicionando ruido
    
    # splines
    fit <- smooth.spline(x, fun_ruido, cv=F)
    
    # MSE e percentual de coeficientes nulos
    residuo <- fun_real - predict(fit)$y
    MSE <- mean(residuo^2)  # sum(\hat{g} - g)/n
    c('MSE' = MSE)
  }
  
  return(data.frame(result))
}
