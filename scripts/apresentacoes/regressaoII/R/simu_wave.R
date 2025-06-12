#' @title Função Para Estudo De Simulação Com Ondaletas
#' 
#' @param fun_real função verdadeira.
#' @param SNR razão sinal-ruído desejada para simulação.
#' @param rep número de replicações.
#' @param policy política de escolha de limiar. Possíveis são 'sure' e 'bayes'
#'   que executará a regra bayesiana.
#' @param a,l alpha da priori e l de hiperparámetro da
#'   logística.
#' @param type tipo de limiar.
#' @param filter.number parâmetro da função wd.
#' @param family parâmetro da função wd.

# library(doFuture)

simu_wave <- function(fun_real, SNR, rep=500, policy=c('sure', 'bayes'), a, l,
                      filter.number=5, family='DaubExPhase') {
  policy <- match.arg(policy)
  if ((missing(a) | missing(l)) & policy == 'bayes') {
    stop('Especifique os parâmetros da priori.')
  }
  if (policy == 'sure') {
    a <- l <- NULL
  }
  
  # criando objetos
  MSE <- numeric(rep)
  PCN <- numeric(rep)
  n <- length(fun_real)
  W <- GenW(n, filter.number=5, family='DaubExPhase')
  
  result <- foreach(i = 1:rep, .combine = 'rbind',
                    .options.future=list(seed=TRUE, globals=TRUE)) %dofuture% {
    # gerando amostra
    sd_ruido <- sd(fun_real)/SNR              # calculando sd(ruido)
    ruido <- rnorm(n, mean=0, sd=sd_ruido)    # erro
    fun_ruido <- fun_real + ruido             # Adicionando ruido
    
    # wavelet
    ywt <- wd(fun_ruido, filter.number=filter.number, family=family)  # DWT
    if (policy == 'sure') {
      ywt_T <- threshold(ywt, policy=policy, value=lambda)
      fun_estimada <- wr(ywt_T)  # IDWT
    } else {
      si <- mad(accessD(ywt, lev=nlevelsWT(ywt)-1))
      l <- 1/si^2 + exp(-si/2)/2
      beta <- max(abs(ywt$D))
      ywt_T <- c(accessC(ywt, lev=0), epanec_shrink(ywt$D, a=a, b=beta, l=l))
      fun_estimada <- W %*% ywt_T # IDWT
    }
    
    # MSE e percentual de coeficientes nulos
    residuo <- fun_real - fun_estimada
    MSE <- mean(residuo^2)  # sum(\hat{g} - g)/n
    PCN <- ifelse(policy == 'sure', 100*mean(ywt_T$D == 0),
                  100*mean(ywt_T[-1] == 0))
    c('MSE' = MSE, 'PCN' = PCN)
  }
  
  return(data.frame(result))
}
