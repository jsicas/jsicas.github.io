#' @title Função De Teste SpaHet
#' 
#' @param n quantidade de pontos


spahet <- function(n=1024) {
  x <- seq(1,n)/n
  spahet <- sqrt(x * (1 - x)) * sin(2 * pi * (1 + 2^(-0.6))/(x + 2^(-0.6)))
  spahet <- spahet/sd(spahet) * 7
  return(spahet)
}
