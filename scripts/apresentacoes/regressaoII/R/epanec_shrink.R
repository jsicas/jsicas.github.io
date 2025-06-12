# -------- Shrinkage da Epanechnikov Generalizada -------- #
# d - coeficientes empiricos.
# a - parâmetro alpha da priori spike and slab.
# b - parâmetro da Epanechnikov generalizada.
# l - lambda (hiperparámetro de sigma^2).

# library(nimble)  # para usar ddexp
epanec_shrink <- function(d, a, b, l){
  num <- (1-a) * (3*sqrt(2*l)/(8*b^3)) * (
    (2*l*b^2 + 3*sqrt(2*l)*b + 3)/(2*l^2) * (exp(-sqrt(2*l)*(b - d)) -
                                               exp(-sqrt(2*l)*(b + d))) +
      ((l*b^2 - 3)*sqrt(2*l)*d - l*sqrt(2*l)*d^3)/(l^2)
  )
  den <- a * sqrt(2*l)/2 * exp(-sqrt(2*l)*abs(d)) + (1-a) * (3*sqrt(2*l)/(8*b^3)) *
    (
      (b/l) * (exp(-sqrt(2*l)*(b + d)) + exp(-sqrt(2*l)*(b - d))) + 2/sqrt(2*l) *
        (b^2 - d^2 - 1/l)
    )
  delta <- num/den
  delta <- ifelse(abs(delta) >= abs(d), d, delta)
  return(delta)
}
