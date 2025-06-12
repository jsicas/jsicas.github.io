# Testes
# Distribuição Epanechnikov Generalizada
# b - beta

depanec <- function(x, b) {
  ifelse(-b < x & x < b, 3/(4*b^3) *(b^2 - x^2), 0)
}

# -------- Shrinkage da Epanechnikov Generalizada -------- #
# d - coeficientes empiricos
# a - alpha
# b - parâmetro da Epanechnikov generalizada
# l - lambda (hiperparámetro de sigma^2)

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



# debug(epanec_shrink)
epanec_shrink(1, 0.7, 3, 2)
undebug(epanec_shrink)
epanec_shrink(1:3, 0.7, 3, 2)

curve(epanec_shrink(x, a=0.7, b=3, l=5), -5,5)




library(extraDistr)
delta_star <- function(d,beta,lambda,alpha) {
  
  num <- (1 - alpha) * (3 * sqrt(2 * lambda) / (8 * beta^3)) * (
    ((2 * lambda * beta^2 + 3 * sqrt(2 * lambda) * beta + 3) / (2 * lambda^2)) *
      (exp(-sqrt(2 * lambda) * (beta - d)) - exp(-sqrt(2 * lambda) * (beta + d))) +
      ((lambda * beta^2 - 3) * sqrt(2 * lambda) * d - lambda * sqrt(2 * lambda) * d^3) / lambda^2
  )
  
  
  den <- alpha * dlaplace(x = d, mu = 0, sigma =1/(sqrt(2*lambda ))) + (1 - alpha) * (3 * sqrt(2 * lambda) / (8 * beta^3)) * (
    (beta / lambda) * (exp(-sqrt(2 * lambda) * (beta + d)) + exp(-sqrt(2 * lambda) * (beta - d))) +
      (2 / sqrt(2 * lambda)) * (beta^2 - d^2 - 1 / lambda)
  )
  
  delta <- num / den
  if(abs(delta)>=abs(d)){
    return(d)
  }
  
  return(delta)
  
}


y <- numeric(5000)
x <- seq(-6, 6, length=5000)
for (i in 1:5000) {
  y[i] <- delta_star(x[i], beta = 5, lambda = 1/s^2 + 1/2 * exp(-s/2), alpha = 0.9)
}
plot(x, y, type='l')

y == epanec_shrink(x, a=0.8, b=2, l=2)

curve(delta_star(x, beta, lambda, alpha), from=-6, to=6)






logis_shrink <- function(d, a, s, t) {
  u <- rnorm(10000)
  delta <- vector(mode='double', length=length(d))
  
  for (i in 1:length(d)) {
    x <- s*u + d[i]
    int1 <- mean(x * dlogis(x, scale=t))
    int2 <- mean(dlogis(x, scale=t))
    delta[i] <- (1-a) * int1/(a * dnorm(d[i], sd=s)/s + (1-a) * int2)
  }
  
  return(delta)
}

# -------- Shrinkage da Beta Generalizada -------- #
beta_gen <- function(x, a, m) {
  ((m^2 - x^2)^(a-1))/((2*m)^(2*a - 1) * beta(a,a))
}

# d - coeficientes empiricos
# alpha - coeficiente alpha
# s - desvio padrão
# a,m - parâmetros da beta
beta_shrink <- function(d, alpha, s, m, a) {
  u <- rnorm(10000)
  delta <- vector(mode='double', length=length(d))
  
  for (i in 1:length(d)) {
    x <- s*u + d[i]
    int1 <- mean(x * dbe(x, ))
    int2 <- mean(dlogis(x, scale=t))
    delta[i] <- (1-a) * int1/(a * dnorm(d[i], sd=s)/s + (1-a) * int2)
  }
}











