# simulação entre proposta bayesiana e splines.

# packages
library(doFuture)
library(wavethresh)
library(splines)

# carregando funções de simulação
source('R/spahet.R')
source('R/epanec_shrink.R')
source('R/simu_wave.R')
source('R/simu_spline.R')

# configurações
set.seed(282829)
W_32 <- GenW(32, filter.number=5, family='DaubExPhase')
W_128 <- GenW(128, filter.number=5, family='DaubExPhase')
s_3 <- 7/3
s_7 <- 7/7
l3 <- 1/s_3^2 + exp(-s_3/2)/2
l7 <- 1/s_7^2 + exp(-s_7/2)/2
plan(multisession, workers=6)


# wavelets ----------------------------------------------------------------
# esturatura do nome do objeto: wave_curva_n_SNR_tipo, com tipo abreviado, onde
# 's' representa sure e 'b' representa o bayesiano.
# SNR = 3 | bumps ===================
wave_bumps_32_3_s <- simu_wave(DJ.EX(n=32)$bumps, SNR=3, rep=10000, policy='sure')
wave_bumps_128_3_s <- simu_wave(DJ.EX(n=128)$bumps, SNR=3, rep=10000, policy='sure')
wave_bumps_32_3_b <- simu_wave(DJ.EX(n=32)$bumps, SNR=3, rep=10000, policy='bayes',
                               a=0.8, l=l3)
wave_bumps_128_3_b <- simu_wave(DJ.EX(n=128)$bumps, SNR=3, rep=10000, policy='bayes',
                                a=0.8, l=l3)

# SNR = 7 | bumps ===================
wave_bumps_32_7_s <- simu_wave(DJ.EX(n=32)$bumps, SNR=7, rep=10000, policy='sure')
wave_bumps_128_7_s <- simu_wave(DJ.EX(n=128)$bumps, SNR=7, rep=10000, policy='sure')
wave_bumps_32_7_b <- simu_wave(DJ.EX(n=32)$bumps, SNR=7, rep=10000, policy='bayes',
                               a=0.8, l=l7)
wave_bumps_128_7_b <- simu_wave(DJ.EX(n=128)$bumps, SNR=7, rep=10000, policy='bayes',
                                a=0.8, l=l7)

# SNR = 3 | spahet ==================
wave_spahet_32_3_s <- simu_wave(spahet(n=32)  , SNR=3, rep=10000, policy='sure')
wave_spahet_128_3_s <- simu_wave(spahet(n=128), SNR=3, rep=10000, policy='sure')
wave_spahet_32_3_b <- simu_wave(spahet(n=32)  , SNR=3, rep=10000, policy='bayes',
                                a=0.8, l=l3)
wave_spahet_128_3_b <- simu_wave(spahet(n=128), SNR=3, rep=10000, policy='bayes',
                                 a=0.8, l=l3)

# SNR = 7 | spahet ==================
wave_spahet_32_7_s <- simu_wave(spahet(n=32)  , SNR=7, rep=10000, policy='sure')
wave_spahet_128_7_s <- simu_wave(spahet(n=128), SNR=7, rep=10000, policy='sure')
wave_spahet_32_7_b <- simu_wave(spahet(n=32)  , SNR=7, rep=10000, policy='bayes',
                                a=0.8, l=l7)
wave_spahet_128_7_b <- simu_wave(spahet(n=128), SNR=7, rep=10000, policy='bayes',
                                 a=0.8, l=l7)


# splines -----------------------------------------------------------------
# esturatura do nome do objeto: spline_curva_n_SNR.
# SNR = 3 | bumps ===================
spline_bumps_32_3 <- simu_spline(DJ.EX(n=32)$bumps  , SNR=3, rep=10000)
spline_bumps_128_3 <- simu_spline(DJ.EX(n=128)$bumps, SNR=3, rep=10000)

# SNR = 7 | bumps ===================
spline_bumps_32_7 <- simu_spline(DJ.EX(n=32)$bumps  , SNR=7, rep=10000)
spline_bumps_128_7 <- simu_spline(DJ.EX(n=128)$bumps, SNR=7, rep=10000)

# SNR = 3 | bumps ===================
spline_spahet_32_3 <- simu_spline(spahet(n=32)  , SNR=3, rep=10000)
spline_spahet_128_3 <- simu_spline(spahet(n=128), SNR=3, rep=10000)

# SNR = 7 | bumps ===================
spline_spahet_32_7 <- simu_spline(spahet(n=32)  , SNR=7, rep=10000)
spline_spahet_128_7 <- simu_spline(spahet(n=128), SNR=7, rep=10000)



# Resultados --------------------------------------------------------------
stats <- \(x) round(c(MSE=mean(x$MSE), sd=sd(x$MSE)), 2)

# bumps n=32
list(wave_bumps_32_3_b, wave_bumps_32_7_b,
     wave_bumps_32_3_s, wave_bumps_32_7_s,
     spline_bumps_32_3, spline_bumps_32_7) |>
  sapply(stats) |> 
  t() |> 
  apply(MARGIN=1, \(x) cat(x[1], ' (', x[2], ')', '\n', sep=''))

# bumps n=128
list(wave_bumps_128_3_b, wave_bumps_128_7_b,
     wave_bumps_128_3_s, wave_bumps_128_7_s,
     spline_bumps_128_3, spline_bumps_128_7) |>
  sapply(stats) |> 
  t() |> 
  apply(MARGIN=1, \(x) cat(x[1], ' (', x[2], ')', '\n', sep=''))


# spahet n=32
list(wave_spahet_32_3_b, wave_spahet_32_7_b,
     wave_spahet_32_3_s, wave_spahet_32_7_s,
     spline_spahet_32_3, spline_spahet_32_7) |>
  sapply(stats) |> 
  t() |> 
  apply(MARGIN=1, \(x) cat(x[1], ' (', x[2], ')', '\n', sep=''))

# spahet n=128
list(wave_spahet_128_3_b, wave_spahet_128_7_b,
     wave_spahet_128_3_s, wave_spahet_128_7_s,
     spline_spahet_128_3, spline_spahet_128_7) |>
  sapply(stats) |> 
  t() |> 
  apply(MARGIN=1, \(x) cat(x[1], ' (', x[2], ')', '\n', sep=''))








list(wave_bumps_32_3_b, wave_bumps_32_3_s, spline_bumps_32_3,
     wave_bumps_32_7_s, wave_bumps_32_7_s, spline_bumps_32_7) |>
  sapply(\(x) x$MSE) |> 
  boxplot(main=expression('Bumps com' ~ n == 32), names=rep(c('Bayes', 'SURE', 'Splines'), 2),
          main='', col=rep(c('green', 'blue'), each=3))
legend('right', legend=c(expression(SNR == 3), expression(SNR == 7)),
       col=c('green', 'blue'), bty='n', pch=16, pt.cex=1.2)


list(wave_bumps_128_3_b, wave_bumps_128_3_s, spline_bumps_128_3,
     wave_bumps_128_7_s, wave_bumps_128_7_s, spline_bumps_128_7) |>
  sapply(\(x) x$MSE) |> 
  boxplot(main=expression('Bumps com' ~ n == 128), names=rep(c('Bayes', 'SURE', 'Splines'), 2),
          main='', col=rep(c('green', 'blue'), each=3))
legend('right', legend=c(expression(SNR == 3), expression(SNR == 7)),
       col=c('green', 'blue'), bty='n', pch=16, pt.cex=1.2)

list(wave_spahet_32_3_b, wave_spahet_32_3_s, spline_spahet_32_3,
     wave_spahet_32_7_s, wave_spahet_32_7_s, spline_spahet_32_7) |>
  sapply(\(x) x$MSE) |> 
  boxplot(main=expression('SpaHet com' ~ n == 32), names=rep(c('Bayes', 'SURE', 'Splines'), 2),
          main='', col=rep(c('green', 'blue'), each=3))
legend('right', legend=c(expression(SNR == 3), expression(SNR == 7)),
       col=c('green', 'blue'), bty='n', pch=16, pt.cex=1.2)

list(wave_spahet_128_3_b, wave_spahet_128_3_s, spline_spahet_128_3,
     wave_spahet_128_7_s, wave_spahet_128_7_s, spline_spahet_128_7) |>
  sapply(\(x) x$MSE) |> 
  boxplot(main=expression('SpaHet com' ~ n == 128), names=rep(c('Bayes', 'SURE', 'Splines'), 2),
          main='', col=rep(c('green', 'blue'), each=3))
legend('right', legend=c(expression(SNR == 3), expression(SNR == 7)),
       col=c('green', 'blue'), bty='n', pch=16, pt.cex=1.2)

