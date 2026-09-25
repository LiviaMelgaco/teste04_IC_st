####### Simulação de Monte Carlo - GLARMA(1,0) Poisson #######
####### Comparando phi = 0.3 e phi = 0.8                #######


#garantindo reprodutibilidade
RNGkind(kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection")



install.packages("pacman")
library("pacman")
pacman::p_load(glarma, renv)

set.seed(222)


##### PARAMETROS #####

n     <- 200          # tamanho de cada série simulada
beta0 <- 1
beta1 <- 0.5
beta2 <- -0.1
phis  <- c(0.3, 0.8)   # cenários de phi a comparar

# número de réplicas de Monte Carlo, um valor DIFERENTE por cenário de phi
Rs <- setNames(c(25, 50), as.character(phis))
# Rs[["0.3"]] == 1000   |   Rs[["0.8"]] == 5000


##### FUNÇAO GERADORA GLARMA ######

## Função: gera UMA série GLARMA(1,0) e ajusta o modelo,
# devolvendo os coeficientes estimados (ou NA se não convergir)

simula_e_ajusta <- function(n, beta0, beta1, beta2, phi) {
  
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- rnorm(n, mean = 0, sd = 1)
  
  z   <- numeric(n)
  eta <- numeric(n)
  mu  <- numeric(n)
  y   <- numeric(n)
  e   <- numeric(n)
  
  z_ant <- 0
  e_ant <- 0
  
  for (t in 1:n) {
    z[t]   <- phi * (z_ant + e_ant)
    eta[t] <- beta0 + beta1 * x1[t] + beta2 * x2[t] + z[t]
    mu[t]  <- exp(eta[t])
    y[t]   <- rpois(1, lambda = mu[t])
    e[t]   <- (y[t] - mu[t]) / sqrt(mu[t])
    
    z_ant <- z[t]
    e_ant <- e[t]
  }
  
  X_matriz <- cbind(intercept = 1, x1 = x1, x2 = x2)
  
  ajuste <- tryCatch(  #tryCatch se der erro, retorna NULL
    suppressWarnings(    
      glarma(y, X_matriz,
             type      = "Poi",
             method    = "FS",
             residuals = "Pearson",
             phiLags   = 1)
    ),
    error = function(cnd) NULL
  )
  
  # se o glarma não convergir/der erro para essa réplica, devolve NA
  if (is.null(ajuste)) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }
  
  coefs <- unlist(coef(ajuste)) #unlist() achata a lista preservando a ordem.
  
  if (length(coefs) != 4) {  #só garantindo que se der algo errado, vai virar NA tb
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }
  
  coefs <- as.numeric(coefs)
  coefs[c(2, 3, 4, 1)]  # reordena: beta0, beta1, beta2, phi
}


##### MONTE CARLO #####

resultados <- list()

for (phi_val in phis) {
  
  R_val <- Rs[[as.character(phi_val)]]  # número de réplicas DESSE phi
  
  cat("Rodando Monte Carlo para phi =", phi_val, "com R =", R_val, "réplicas\n")
  
  estimativas <- matrix(NA_real_, nrow = R_val, ncol = 4)
  colnames(estimativas) <- c("beta0", "beta1", "beta2", "phi")
  
  for (r in 1:R_val) {
    estimativas[r, ] <- simula_e_ajusta(n, beta0, beta1, beta2, phi_val)
    if (r %% 50 == 0) cat("  réplica", r, "de", R_val, "\n") #tipo uma barra de progressão, a cada 50 réplicas
  }
  
  resultados[[as.character(phi_val)]] <- estimativas
}

##### EQM E MEDIAS #####

tabela_larga <- data.frame(parametro = c("beta0", "beta1", "beta2", "phi"))

for (phi_val in phis) {
  
  R_val <- Rs[[as.character(phi_val)]]
  est   <- resultados[[as.character(phi_val)]]
  vt    <- c(beta0, beta1, beta2, phi_val)  # valor verdadeiro de cada parâmetro nesse cenário
  
  #filtrando as NA's do calculo e da tabela
  est_validas <- est[stats::complete.cases(est), , drop = FALSE] 
  n_validas   <- nrow(est_validas)
  
  medias <- colMeans(est_validas)
  eqm    <- colMeans(sweep(est_validas, 2, vt, FUN = "-")^2) 
  
  tabela_larga[[paste0("media_phi", phi_val)]] <- round(medias, 4)
  tabela_larga[[paste0("eqm_phi",   phi_val)]] <- round(eqm, 4)
  
  cat("phi =", phi_val, "->", n_validas, "de", R_val, "réplicas convergiram\n")
}

##### TABELA E GRAFICO FINAL #####

## tabela
cat("\nValores verdadeiros dos betas: beta0 =", beta0,
    ", beta1 =", beta1, ", beta2 =", beta2, "\n\n")

print(tabela_larga)

##boxplot

par(mfrow = c(1, 3))
for (i in 1:3) {
  boxplot(
    resultados[["0.3"]][, i], resultados[["0.8"]][, i],
    names = c("phi = 0.3", "phi = 0.8"),
    main = colnames(estimativas)[i]
  )
  abline(h = c(beta0, beta1, beta2)[i], col = "red", lty = 2)
}


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")