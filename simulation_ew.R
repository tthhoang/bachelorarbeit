### Liste von Matrizen erstellen ###
eintraege = c(-4, -3, -2, -1, -0.5, -1/3, -1/4, 1/4, 1/3, 0.5, 1, 2, 3, 4)

matrix_liste <- function(werte = c("kleiner", "groesser"), anzahl, eintraege, det_range = NULL, s_max = NULL, max_tries = 100000) {
  werte <- match.arg(werte)
  liste <- vector("list", anzahl)
  n_found <- 0
  tries <- 0
  
  while (n_found < anzahl && tries < max_tries) {
    tries <- tries + 1
    M <- matrix(sample(eintraege, 4, replace = TRUE), 2, 2)
    ew <- eigen(M, only.values = TRUE)$values
    betrag <- Mod(ew)
    if ((werte == "kleiner" && !all(betrag < 1)) ||
        (werte == "groesser" && !all(betrag > 1))) next
    
    # Determinantenbedingung
    if (!is.null(det_range)) {
      d <- det(M)
      if (d < det_range[1] || d > det_range[2]) next
    }
    
    # Singulärwertbedingung
    if (!is.null(s_max)) {
      s <- svd(M, nu = 0, nv = 0)$d
      if (max(s) > s_max) next
    }
    n_found <- n_found + 1
    liste[[n_found]] <- M
  }
  return(liste)
}

# 300 Matrizen mit |Eigenwerten| < 1 
kleiner_liste <- matrix_liste("kleiner", anzahl = 300, eintraege = eintraege, det_range = c(-2, 2), s_max = 4)
# 300 Matrizen mit |Eigenwerten| > 1 
groesser_liste <- matrix_liste("groesser", anzahl = 300, eintraege = eintraege, det_range = c(-2, 2), s_max = 4)
zero <- list(matrix(c(0,0,0,0), 2, 2)) # Nullmatrix für den Fall mit 0



### Simulation ###
x0 <- c(1, 1)           # Startpunkt
tau <- 4                # äußerer Absorptionsrand
a = 1e-12               # 0
n_trials <- 100000      # Anzahl der Simulationen

p_grid <- seq(0.01, 0.99, by = 0.025)       # p für den Fall ohne 0    
p_grid_0 <- seq(0.00, 1, by = 0.05)         # p für den Fall mit 0
q_grid_0 <- seq(0.00, 1, by = 0.05)         # q für den Fall mit 0


# Fuktion der Simulation für den Fall ohne 0
simulation <- function(p, kleiner_liste, groesser_liste, n_trials, x0, a, tau) {
  erfolg_count <- 0
  total_steps   <- 0
  
  for (i in seq_len(n_trials)) {
    x <- x0
    steps_count <- 0
    while (!(all(abs(x) < a)) && all(abs(x) <= tau)) {
      if (runif(1) < p) {
        M <- groesser_liste[[ sample.int(length(groesser_liste), 1) ]]
      } else {
        M <- kleiner_liste[[ sample.int(length(kleiner_liste), 1) ]]
      }
      x <- as.vector(M %*% x)
      steps_count <- steps_count + 1
    }
    if (any(abs(x) > tau)) {
      erfolg_count <- erfolg_count + 1
    }
    total_steps <- total_steps + steps_count
  }
  list(
    erfolg_wkt = erfolg_count / n_trials,
    erwartete_absorptionszeit = total_steps / n_trials
  )
}


# Fuktion der Simulation für den Fall mit 0
simulation_0 <- function(p, q, kleiner_liste, groesser_liste, n_trials, zero, x0, a, tau) {
  erfolg_count <- 0L
  total_steps  <- 0L
  
  for (i in seq_len(n_trials)) {
    x <- x0
    steps_count <- 0L
    while (!(all(abs(x) < a)) && all(abs(x) <= tau)) {
      u <- runif(1)
      if (u < p) {
        M <- groesser_liste[[ sample.int(length(groesser_liste), 1) ]]
      } else if (u < p + q) {
        M <- kleiner_liste[[ sample.int(length(kleiner_liste), 1) ]]
      } else {
        M <- zero[[1L]]  # Nullmatrix
      }
      x <- as.vector(M %*% x)
      steps_count <- steps_count + 1
    }
    if (any(abs(x) > tau)) {
      erfolg_count <- erfolg_count + 1
    }
    total_steps <- total_steps + steps_count
  }
  list(
    erfolg_wkt = erfolg_count/n_trials,
    erwartete_absorptionszeit = total_steps/n_trials
  )
}


### Fall ohne 0 ###
absorptionswkt <- numeric(length(p_grid))
absorptionszeit <- numeric(length(p_grid))

# Die Simulation
system.time({
  for (i in seq_along(p_grid)) {
    p <- p_grid[i]
    res <- simulation(p, kleiner_liste, groesser_liste, n_trials, x0, a, tau)
    absorptionswkt[i] <- res$erfolg_wkt
    absorptionszeit[i] <- res$erwartete_absorptionszeit
    print(p)
  }
})

## Plots
par(mar = c(5, 5, 4, 2) + 0.1) 

# Absorptionswahrscheinlichkeit
plot(p_grid, absorptionswkt, type = "l", lwd = 2,
     xlab = expression(p), ylab = expression(P(S[T] == b)),
     col = "deepskyblue3",bty = "l", cex.axis = 1.5, cex.lab  = 1.4)

# Erwartete Absorptionszeit
plot(p_grid, absorptionszeit, type = "l", lwd = 2,
     xlab = expression(p), ylab = expression(E(T)),
     col = "deepskyblue3", bty = "l", cex.axis = 1.5, cex.lab  = 1.4)




### Fall mit 0 ###
absorptionswkt_mit_0  <- matrix(NA_real_, nrow = length(q_grid_0), ncol = length(p_grid_0),
                                dimnames = list(q = q_grid_0, p = p_grid_0))                                  
absorptionszeit_mit_0 <- matrix(NA_real_, nrow = length(q_grid_0), ncol = length(p_grid_0),
                                dimnames = list(q = q_grid_0, p = p_grid_0))

# die Simulation
system.time({
  for (ip in seq_along(p_grid_0)) {
    p <- p_grid_0[ip]
    for (iq in seq_along(q_grid_0)) {
      q <- q_grid_0[iq]
      if (p + q <= 1) {
        res_0 <- simulation_0(p, q, kleiner_liste, groesser_liste, n_trials, zero, x0, a, tau)
        absorptionswkt_mit_0[ip, iq]  <- res_0$erfolg_wkt
        absorptionszeit_mit_0[ip, iq] <- res_0$erwartete_absorptionszeit
      }
    }
    print(p)
  }
})

## 3D-Plots
par(mar = c(3, 3, 2, 1) + 0.1)

# Absorptionswahrscheinlichkeit
persp(x = p_grid_0, y = q_grid_0, z = absorptionswkt_mit_0,
      theta = 20, phi = 30, expand = 0.8,
      xlab = "p", ylab = "q", zlab = "P(S_T=b)",
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.3, cex.lab = 1.3)

# Erwartete Absorptionszeit
persp(x = p_grid_0, y = q_grid_0, z = absorptionszeit_mit_0,
      theta = 20, phi = 20, expand = 0.8,
      xlab = "p", ylab = "q", zlab = expression(E(T)),
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.3, cex.lab = 1.3)