S_0 <- 0                # Startpunkt
b <- log(4)             # oberer Absorptionsrand
a <- log(1e-12)         # unterer Absorptionsrand (nahe ln(0))
n_trials <- 100000      # Anzahl der Simulationen

p_grid <- seq(0.01, 0.99, by = 0.025)      # p für Simulationen ohne 0
p_grid_0 <- seq(0.00, 1.00, by = 0.05)     # p für Simulationen mit 0
q_grid_0 <- seq(0.00, 1.00, by = 0.05)     # q für Simulationen mit 0

## Schritte für den Fall ohne 0 (immer ein auswählen und den Rest auskommentieren, das Gleiche mit wkt/wkt_0)
schritte <- c(log(0.5), log(2))                                            # Beispiel 3.4.1
# schritte <- c(log(1/3),log(0.5), log(2),log(3))                          # Beispiel 3.4.2
# schritte <- c(log(1/5),log(1/4),log(1/3),log(0.5), log(2),log(3))        # Beispiel 3.4.3
# schritte <- c(log(1/5), log(2),log(3), log(4))                           # Beispiel 3.4.4

## Schritte für den Fall mit 0 (immer ein auswählen und den Rest auskommentieren)
schritte_0 <- c(log(0), log(0.5), log(2))                                      # Beispiel 3.4.5
# schritte_0 <- c(log(0),log(1/5),log(1/4),log(1/3),log(0.5), log(2),log(3))   # Beispiel 3.4.6


# Funktion der Simulation
simulation <- function(schritte, wkt, n_trials, S_0, a, b){
  erfolg_count <- 0
  total_steps <- 0

  for(i in seq_len(n_trials)){
    S_n <- S_0
    steps_count <- 0
    while((S_n > a && S_n < b)){
      X_j <- sample(schritte, size = 1, prob = wkt)
      S_n <- S_n + X_j
      steps_count <- steps_count + 1
    }
    if (S_n >= b) erfolg_count <- erfolg_count + 1
    total_steps <- total_steps + steps_count
  }
  list(
    erfolg_wkt = erfolg_count/n_trials,
    erwartete_absorptionszeit = total_steps/n_trials
  )
}


### Fall ohne 0 ###
absorptionswkt <- numeric(length(p_grid))      # Vektor zur Speicherung der Absorptionswahrscheinlichkeiten
absorptionszeit <- numeric(length(p_grid))     # Vektor zur Speicherung der Absorptionszeiten

# Die Simulation
system.time({
  for(i in seq_along(p_grid)){
    p <- p_grid[i]
    wkt <- c(1-p,p)                                                         # Beispiel 3.4.1
    # wkt <- c((1 - p)/2, (1 - p)/2, p/2, p/2)                              # Beispiel 3.4.2
    # wkt <- c((1 - p)/4, (1 - p)/4, (1 - p)/4, (1 - p)/4, p/2, p/2)        # Beispiel 3.4.3
    # wkt <- c(1 - p, p/3, p/3, p/3)                                        # Beispiel 3.4.4
    res <- simulation(schritte, wkt, n_trials, S_0, a, b)
    absorptionswkt[i] <- res$erfolg_wkt
    absorptionszeit[i] <- res$erwartete_absorptionszeit
    print(p)
  }
})

## Plots
par(mar = c(5, 5.5, 2, 2) + 0.1) 

# Absorptionswahrscheinlichkeit
plot(p_grid, absorptionswkt, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(P(S[T] == b)),
     bty = "l",cex.axis = 1.5, cex.lab = 1.4)

# Erwartete Absorptionszeit
plot(p_grid, absorptionszeit, type = "l", lwd = 2, col = "deepskyblue3",
     xlab = expression(p), ylab = expression(E(T)),
     bty = "l",cex.axis = 1.5, cex.lab = 1.4)


### Fall mit 0 ###
absorptionswkt_mit_0  <- matrix(NA_real_, nrow = length(q_grid_0), ncol = length(p_grid_0),
                     dimnames = list(q = q_grid_0, p = p_grid_0))                                  
absorptionszeit_mit_0 <- matrix(NA_real_, nrow = length(q_grid_0), ncol = length(p_grid_0),
                     dimnames = list(q = q_grid_0, p = p_grid_0))

# Die Simulation
system.time({
  for(ip in seq_along(p_grid_0)){
    p <- p_grid_0[ip]
    for(iq in seq_along(q_grid_0)){
      q <- q_grid_0[iq]
      if (p+q > 1){next}
      wkt_0 <- c(1-p-q, q, p)                             # Beispiel 3.4.5
      # wkt_0 <- c(1-p-q, q/4, q/4, q/4, q/4, p/2, p/2)   # Beispiel 3.4.6
      if (any(wkt_0 < -1e-12)) {next}
      wkt_0[wkt_0 < 0] <- 0
      wkt_0 <- wkt_0/sum(wkt_0)
      res <- simulation(schritte_0, wkt_0, n_trials, S_0, a, b)
      absorptionswkt_mit_0[ip,iq] <- res$erfolg_wkt
      absorptionszeit_mit_0[ip,iq] <- res$erwartete_absorptionszeit
    }
  print(p)
  }
})

## 3D-Plots
par(mar = c(3, 3, 2, 1) + 0.1)

# Absorptionswahrscheinlichkeit
persp(x = p_grid_0, y = q_grid_0, z = absorptionswkt_mit_0,
      theta = 20, phi = 25, expand = 0.7,
      xlab = "p", ylab = "q", zlab = "P(S_T=b)",
      ticktype = "detailed", col = "deepskyblue3",
      cex.axis = 1.5, cex.lab = 1.4)

# Erwartete Absorptionszeit
persp(x = p_grid_0, y = q_grid_0, z = absorptionszeit_mit_0,
      theta = 15, phi = 25, expand = 0.7,
      xlab = "p", ylab = "q", zlab = "E(T)",
      ticktype = "detailed", col = "deepskyblue3",cex.axis = 1.5, cex.lab = 1.4)


