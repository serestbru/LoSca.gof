#############################
#### EXAMPLE - POWER ANALYSIS
#############################

Tabla <- mapply(
  function(x1, x2) {
    potencia.cvm(
      n = 30,
      str.distH0   = dist.map[[x1]],
      str.dist.used = dist.map[[x2]],
      params.H0    = parametros[[x1]],
      params.used  = parametros[[x2]],
      method = "MLE"
    )
  },
  x1 = rep(d1, each = length(d2)),
  x2 = rep(d2, times = length(d1))
)

dim(Tabla) <- c(length(d2), length(d1))
rownames(Tabla) <- d2
colnames(Tabla) <- d1
print(Tabla)
print(xtable(Tabla, digits = c(0, rep(4, ncol(Tabla)))))
heatmap.LoScF(Tabla, "cvm30mle.png")