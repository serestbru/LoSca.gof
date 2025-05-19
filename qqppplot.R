# Librerías 
library(ggplot2)
library(gridExtra)

# 1. Generar datos de ejemplo, normLS
set.seed(123)  # Para reproducibilidad
datos <- rnormLS(n = 100, mu = 5, sigma = 2)

# Parámetrosmáxima verosimilitud
mu_est <- mean(datos)
sigma_est <- sd(datos)

# PP-plot (Probabilidad vs Probabilidad teórica)
pp_data <- data.frame(
  teorica = ppoints(length(datos)),  # Probabilidades teóricas
  empirica = sort(pnorm(datos, mu_est, sigma_est))  # Probabilidades empíricas
)

pp_plot <- ggplot(pp_data, aes(x = teorica, y = empirica)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "PP-plot: Normalidad",
       x = "Probabilidad teórica",
       y = "Probabilidad empírica") +
  theme_minimal()

# QQ-plot (Cuantiles empíricos vs teóricos)
qq_plot <- ggplot(data.frame(datos), aes(sample = datos)) +
  geom_qq(distribution = qnorm, dparams = list(mean = mu_est, sd = sigma_est)) +
  geom_qq_line(distribution = qnorm, dparams = list(mean = mu_est, sd = sigma_est), 
               color = "red", linetype = "dashed") +
  labs(title = "QQ-plot: Normalidad",
       x = "Cuantiles teóricos",
       y = "Cuantiles empíricos") +
  theme_minimal()

grid.arrange(pp_plot, qq_plot, ncol = 2)
