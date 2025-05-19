library(ggplot2)
library(gridExtra)

# Valores para el eje x
x <- seq(-5, 5, length.out = 1000)

plot_distribution <- function(title, func, params, x_range = c(-5, 5)) {
  df <- data.frame()
  for (i in seq_along(params$mu)) {
    y <- do.call(func, list(x, params$mu[i], params$sigma[i]))
    df <- rbind(df, data.frame(x = x, y = y, 
                               group = paste0("μ=", params$mu[i], ", σ=", params$sigma[i])))
  }
  
  ggplot(df, aes(x = x, y = y, color = group)) +
    geom_line(linewidth = 1) +
    labs(title = title,
         x = "x",
         y = "F(x)",
         color = "Parámetros") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    coord_cartesian(xlim = x_range) +
    theme(legend.position = "bottom")
}

## Parámetros
params_norm <- list(mu = c(0, 0, 0, 1, -1), 
                    sigma = c(2, 1, 0.5, 1, 0.5))
params_logis <- list(mu = c(0, 0, 0, 1, -1), 
                     sigma = c(2, 1, 0.5, 1, 0.5))
params_exp <- list(mu = c(0, 0, 0, 1, -1), 
                   sigma = c(2, 1, 0.5, 1, 0.5))

params_laplace <- list(mu = c(0, 0, 0, 1, -1), 
                       sigma = c(2, 1, 0.5, 1, 0.5))

## Gráficos
p_norm <- plot_distribution("LoScF Normal", pnormLS, params_norm)
p_logis <- plot_distribution("LoScF Logística", plogisLS, params_logis)
p_exp <- plot_distribution("LoScF Exponencial", pexpLS, params_exp)
p_laplace<- plot_distribution("LoScF de Laplace", plaplaceLS, params_laplace)

grid.arrange(p_norm, p_logis, p_laplace, p_exp, ncol = 1)

