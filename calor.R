#######################################################################
##### heatmap.LoScF
#######################################################################
# Personalised heatmap function to enhance the visualization of the power analysis

heatmap.LoScF <- function(datos, nombre){
  datos <- as.data.frame(datos)
  
  # INIT -------------------------------------------------------------
  datos$Hipotesis <- rownames(datos)
  datos_melt <- melt(datos, id.vars = "Hipotesis")
  colnames(datos_melt) <- c("H0", "H1", "Valor")
  orden <- rownames(datos)
  
  datos_melt$H0 <- factor(datos_melt$H0, levels = rev(orden))  # Y AXIS
  datos_melt$H1 <- factor(datos_melt$H1, levels = orden)       # X AXIS
  # COLOURS PALETTE -------------------------------------------
  escala_colores <- c(
    "#1a9850",  
    "#fee08b",  
    "#d7191c"    
  )
  
  puntos_corte <- c(0, 0.05, 0.3, 1)  # CRITICAL VALUES 
  
  # HEATMAP --------------------------------------------------------
  heatmap_premium <- ggplot(datos_melt, aes(x = H1, y = H0, fill = Valor)) +
    
    geom_tile(color = "white", size = 1.2, width = 0.95, height = 0.95) +
    
    scale_fill_gradientn(
      colours = escala_colores,
      values = scales::rescale(puntos_corte),
      limits = c(0, 1),
      na.value = "gray90",
      breaks = seq(0, 1, 0.2),
      guide = guide_colorbar(
        title = "Potencia",
        title.position = "top",
        title.theme = element_text(size = 12, face = "bold"),
        barwidth = unit(8, "cm"),
        barheight = unit(0.4, "cm"),
        frame.colour = "black",
        ticks = TRUE,
        direction = "horizontal"
      )
    ) +
    
    geom_text(aes(label = ifelse(Valor < 0.001, "<0.001",
                                 ifelse(Valor < 0.01, sprintf("%.3f", Valor),
                                        ifelse(Valor < 0.05, sprintf("%.3f", Valor),
                                               sprintf("%.3f", Valor))))),
              color = "black", 
              size = 5.5,
              fontface = "bold") +
    
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
      axis.title = element_blank(),
      axis.text = element_text(face = "bold", color = "black", size = 13),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      legend.position = "bottom",
      legend.title.align = 0.5,
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      aspect.ratio = 1
    ) +
    
    labs(
      title = "ANALISIS DE POTENCIAS",
      subtitle = "Mapa de calor: potencias obtenidas en las comparaciones por pares",
      caption = ""
    ) +
    
    coord_fixed()
  
  # VISUALIZATION -----------------------------------------------------------
  print(heatmap_premium)
  ggsave(nombre, plot = heatmap_premium,
         width = 9, height = 7, dpi = 600, bg = "white")
}
