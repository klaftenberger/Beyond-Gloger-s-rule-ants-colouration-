#### Title: Script Figure 2 for the study: "Beyond Gloger’s rule: multiple biogeographic drivers of red and dark colouration in ants"
### Authors: XX
## Analyses/script: XX
# Last update: 17/03/2025

library(ggplot2)   #‘3.5.2’
library(biscale)   #v‘1.1.0’
library(patchwork) #v1.3.0’
library(sf)        #v‘1.0.16’


eco_bi <- bi_class(
      Ecoreg_4_map,
      x = NMDS1,
      y = NMDS2,
      style = "fisher",
      dim = 3
)

p_map <- ggplot() +
      geom_sf(data = ocean, fill = "lightblue1", color = NA, alpha = 0.2) +
      geom_sf(data = eco_bi, aes(fill = bi_class), colour = NA) +
      bi_scale_fill(pal = "DkViolet2", dim = 3,    
                    flip_axes = TRUE,
                    #rotate_pal = TRUE
                    ,) +
      theme_minimal(base_family = "Arial") +
      theme(
            text = element_text(family = "Arial"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none"
      )

p_legend <- bi_legend(
      pal = "DkViolet2",
      dim = 3,
      xlab = "NMDS1",
      ylab = "NMDS2",
      flip_axes = TRUE
      
      #rotate_pal = TRUE,
)

p_map + p_legend +
      plot_layout(widths = c(10, 1))




