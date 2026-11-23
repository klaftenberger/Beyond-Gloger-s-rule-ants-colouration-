#### Title: SCript Figure 2 for the study: "Beyond Gloger’s rule: multiple biogeographic drivers of red and dark colouration in ants"
### Authors: XX
## Analyses/script: XX
# Last update: 17/03/2025

library(dplyr)     #v‘1.1.4’
library(ggplot2)   #v‘3.5.2’
library(gtable)    #v‘0.3.6’
library(grid)      #v‘4.4.1’
library(patchwork) #v‘1.3.0’

# PARAMETERS (same limits as in the PCA from the main script)

xlim_use <- c(-6, 5)
ylim_use <- c(-6, 5)

# Pre-compute density ranges to avoid clipping
#    (we fix axis limits with a small margin)

d1 <- density(ind_df$Dim1, na.rm = TRUE)
d2 <- density(ind_df$Dim2, na.rm = TRUE)
ymax1 <- max(d1$y) * 1.10
ymax2 <- max(d2$y) * 1.10

# TOP SMOOTH CURVE: density of PCA axis 1 (Dim1)

dens_plot_x <- ggplot(ind_df, aes(x = Dim1)) +
      geom_density(linewidth = 0.8, color = "black", fill = NA, adjust = 1) +
      coord_cartesian(
            xlim = xlim_use,
            ylim = c(0, ymax1),
            expand = FALSE,
            clip = "off"
      ) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

# RIGHT SMOOTH CURVE: density of PCA axis 2 (Dim2)
#    -> Dim2 is shown on the Y axis, density on the X axis

dens_plot_y <- ggplot(ind_df, aes(y = Dim2)) +
      geom_density(linewidth = 0.8, color = "black", fill = NA, adjust = 1) +
      coord_cartesian(
            ylim = ylim_use,
            xlim = c(0, ymax2),
            expand = FALSE,
            clip = "off"
      ) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

#PCA plot

txt_ticks <- 14   # <-- taille chiffres axes
txt_axis  <- 16   # <-- taille noms des axes

theme_axes <- theme(
      axis.text  = element_text(size = txt_ticks, color = "black"),
      axis.title = element_text(size = txt_axis,  color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.18, "cm")
)

PCA_plot <- ggplot() +
      geom_point(
            data = ind_df,
            aes(x = Dim1, y = Dim2),
            color = ind_df$color,
            size = 2.5,
            shape = 16,
            stroke = 0,
            alpha = 1
      ) +
      geom_density_2d(
            data = ind_df,
            aes(x = Dim1, y = Dim2),
            contour_var = "ndensity",
            breaks = c(0.25, 0.5, 0.75, 0.95),
            color = "black",
            linewidth = 1
      ) +
      geom_segment(
            data = var_df,
            aes(x = 0, y = 0, xend = Dim1 * r, yend = Dim2 * r),
            arrow = arrow(length = unit(0.2, "cm")),
            linewidth = 0.7,
            color = "black"
      ) +
      geom_text(
            data = var_df,
            aes(x = Dim1 * r * 1.15, y = Dim2 * r * 1.1, label = varname),
            color = "black",
            size = 4.5
      ) +
      xlab("PCA1 (73.05%)") +
      ylab("PCA2 (23.03%)") +
      #coord_equal(xlim = xlim_use, ylim = ylim_use) +
      coord_cartesian(xlim = xlim_use, ylim = ylim_use)+

      scale_x_continuous(breaks = seq(-6, 5, 1)) +
      scale_y_continuous(breaks = seq(-6, 5, 1)) +
      theme_classic() +
      theme(
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      theme_axes

# 4) INSERT THE TWO DENSITY CURVES INTO THE PCA GTABLE
#    so that they are perfectly aligned with the PCA panel

g_pca    <- ggplotGrob(PCA_plot)
g_densx  <- ggplotGrob(dens_plot_x)
g_densy  <- ggplotGrob(dens_plot_y)

# Position of the PCA panel
p_pca <- g_pca$layout[g_pca$layout$name == "panel", ]

# Extract density panels
panel_densx <- g_densx$grobs[[which(g_densx$layout$name == "panel")]]
panel_densy <- g_densy$grobs[[which(g_densy$layout$name == "panel")]]

# Size of the added panels (can be adjusted)
h_top   <- unit(1.0, "in")
w_right <- unit(1.0, "in")

# add a row above the PCA panel
g_pca <- gtable_add_rows(g_pca, heights = h_top, pos = p_pca$t - 1)
p_pca2 <- g_pca$layout[g_pca$layout$name == "panel", ]

g_pca <- gtable_add_grob(
      g_pca, panel_densx,
      t = p_pca2$t - 1, b = p_pca2$t - 1,
      l = p_pca2$l,     r = p_pca2$r,
      clip = "off"   

# Border around top density panel

g_pca <- gtable_add_grob(
      g_pca,
      rectGrob(gp = gpar(col = "black", fill = NA, lwd = 1)),
      t = p_pca2$t - 1, b = p_pca2$t - 1,
      l = p_pca2$l,     r = p_pca2$r,
      clip = "off"
)

# add a column to the right of the PCA panel
g_pca <- gtable_add_cols(g_pca, widths = w_right, pos = p_pca2$r)
p_pca3 <- g_pca$layout[g_pca$layout$name == "panel", ]

g_pca <- gtable_add_grob(
      g_pca, panel_densy,
      t = p_pca3$t, b = p_pca3$b,
      l = p_pca3$r + 1, r = p_pca3$r + 1,
      clip = "off"  
)

# Border around right density pane
g_pca <- gtable_add_grob(
      g_pca,
      rectGrob(gp = gpar(col = "black", fill = NA, lwd = 1)),
      t = p_pca3$t, b = p_pca3$b,
      l = p_pca3$r + 1, r = p_pca3$r + 1,
      clip = "off"
)

grid.draw(g_pca)



####### Histograms 

# Histogram bin width along PCA1

binwidth_x <- 1

breaks_x <- seq(xlim_use[1], xlim_use[2], by = binwidth_x)


ind_x <- ind_df %>%
      mutate(bin = cut(Dim1, breaks = breaks_x, include.lowest = TRUE, right = FALSE)) %>%
      filter(!is.na(bin))

summary(ind_x$Dim2)
quantile(ind_df$Dim2, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)

# Histogram coloured using the species colour closest to the 5th percentile of PCA2
hist_x_5<- ind_x %>%
      group_by(bin) %>%
      summarise(
            count = n(),
            #fill  = color[which.max(Dim2)],
            #fill  = color[which.min(abs(Dim2))],
            #fill  = color[which.min(Dim2)],
            
            fill = color[which.min(abs(Dim2 - (-2.06)))], #qauntile 5%
            #fill = color[which.min(abs(Dim2 - 2.14))], # quantile 95%
            #fill = color[which.min(abs(Dim2 - 0.027))], #median
            
            
            
            
            .groups = "drop"
      ) %>%
      mutate(
            bin_id  = as.integer(bin),
            x_left  = breaks_x[bin_id],
            x_right = breaks_x[bin_id + 1],
            x_mid   = (x_left + x_right) / 2
      )

hist_plot_x_5 <- ggplot(hist_x_5, aes(x = x_mid, y = count)) +
      geom_col(aes(fill = fill), width = binwidth_x, color = "black", linewidth = 0.25) +
      scale_fill_identity() +
      coord_cartesian(xlim = xlim_use, expand = FALSE) +
      theme_classic() +
      theme(plot.margin = margin(0, 0, 0, 0))

hist_plot_x_5




# Histogram coloured using the species colour closest to the median of PCA2

hist_x_50<- ind_x %>%
      group_by(bin) %>%
      summarise(
            count = n(),
            #fill  = color[which.max(Dim2)],
            #fill  = color[which.min(abs(Dim2))],
            #fill  = color[which.min(Dim2)],
            
            #fill = color[which.min(abs(Dim2 - (-2.06)))], #qauntile 5%
            #fill = color[which.min(abs(Dim2 - 2.14))], # quantile 95%
            fill = color[which.min(abs(Dim2 - 0.027))], #median
            
            
            
            
            .groups = "drop"
      ) %>%
      mutate(
            bin_id  = as.integer(bin),
            x_left  = breaks_x[bin_id],
            x_right = breaks_x[bin_id + 1],
            x_mid   = (x_left + x_right) / 2
      )

hist_plot_x_50 <- ggplot(hist_x_50, aes(x = x_mid, y = count)) +
      geom_col(aes(fill = fill), width = binwidth_x, color = "black", linewidth = 0.25) +
      scale_fill_identity() +
      coord_cartesian(xlim = xlim_use, expand = FALSE) +
      theme_classic() +
      theme(plot.margin = margin(0, 0, 0, 0))

hist_plot_x_50


# Histogram coloured using the species colour closest to the 95th percentile of PCA2

hist_x_95<- ind_x %>%
      group_by(bin) %>%
      summarise(
            count = n(),
            #fill  = color[which.max(Dim2)],
            #fill  = color[which.min(abs(Dim2))],
            #fill  = color[which.min(Dim2)],
            
            #fill = color[which.min(abs(Dim2 - (-2.06)))], #qauntile 5%
            fill = color[which.min(abs(Dim2 - 2.14))], # quantile 95%
            #fill = color[which.min(abs(Dim2 - 0.027))], #median
            
            
            
            
            .groups = "drop"
      ) %>%
      mutate(
            bin_id  = as.integer(bin),
            x_left  = breaks_x[bin_id],
            x_right = breaks_x[bin_id + 1],
            x_mid   = (x_left + x_right) / 2
      )

hist_plot_x_95 <- ggplot(hist_x_95, aes(x = x_mid, y = count)) +
      geom_col(aes(fill = fill), width = binwidth_x, color = "black", linewidth = 0.25) +
      scale_fill_identity() +
      coord_cartesian(xlim = xlim_use, expand = FALSE) +
      theme_classic() +
      theme(plot.margin = margin(0, 0, 0, 0))

hist_plot_x_95


# Histogram coloured using the species with the maximum PCA2 value in each PCA1 bin

hist_x_max<- ind_x %>%
      group_by(bin) %>%
      summarise(
            count = n(),
            fill  = color[which.max(Dim2)],
            #fill  = color[which.min(abs(Dim2))],
            #fill  = color[which.min(Dim2)],
            
            #fill = color[which.min(abs(Dim2 - (-2.06)))], #qauntile 5%
            #fill = color[which.min(abs(Dim2 - 2.14))], # quantile 95%
            #fill = color[which.min(abs(Dim2 - 0.027))], #median
            
            
            
            
            .groups = "drop"
      ) %>%
      mutate(
            bin_id  = as.integer(bin),
            x_left  = breaks_x[bin_id],
            x_right = breaks_x[bin_id + 1],
            x_mid   = (x_left + x_right) / 2
      )

hist_plot_x_max <- ggplot(hist_x_max, aes(x = x_mid, y = count)) +
      geom_col(aes(fill = fill), width = binwidth_x, color = "black", linewidth = 0.25) +
      scale_fill_identity() +
      coord_cartesian(xlim = xlim_use, expand = FALSE) +
      theme_classic() +
      theme(plot.margin = margin(0, 0, 0, 0))

hist_plot_x_max

(hist_plot_x_5 / hist_plot_x_50 / hist_plot_x_95 / hist_plot_x_max) +
      plot_layout(ncol = 1, heights = c(1,1,1,1))


# breaks 
x_breaks_int <- -6:5

p1 <- hist_plot_x_5 +
      scale_x_continuous(
            limits = xlim_use,
            breaks = x_breaks_int,
            labels = function(x) sprintf("%d", as.integer(x))
      ) +
      xlab("PCA1 (73.05%)") +
      ylab("Number of species")

p2 <- hist_plot_x_50 +
      ylab(NULL) +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())

p3 <- hist_plot_x_95 +
      ylab(NULL) +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())

p4 <- hist_plot_x_max +
      ylab(NULL) +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())

# Rebuild the final histogram panel
hist_final <- (p4 / p3 / p2 / p1) +
      plot_layout(ncol = 1)

# Styling (frame + text sizes)
hist_final <- hist_final &
      theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
      ) &
      theme_axes

# Make the Y axis slightly more discreet
hist_final <- hist_final & theme(
      axis.text.y = element_text(size = 14)
)


hist_final

# ### save the two plots
# 
# side_pca <- 8            
# hist_width_pct <- 90    
# 
# side_h <- side_pca
# side_hist <- side_pca * (hist_width_pct / 100)
# 
# pdf("Hist_final_90.pdf", width = side_hist, height = side_h, onefile = TRUE)
# print(hist_final)
# dev.off()
# 
# pdf("PCA_final.pdf", width = side_pca, height = side_pca, onefile = TRUE)
# grid.newpage()
# grid.draw(g_pca)
# dev.off()
# 

