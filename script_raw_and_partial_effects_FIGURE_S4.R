#### Title: Script Figure S4 for the study: "Beyond Gloger’s rule: multiple biogeographic drivers of red and dark colouration in ants"
### Authors: XX
## Analyses/script: XX
# Last update: 17/03/2025

library(visreg)  #v‘2.8.0’
library(ggplot2) #v‘3.5.2’

# =========================
# PCA1 
# =========================

# # Get PCA individual coordinates (axes 1 and 2)
# ind_coord <- factoextra::get_pca_ind(PCA_raw_colours)$coord[, 1:2]
# 
# ind_df <- data.frame(
#       Dim1 = ind_coord[, 1],
#       Dim2 = ind_coord[, 2],
#       color = hex_col
# )
# 
# # For each row in data_95_full, assign the nearest PCA color in (Dim1, Dim2) space
# data_95_full$hex_col_PCA12 <- apply(
#       data_95_full[, c("PCA1_mean_95rich", "PCA2_mean_95rich")],
#       1,
#       function(v) {
#             d <- (ind_df$Dim1 - v[1])^2 + (ind_df$Dim2 - v[2])^2
#             ind_df$color[which.min(d)]
#       }
# )
# 
# data_95_full$hex_col_PCA12

### data_95_full is provided directly as an .xlsx file 
### (.../Data/data_95_full.xlsx) at this stage to avoid having to rerun the entire pipeline,
### especially the environmental and climatic data extraction steps.


### Above code is already included "data_95_full"
data_95_full <- read_excel(".../Data/data_95_full.xlsx")

# ---- mean_UVB1_scaled ----
v_UVB1<- visreg(
      mod_sdm, # mod_sdm = model for PCA1_mean (darkness) row 1021 in the main script Beyond_gloger_ants_main_script.R
      "mean_UVB1_scaled",
      type = "conditional",
      plot = FALSE
)

v_UVB1$res$hex_col   <- data_95_full$hex_col_PCA12
v_UVB1$res$PCA1_mean_95rich <- data_95_full$PCA1_mean_95rich

size_point <- 2.5

# Partial residual plot 
p_partiel_UVB1 <- ggplot() +
      
      geom_point(
            data = v_UVB1$res,
            aes(x = mean_UVB1_scaled,
                y = visregRes,
                colour = hex_col),
            alpha = 1,
            size  = size_point
      ) +
      
      # Confidence band for the conditional effect
      geom_ribbon(
            data = v_UVB1$fit,
            aes(x = mean_UVB1_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "darkorchid4"
      ) +
      # Conditional effect line from the model
      geom_line(
            data  = v_UVB1$fit,
            aes(x = mean_UVB1_scaled, y = visregFit),
            colour = "darkorchid4",
            size   = 1.3
      ) +
      
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_UVB1_scaled",
            y = "Partial residuals",
            title = ""
      )

p_partiel_UVB1


# Raw response vs predictor 
p_brut_UVB1 <- ggplot() +
      
      geom_point(
            data = data_95_full,
            aes(x = mean_UVB1_scaled,
                y = PCA1_mean_95rich,
                colour = hex_col_PCA12),
            alpha = 1,
            size  = size_point
      ) +
      
      # Same confidence band as the partial effect
      geom_ribbon(
            data = v_UVB1$fit,
            aes(x = mean_UVB1_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "darkorchid4"
      ) +
      # Same conditional effect line as above
      geom_line(
            data  = v_UVB1$fit,
            aes(x = mean_UVB1_scaled, y = visregFit),
            colour = "darkorchid4",
            size   = 1.3
      ) +
      
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_UVB1_scaled",
            y = "PCA1_mean (brut)",
            title = ""
      )

p_brut_UVB1


# ---- mean_bio17_scaled ----
v_bio17 <- visreg(
      mod_sdm,
      "mean_bio17_scaled",
      type = "conditional",
      plot = FALSE
)

v_bio17$res$hex_col   <- data_95_full$hex_col_PCA12
v_bio17$res$PCA1_mean <- data_95_full$PCA1_mean_95rich

# Partial residual plot 
p_partiel_bio17 <- ggplot() +
      
      geom_point(
            data = v_bio17$res,
            aes(x = mean_bio17_scaled,
                y = visregRes,
                colour = hex_col),
            alpha = 1,
            size  = size_point
      ) +
      
      geom_ribbon(
            data = v_bio17$fit,
            aes(x = mean_bio17_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "cornflowerblue"
      ) +
      geom_line(
            data  = v_bio17$fit,
            aes(x = mean_bio17_scaled, y = visregFit),
            colour = "cornflowerblue",
            size   = 1.3
      ) +
      
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_bio17_scaled",
            y = "Partial residuals",
            title = ""
      )

p_partiel_bio17


# Raw response vs predictor 
p_brut_bio17 <- ggplot() +
      
      geom_point(
            data = data_95_full,
            aes(x = mean_bio17_scaled,
                y = PCA1_mean_95rich,
                colour = hex_col_PCA12),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_bio17$fit,
            aes(x = mean_bio17_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "cornflowerblue"
      ) +
      geom_line(
            data  = v_bio17$fit,
            aes(x = mean_bio17_scaled, y = visregFit),
            colour = "cornflowerblue",
            size   = 1.3
      ) +
      
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_bio17_scaled",
            y = "PCA1_mean (brut)",
            title = ""
      )

p_brut_bio17


# ---- mean_canopy_scaled ----
v_canopy <- visreg(
      mod_sdm,
      "mean_canopy_scaled",   
      type = "conditional",
      plot = FALSE
)

# Attach colors to partial residuals
v_canopy$res$hex_col   <- data_95_full$hex_col_PCA12
v_canopy$res$PCA1_mean <- data_95_full$PCA1_mean_95rich

# Partial residual plot 
p_partiel_canopy <- ggplot() +
      geom_point(
            data = v_canopy$res,
            aes(x = mean_canopy_scaled,
                y = visregRes,
                colour = hex_col),
            alpha = 1,
            size  = size_point
      ) +
      
      geom_ribbon(
            data = v_canopy$fit,
            aes(x = mean_canopy_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "#1B5E20"
      ) +
      geom_line(
            data  = v_canopy$fit,
            aes(x = mean_canopy_scaled, y = visregFit),
            colour = "#1B5E20",
            size   = 1.3
      ) +
      
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_canopy_scaled",
            y = "Partial residuals",
            title = ""
      )

p_partiel_canopy


# Raw response vs predictor 
p_brut_canopy <- ggplot() +
      geom_point(
            data = data_95_full,
            aes(x = mean_canopy_scaled,
                y = PCA1_mean_95rich,
                colour = hex_col_PCA12),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_canopy$fit,
            aes(x = mean_canopy_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "#1B5E20"
      ) +
      geom_line(
            data  = v_canopy$fit,
            aes(x = mean_canopy_scaled, y = visregFit),
            colour = "#1B5E20",
            size   = 1.3
      ) +
      
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_canopy_scaled",
            y = "PCA1_mean (brut)",
            title = ""
      )

p_brut_canopy


# Compute a common y-range across partial residual plots
all_res_partials <- c(
      v_UVB1$res$visregRes,
      v_bio17$res$visregRes,
      v_canopy$res$visregRes
)

ylims_partiel <- range(all_res_partials, na.rm = TRUE)
ylims_partiel

# Apply the same y-limits to each partial residual plot
p_partiel_UVB1 <- p_partiel_UVB1 +
      coord_cartesian(ylim = ylims_partiel)

p_partiel_bio17 <- p_partiel_bio17 +
      coord_cartesian(ylim = ylims_partiel)

p_partiel_canopy <- p_partiel_canopy +
      coord_cartesian(ylim = ylims_partiel)

# Combine raw plots (top) and partial residual plots (bottom)
p_all_PCA1 <- (p_brut_UVB1 | p_brut_canopy | p_brut_bio17 ) /
      (p_partiel_UVB1 | p_partiel_canopy | p_partiel_bio17 )

p_all_PCA1


# =========================
# PCA2 
# =========================

size_point <- 2.5


# ---- NMDS2_scaled ----
v_NMDS2 <- visreg(
      mod_sdm2, # mod_sdm2 = model for PCA2_mean (redness) row 1353 in the main script Beyond_gloger_ants_main_script.R
      "NMDS2_scaled",
      type = "conditional",
      plot = FALSE
)

# Attach colors to partial residuals
v_NMDS2$res$hex_col <- data_95_full$hex_col_PCA12
v_NMDS2$res$PCA2_mean_95rich <- data_95_full$PCA2_mean_95rich

# Partial residual plot (colored by PCA palette)
p_partiel_NMDS2 <- ggplot() +
      geom_point(
            data = v_NMDS2$res,
            aes(x = NMDS2_scaled,
                y = visregRes,
                colour = hex_col),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_NMDS2$fit,
            aes(x = NMDS2_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "lightgoldenrod2"
      ) +
      geom_line(
            data  = v_NMDS2$fit,
            aes(x = NMDS2_scaled, y = visregFit),
            colour = "lightgoldenrod2",
            size   = 1.3
      ) +
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "NMDS2_scaled",
            y = "Partial residuals",
            title = ""
      )

p_partiel_NMDS2

# Raw response vs predictor 
p_brut_NMDS2 <- ggplot() +
      geom_point(
            data = data_95_full,
            aes(x = NMDS2_scaled,
                y = PCA2_mean_95rich,
                colour = hex_col_PCA12),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_NMDS2$fit,
            aes(x = NMDS2_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "lightgoldenrod2"
      ) +
      geom_line(
            data  = v_NMDS2$fit,
            aes(x = NMDS2_scaled, y = visregFit),
            colour = "lightgoldenrod2",
            size   = 1.3
      ) +
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "NMDS2_scaled",
            y = "PCA2_mean (brut)",
            title = ""
      )

p_brut_NMDS2



# ---- mean_bio10_scaled ----
v_bio10 <- visreg(
      mod_sdm2,
      "mean_bio10_scaled",
      type = "conditional",
      plot = FALSE
)

# Attach colors to partial residuals
v_bio10$res$hex_col <- data_95_full$hex_col_PCA12
v_bio10$res$PCA2_mean_95rich <- data_95_full$PCA2_mean_95rich

# Partial residual plot 
p_partiel_bio10 <- ggplot() +
      geom_point(
            data = v_bio10$res,
            aes(x = mean_bio10_scaled,
                y = visregRes,
                colour = hex_col),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_bio10$fit,
            aes(x = mean_bio10_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "firebrick"
      ) +
      geom_line(
            data  = v_bio10$fit,
            aes(x = mean_bio10_scaled, y = visregFit),
            colour = "firebrick",
            size   = 1.3
      ) +
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_bio10_scaled",
            y = "Partial residuals",
            title = ""
      )

p_partiel_bio10

# Raw response vs predictor
p_brut_bio10 <- ggplot() +
      geom_point(
            data = data_95_full,
            aes(x = mean_bio10_scaled,
                y = PCA2_mean_95rich,
                colour = hex_col_PCA12),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_bio10$fit,
            aes(x = mean_bio10_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "firebrick"
      ) +
      geom_line(
            data  = v_bio10$fit,
            aes(x = mean_bio10_scaled, y = visregFit),
            colour = "firebrick",
            size   = 1.3
      ) +
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_bio10_scaled",
            y = "PCA2_mean (brut)",
            title = ""
      )

p_brut_bio10



# ---- mean_canopy_scaled ----
v_canopy2 <- visreg(
      mod_sdm2,
      "mean_canopy_scaled",
      type = "conditional",
      plot = FALSE
)

# Attach colors to partial residuals
v_canopy2$res$hex_col <- data_95_full$hex_col_PCA12
v_canopy2$res$PCA2_mean_95rich <- data_95_full$PCA2_mean_95rich

# Partial residual plot 
p_partiel_canopy2 <- ggplot() +
      geom_point(
            data = v_canopy2$res,
            aes(x = mean_canopy_scaled,
                y = visregRes,
                colour = hex_col),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_canopy2$fit,
            aes(x = mean_canopy_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "#1B5E20"
      ) +
      geom_line(
            data  = v_canopy2$fit,
            aes(x = mean_canopy_scaled, y = visregFit),
            colour = "#1B5E20",
            size   = 1.3
      ) +
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_canopy_scaled",
            y = "Partial residuals",
            title = ""
      )

p_partiel_canopy2

# Raw response vs predictor 
p_brut_canopy2 <- ggplot() +
      geom_point(
            data = data_95_full,
            aes(x = mean_canopy_scaled,
                y = PCA2_mean_95rich,
                colour = hex_col_PCA12),
            alpha = 1,
            size  = size_point
      ) +
      geom_ribbon(
            data = v_canopy2$fit,
            aes(x = mean_canopy_scaled, ymin = visregLwr, ymax = visregUpr),
            alpha = 0.4,
            fill  = "#1B5E20"
      ) +
      geom_line(
            data  = v_canopy2$fit,
            aes(x = mean_canopy_scaled, y = visregFit),
            colour = "#1B5E20",
            size   = 1.3
      ) +
      scale_color_identity() +
      theme_minimal() +
      labs(
            x = "mean_canopy_scaled",
            y = "PCA2_mean (brut)",
            title = ""
      )

p_brut_canopy2



# ---- use the same y-limits for all partial residual plots ----
all_res_partials <- c(
      v_NMDS2$res$visregRes,
      v_bio10$res$visregRes,
      v_canopy2$res$visregRes
)

ylims_partiel <- range(all_res_partials, na.rm = TRUE)
ylims_partiel

p_partiel_NMDS2 <- p_partiel_NMDS2 + coord_cartesian(ylim = ylims_partiel)
p_partiel_bio10 <- p_partiel_bio10 + coord_cartesian(ylim = ylims_partiel)
p_partiel_canopy2 <- p_partiel_canopy2 + coord_cartesian(ylim = ylims_partiel)

# Combine raw plots (top) and partial residual plots (bottom)
p_all_PCA2 <- (p_brut_NMDS2 | p_brut_bio10 | p_brut_canopy2) /
      (p_partiel_NMDS2 | p_partiel_bio10 | p_partiel_canopy2)

p_all_PCA2
