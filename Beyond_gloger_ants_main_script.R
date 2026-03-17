#### Title: Main analysis script for the study: "Beyond Gloger’s rule: multiple biogeographic drivers of red and dark colouration in ants"
### Authors: Tristan Klaftenberger, Brian L. Fisher, Eddie Pérochon, Cleo Bertelsmeier*, Sébastien Ollier*
## Analyses/script: Tristan Klaftenberger
# Last update: 17/03/2025

library(sf)             #v‘1.0.16’
library(dplyr)          #v'‘1.1.4
library(ade4)           #v‘1.7.22’
library(ggplot2)        #v‘3.5.2’
library(rnaturalearth)  #v‘1.0.1’
library(mgcv)           #v‘1.9.1’
library(terra)          #v‘1.7.78’
library(ggpubr)         #v‘0.6.0’
library(vegan)          #v‘2.6.6.1’
library(corrplot)       #v‘0.92’
library(sdmTMB)         #v‘0.7.4’
library(performance)    #v‘0.12.2’
library(gt)             #v‘1.0.0’
#remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE) #‘0.0.5’
library(sdmTMBextra)    #v‘0.0.5’
library(patchwork)      #v‘1.3.0’
library(spdep)          #v‘1.3.8’
library(ggrepel)        #v‘0.9.6’

### Load the occurances data from Kass et al 2022 
### Load the occurances data from Kass et al 2022
### Load the occurances data from Kass et al 2022
#DOI: 10.1126/sciadv.abp9908

KASS_occ_species <- read.csv(".../KASS_occ_species.csv") # "..." replace with your path
table(is.na(KASS_occ_species$lat))
table(is.na(KASS_occ_species$lon))

dim(KASS_occ_species)
length(unique(KASS_occ_species$valid_species_name)) # 14347 unique species

KASS_occ_species$valid_species_name <- gsub("\\.","_", KASS_occ_species$valid_species_name)

colnames(KASS_occ_species) <- c("lon","lat","valid_species_name","gabi_acc_number")


######  Ecoregions  ###### 
######  Ecoregions  ###### 
######  Ecoregions  ###### 

# Dinerstein et al. 2017; https://doi.org/10.1093/biosci/bix014

### dowload the ecoregions ###
### dowload the ecoregions ###
### dowload the ecoregions ###

Ecoreg <- st_read(".../Ecoregions2017.shp") # 847 ecoregions ; "..." replace with your path
st_crs(Ecoreg) # 4326

ants_occ_sf <- st_as_sf(KASS_occ_species, coords = c("lon","lat"), crs = 4326)
st_crs(ants_occ_sf) #4326

st_crs(Ecoreg) == st_crs(ants_occ_sf) # TRUE

Ecoreg <- st_make_valid(Ecoreg)

ants_occ_sf_intersect <- st_join(ants_occ_sf, Ecoreg, join = st_intersects)

table(is.na(ants_occ_sf_intersect$ECO_NAME)) #18116 NA --> probably the ones falling in the sea 

id <- which(is.na(ants_occ_sf_intersect$ECO_NAME))

ants_occ_sf_intersect <- ants_occ_sf_intersect[-id,]
dim(ants_occ_sf_intersect) #1461181      18
length(unique(ants_occ_sf_intersect$valid_species_name)) # 14317 unique species
length(unique(ants_occ_sf_intersect$ECO_NAME))

### Ecoregions considered as well sampled, dowoladed from Pérochon et al. 2026 https://doi.org/10.1073/pnas.2530826123 

load(".../sampling_bias_ecoregions.RDATA") # provided but from https://github.com/EddiePerochon/Ant_Sociality_Env/blob/main/outputs/sampling_bias_ecoregions.RDATA
l_eco_eddie_1 <- ecor50$ECO_NAME
rm(ecor50)
l_eco_eddie_1<- gsub(" ", "_", l_eco_eddie_1)

########## Dowload the colour data (thorax) ########## 
########## Dowload the colour data (thorax) ########## 
########## Dowload the colour data (thorax) ########## 

# We are sharing specimen colour data for 50 species (see Thorax_colour_50_species.xlsx) to allow verification of the replicability of mean species colour estimates and to enable users to compare these with their own colour extraction methods.
# We are providing the average colour for 10,400 species used for the analyses in the study
# Please contact us if you would like to work with the full raw dataset.


temp_env <- new.env()
load(".../data_clean_colour.RData", envir = temp_env) 
Thorax_colour <- temp_env[["Thorax_HSV_noNA"]]
rm(temp_env)

dim(Thorax_colour) #34331

table(is.na(Thorax_colour$H_hsv))
table(is.na(Thorax_colour$S_hsv))
table(is.na(Thorax_colour$V_hsv))

table(is.na(Thorax_colour$R))
table(is.na(Thorax_colour$G))
table(is.na(Thorax_colour$B))

table(is.na(Thorax_colour$L_luv))
table(is.na(Thorax_colour$U_luv))
table(is.na(Thorax_colour$V_luv))

summary(as.numeric(table(Thorax_colour$species))) # median = 2 specimens / species

freq_species <- table(Thorax_colour$species)

n_species_single <- sum(freq_species == 1) #3436
n_species_two <- sum(freq_species == 2) #2671
n_species_3 <- sum(freq_species == 3) #1340
n_species_total <- length(freq_species) #10400 species

percent_single_image <- (n_species_single / n_species_total) * 100 #33.0
percent_two_images <- (n_species_two / n_species_total) * 100 #25.68
percent_1_2 <- ((n_species_single+n_species_two)/n_species_total) * 100 #58.7
percent_1_2_3 <- ((n_species_single+n_species_two+n_species_3)/n_species_total) * 100 #71.6

####### mean colour per species ####### 
####### mean colour per species ####### 
####### mean colour per species ####### 

# Just below it's how we computed the mean_per_speices dataset but we are providing it as .xlsx --> mean_per_species.xlsx
mean_per_species <- Thorax_colour %>%
      group_by(species) %>%
      summarize(
            H_hsv_mean = mean(H_hsv, na.rm = TRUE),
            S_hsv_mean = mean(S_hsv, na.rm = TRUE),
            V_hsv_mean = mean(V_hsv, na.rm = TRUE),
            
            L_luv_mean = mean(L_luv, na.rm = TRUE),
            U_luv_mean = mean(U_luv, na.rm = TRUE),
            V_luv_mean = mean(V_luv, na.rm = TRUE),
            
            R_mean = mean(R, na.rm = TRUE),
            G_mean = mean(G, na.rm = TRUE),
            B_mean = mean(B, na.rm = TRUE),
            
            .groups = 'drop')

# adding the corresponding species mean colour to all the specimens
l_spe <- unique(mean_per_species$species) #10400 species

Thorax_colour[c("H_hsv_mean", "S_hsv_mean", "V_hsv_mean",
                "L_luv_mean", "U_luv_mean", "V_luv_mean",
                "R_mean", "G_mean", "B_mean")] <- NA
z=0

for (i in 1:length(l_spe)){
      
      id <- which(mean_per_species$species == l_spe[i])
      
      id2 <- which(Thorax_colour$species == l_spe[i])
      
      val_H_hsv <- mean_per_species [id, ][["H_hsv_mean"]]
      val_S_hsv <- mean_per_species [id, ][["S_hsv_mean"]]
      val_V_hsv <- mean_per_species [id, ][["V_hsv_mean"]]
      
      val_L_luv <- mean_per_species [id, ][["L_luv_mean"]]
      val_U_luv <- mean_per_species [id, ][["U_luv_mean"]]
      val_V_luv <- mean_per_species [id, ][["V_luv_mean"]]
      
      val_R <- mean_per_species [id, ][["R_mean"]]
      val_G <- mean_per_species [id, ][["G_mean"]]
      val_B <- mean_per_species [id, ][["B_mean"]]
      
      Thorax_colour[id2,"H_hsv_mean"] <- val_H_hsv
      Thorax_colour[id2,"S_hsv_mean"] <- val_S_hsv
      Thorax_colour[id2,"V_hsv_mean"] <- val_V_hsv
      
      Thorax_colour[id2,"L_luv_mean"] <- val_L_luv
      Thorax_colour[id2,"U_luv_mean"] <- val_U_luv
      Thorax_colour[id2,"V_luv_mean"] <- val_V_luv
      
      Thorax_colour[id2,"R_mean"] <- val_R
      Thorax_colour[id2,"G_mean"] <- val_G
      Thorax_colour[id2,"B_mean"] <- val_B
      
      z=z+1
      print(z)
      
}

ants_occ_sf_intersect[c("H_hsv_mean", "S_hsv_mean", "V_hsv_mean",
                        "L_luv_mean", "U_luv_mean", "V_luv_mean",
                        "R_mean", "G_mean", "B_mean")] <- NA

z=0

for (i in 1:length(l_spe)){
      
      id <- which(mean_per_species$species == l_spe[i])
      id2 <- which(ants_occ_sf_intersect$valid_species_name == l_spe[i])
      
      val_H_hsv <- mean_per_species [id, ][["H_hsv_mean"]]
      val_S_hsv <- mean_per_species [id, ][["S_hsv_mean"]]
      val_V_hsv <- mean_per_species [id, ][["V_hsv_mean"]]
      
      val_L_luv <- mean_per_species [id, ][["L_luv_mean"]]
      val_U_luv <- mean_per_species [id, ][["U_luv_mean"]]
      val_V_luv <- mean_per_species [id, ][["V_luv_mean"]]
      
      val_R <- mean_per_species [id, ][["R_mean"]]
      val_G <- mean_per_species [id, ][["G_mean"]]
      val_B <- mean_per_species [id, ][["B_mean"]]
      
      ants_occ_sf_intersect[id2,"H_hsv_mean"] <- val_H_hsv
      ants_occ_sf_intersect[id2,"S_hsv_mean"] <- val_S_hsv
      ants_occ_sf_intersect[id2,"V_hsv_mean"] <- val_V_hsv
      
      ants_occ_sf_intersect[id2,"L_luv_mean"] <- val_L_luv
      ants_occ_sf_intersect[id2,"U_luv_mean"] <- val_U_luv
      ants_occ_sf_intersect[id2,"V_luv_mean"] <- val_V_luv
      
      ants_occ_sf_intersect[id2,"R_mean"] <- val_R
      ants_occ_sf_intersect[id2,"G_mean"] <- val_G
      ants_occ_sf_intersect[id2,"B_mean"] <- val_B
      
      z=z+1
      print(z)
      
}


table(is.na(ants_occ_sf_intersect$V_hsv_mean))
length(unique(ants_occ_sf_intersect$valid_species_name))#14317

id <- which(is.na(ants_occ_sf_intersect$V_hsv_mean))

ants_occ_sf_intersect_noNA <- ants_occ_sf_intersect[-id,]
length(unique(ants_occ_sf_intersect_noNA$valid_species_name)) ### 9859

ants_occ_sf_intersect_noNA$ECO_NAME <- gsub(" ","_",ants_occ_sf_intersect_noNA$ECO_NAME)

df_result <- ants_occ_sf_intersect_noNA %>%
      group_by(ECO_NAME) %>%
      summarize(count_species = n_distinct(valid_species_name),
                n_speciemens = n())

summary(df_result$count_species) # median 64 species / ecoregion
summary(df_result$n_speciemens)  # median 318  specimens / ecoregion

plot(log(df_result$count_species),log(df_result$n_speciemens))
abline(a=0,b=1)
summary(lm(log(df_result$count_species)~log(df_result$n_speciemens))) #R²= 0.80

ants_all_eco_name <- ants_occ_sf_intersect_noNA

length(unique(ants_all_eco_name$valid_species_name)) # 9854 unique species
length(unique(ants_all_eco_name$ECO_NAME)) #753 ecoregions with ants

uni_eco <- unique(ants_all_eco_name$ECO_NAME)
table(is.na(uni_eco))

id_e1 <- which(ants_all_eco_name$ECO_NAME %in% l_eco_eddie_1)
ants_all_eco_name_e1 <- ants_all_eco_name[id_e1,]

length(unique(ants_all_eco_name_e1$ECO_NAME)) #633

setdiff(l_eco_eddie_1, unique(ants_all_eco_name_e1$ECO_NAME)) #28 ecoregions without ant coccurences

obj <- setdiff(l_eco_eddie_1, unique(ants_all_eco_name_e1$ECO_NAME))

##### PCA #####
##### PCA #####
##### PCA #####

# PCA across 10400 ant species

dim(mean_per_species) #10400 10

mean_per_species_PCA <- mean_per_species

mean_per_species_PCA$J_mean <- mean_per_species_PCA$R_mean + mean_per_species_PCA$G_mean #Yellow channel
hex_col <- with(mean_per_species_PCA,
                rgb(R_mean, G_mean, B_mean, maxColorValue = 255))


mean_per_species_PCA$hex_col <- hex_col

PCA_raw_colours <- dudi.pca(mean_per_species_PCA[,c(
      "S_hsv_mean",
      "V_hsv_mean",
      
      "L_luv_mean",
      "U_luv_mean",
      "V_luv_mean",
      
      "J_mean",
      "B_mean"
)])

summary(PCA_raw_colours)
s.arrow(PCA_raw_colours$c1)
s.corcircle(PCA_raw_colours$co)

# Extract the first two axes of the PCA

ind_coord <- factoextra::get_pca_ind(PCA_raw_colours)$coord[, 1:2]

ind_df <- data.frame(
      Dim1 = ind_coord[, 1],
      Dim2 = ind_coord[, 2],
      color = hex_col
)

var_coord <- factoextra::get_pca_var(PCA_raw_colours)$coord[, 1:2]

var_df <- data.frame(
      Dim1    = var_coord[, 1],
      Dim2    = var_coord[, 2],
      varname = rownames(var_coord)
)

rang_ind1 <- max(ind_df$Dim1) - min(ind_df$Dim1)
rang_ind2 <- max(ind_df$Dim2) - min(ind_df$Dim2)
rang_var1 <- max(var_df$Dim1) - min(var_df$Dim1)
rang_var2 <- max(var_df$Dim2) - min(var_df$Dim2)

r <- min(rang_ind1 / rang_var1, rang_ind2 / rang_var2) * 0.7

# PCA plot (Figure 2a) prior to visual enhancement in Inkscape ; but see script_PCA_density_histo_curves for Figure 2 a and b before enhancement in Inkscape 

PCA_plot <- ggplot() +
      geom_point(
            data = ind_df,
            aes(x = Dim1, y = Dim2),
            color = ind_df$color,
            size = 2.5,
            shape = 16,
            stroke = 0
      ) +
      geom_segment(
            data = var_df,
            aes(x = 0, y = 0, xend = Dim1 * r, yend = Dim2 * r),
            arrow = arrow(length = unit(0.2, "cm")),
            size = 0.7, color = "black"
      ) +
      geom_text(
            data = var_df,
            aes(x = Dim1 * r * 1.15, y = Dim2 * r * 1.1, label = varname),
            color = "black", size = 4.5
      ) +
      xlab("PCA1 (73.05%)") +
      ylab("PCA2 (23.03%)") +
      theme_classic() +
      theme(panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank()) +
      coord_equal(xlim = c(-6, 5), ylim = c(-6, 5)) +
      scale_x_continuous(breaks = seq(-6, 5, 1)) +  
      scale_y_continuous(breaks = seq(-6, 5, 1))

PCA_plot

PCA1 <- PCA_raw_colours$li[, 1]
PCA2 <- PCA_raw_colours$li[, 2]

mean_per_species_PCA$PCA1 <- PCA1
mean_per_species_PCA$PCA2 <- PCA2

head(mean_per_species_PCA)
head(ants_all_eco_name)

ants_with_PCA <- ants_all_eco_name %>%
      left_join(
            mean_per_species_PCA %>% select(species, PCA1, PCA2, hex_col), 
            by = c("valid_species_name" = "species")
      )

# Cpmpute the minimum, maximum and average values of PCA1 and PCA2 for each ecoregion

summary_by_ecoregion <- ants_with_PCA %>%
      group_by(ECO_NAME) %>%
      summarize(
            
            PCA1_min   = min(PCA1,   na.rm = TRUE),
            PCA1_max   = max(PCA1,   na.rm = TRUE),
            PCA1_mean  = mean(PCA1,  na.rm = TRUE),
            PCA1_med  = median(PCA1,  na.rm = TRUE),
            
            
            PCA2_min   = min(PCA2,   na.rm = TRUE),
            PCA2_max   = max(PCA2,   na.rm = TRUE),
            PCA2_mean  = mean(PCA2,  na.rm = TRUE),
            PCA2_med  = median(PCA2,  na.rm = TRUE),
            
            .groups = "drop"
      )


summary_e1 <- summary_by_ecoregion %>%
      filter(ECO_NAME %in% l_eco_eddie_1)


#### MAPS #####
#### MAPS #####
#### MAPS #####

# MAP PCA1

Ecoreg_m <- st_transform(Ecoreg, crs = "+proj=moll +datum=WGS84 +units=m") #Mollweide projection
st_crs(Ecoreg_m)

Ecoreg_m$ECO_NAME <-gsub(" ", "_", Ecoreg_m$ECO_NAME)

Ecoreg_m <- Ecoreg_m %>%
      left_join(
            df_result %>%
                  st_drop_geometry() %>%
                  select(ECO_NAME, count_species),
            by = "ECO_NAME"
      )


set.seed(666)

centroids_m <- st_centroid(Ecoreg_m)

coords_m <- st_coordinates(centroids_m)

df_centroids_m <- data.frame(
      ECO_NAME = Ecoreg_m$ECO_NAME,
      lon_centroid_m = coords_m[, "X"],
      lat_centroid_m = coords_m[, "Y"],
      stringsAsFactors = FALSE
)

Ecoreg_4_map <- Ecoreg_m %>%
      left_join(
            summary_e1 %>%
                  st_drop_geometry() %>%
                  select(ECO_NAME, PCA1_mean, PCA2_mean),
            by = "ECO_NAME"
      )

Ecoreg_4_map <- Ecoreg_4_map %>%
      left_join(
            df_centroids_m %>%
                  st_drop_geometry() %>%
                  select(ECO_NAME, lon_centroid_m, lat_centroid_m),
            by = "ECO_NAME"
      )

PCA1_range <- range(Ecoreg_4_map$PCA1_mean, na.rm = TRUE) #-2.809498  3.671034

tropic_lat <- c(-23.5, 0, 23.5)

tropic_lines_m <- st_sfc(
      lapply(tropic_lat, function(y) {
            st_linestring(matrix(c(-180, y, 180, y), ncol = 2, byrow = TRUE))
      }),
      crs = 4326
) |> st_transform(st_crs(Ecoreg_4_map))


ocean <- ne_download(scale = "medium", type = "ocean",
                     category = "physical", returnclass = "sf") |>
      st_transform(st_crs(Ecoreg_4_map))


# remoove 5% of the ecoregions with the lowest number of species

eco_e1 <- Ecoreg_4_map %>%
      filter(ECO_NAME %in% summary_e1$ECO_NAME)

q5_e1 <- quantile(eco_e1$count_species, probs = 0.05, na.rm = TRUE)
q5_e1

eco_e1$low_richness_5pct <- eco_e1$count_species <= q5_e1
sum(eco_e1$low_richness_5pct, na.rm = TRUE)

eco_to_na <- eco_e1$ECO_NAME[eco_e1$low_richness_5pct]

Ecoreg_4_map$PCA1_mean_95rich <- Ecoreg_4_map$PCA1_mean
Ecoreg_4_map$PCA2_mean_95rich <- Ecoreg_4_map$PCA2_mean

id_na <- Ecoreg_4_map$ECO_NAME %in% eco_to_na

Ecoreg_4_map$PCA1_mean_95rich[id_na] <- NA
Ecoreg_4_map$PCA2_mean_95rich[id_na] <- NA

table(is.na(Ecoreg_4_map$PCA1_mean_95rich)) #598 ecoregions OK
table(is.na(Ecoreg_4_map$PCA2_mean_95rich)) #598 ecoregions OK

PCA1_range_95 <- range(Ecoreg_4_map$PCA1_mean_95rich, na.rm = TRUE)
PCA1_range_95

lat_min_95 <- min(Ecoreg_4_map$lat_centroid_m[!is.na(Ecoreg_4_map$PCA1_mean_95rich)], na.rm = TRUE)
lat_max_95 <- max(Ecoreg_4_map$lat_centroid_m[!is.na(Ecoreg_4_map$PCA1_mean_95rich)], na.rm = TRUE)


map_plot_PCA1_mean_95 <- ggplot() +
      
      geom_sf(data = ocean, fill = "lightblue1", color = NA, alpha = 0.2) +
      
      geom_sf(data = tropic_lines_m, color = "gray30", linetype = "solid", linewidth = 1) +
      
      geom_sf(data = Ecoreg_4_map, aes(fill = PCA1_mean_95rich), color = NA) +
      
      # scale_fill_viridis_c(option = "B", direction = 1, name = "Colour diversity index (%)", na.value = "grey50") +
      
      scale_fill_viridis_c(
            option = "B", direction = -1,
            limits = PCA1_range_95,
            name = "PCA1_mean",
            na.value = "grey50"
      ) +
      
      
      scale_y_continuous(limits = c(lat_min_95, lat_max_95)) +
      
      theme_minimal(base_family = "Arial") +
      # labs(title = "Colour diversity of ant assemblages across ecoregions") +
      theme(
            text = element_text(family = "Arial"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.position = "left"
      )


data_e1_95 <- Ecoreg_4_map %>%
      dplyr::filter(ECO_NAME %in% summary_e1$ECO_NAME) %>%
      dplyr::filter(!is.na(PCA1_mean_95rich)) %>%
      sf::st_drop_geometry()

gam_lat_PCA1_mean_95 <- gam(
      PCA1_mean_95rich ~ s(lat_centroid_m, bs = "tp"),
      data = data_e1_95
)

summary(gam_lat_PCA1_mean_95)

pred_95 <- predict(gam_lat_PCA1_mean_95, newdata = data_e1_95, se.fit = TRUE)

data_e1_95 <- data_e1_95 |>
      dplyr::mutate(
            gam_pred_PCA1_mean_lat = pred_95$fit,
            gam_pred_PCA1_mean_lat_lower_ci = pred_95$fit - 1.96 * pred_95$se.fit,
            gam_pred_PCA1_mean_lat_upper_ci = pred_95$fit + 1.96 * pred_95$se.fit
      )

PCA1_range_95 <- range(Ecoreg_4_map$PCA1_mean_95rich, na.rm = TRUE)

PCA1_range_95 <- range(Ecoreg_4_map$PCA1_mean_95rich, na.rm = TRUE)

gam_plot_PCA1_mean_95 <- ggplot(data_e1_95, aes(x = lat_centroid_m, y = PCA1_mean_95rich)) +
      geom_vline(xintercept = latitudes_proj_y, linetype = "solid", color = "gray30", linewidth = 1) +
      geom_point(aes(fill = PCA1_mean_95rich), shape = 21, color = "black", stroke = 0.1, size = 2, alpha = 1) +
      geom_ribbon(aes(ymin = gam_pred_PCA1_mean_lat_lower_ci, ymax = gam_pred_PCA1_mean_lat_upper_ci),
                  fill = "grey30", alpha = 0.45) +
      geom_line(aes(y = gam_pred_PCA1_mean_lat), color = "black", linewidth = 1) +
      coord_flip() +
      scale_fill_viridis_c(option = "B", direction = -1, limits = PCA1_range_95, name = "") +
      
      scale_x_continuous(limits = c(lat_min_95, lat_max_95), 
                         labels = scales::label_comma(),
                         sec.axis = dup_axis(name = "Latitude"))+
      
      scale_y_continuous(name = "") +
      theme_minimal() +
      labs(x = NULL) +
      theme(
            axis.title.y.left  = element_blank(),
            axis.text.y.left   = element_blank(),
            axis.text.y.right  = element_text(size = 9),
            axis.title.y.right = element_text(size = 12),
            axis.text.x        = element_text(size = 9),
            plot.title         = element_text(hjust = 0.5),
            aspect.ratio       = 1.2,
            legend.position    = "none"
      )


PCA1_mean_map_lat_95 <-  map_plot_PCA1_mean_95 | gam_plot_PCA1_mean_95
PCA1_mean_map_lat_95

# MAP PCA2

PCA2_range_95 <- range(Ecoreg_4_map$PCA2_mean_95rich, na.rm = TRUE)

lat_min_95_PCA2 <- min(Ecoreg_4_map$lat_centroid_m[!is.na(Ecoreg_4_map$PCA2_mean_95rich)], na.rm = TRUE)
lat_max_95_PCA2 <- max(Ecoreg_4_map$lat_centroid_m[!is.na(Ecoreg_4_map$PCA2_mean_95rich)], na.rm = TRUE)

map_plot_PCA2_mean_95 <- ggplot() +
      geom_sf(data = ocean, fill = "lightblue1", color = NA, alpha = 0.2) +
      geom_sf(data = tropic_lines_m, color = "gray30", linetype = "solid", linewidth = 1) +
      geom_sf(data = Ecoreg_4_map, aes(fill = PCA2_mean_95rich), color = NA) +
      scale_fill_viridis_c(
            option = "B", direction = -1,
            limits = PCA2_range_95,
            name = "PCA2_mean",
            na.value = "grey50"
      ) +
      scale_y_continuous(limits = c(lat_min_95_PCA2, lat_max_95_PCA2)) +
      theme_minimal(base_family = "Arial") +
      theme(
            text = element_text(family = "Arial"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.position = "left"
      )


data_e1_95_PCA2 <- Ecoreg_4_map %>%
      dplyr::filter(ECO_NAME %in% summary_e1$ECO_NAME) %>%
      dplyr::filter(!is.na(PCA2_mean_95rich)) %>%
      sf::st_drop_geometry()

gam_lat_PCA2_mean_95 <- gam(
      PCA2_mean_95rich ~ s(lat_centroid_m, bs = "tp"),
      data = data_e1_95_PCA2
)
summary(gam_lat_PCA2_mean_95)

pred_95_PCA2 <- predict(gam_lat_PCA2_mean_95, newdata = data_e1_95_PCA2, se.fit = TRUE)

data_e1_95_PCA2 <- data_e1_95_PCA2 |>
      dplyr::mutate(
            gam_pred_PCA2_mean_lat = pred_95_PCA2$fit,
            gam_pred_PCA2_mean_lat_lower_ci = pred_95_PCA2$fit - 1.96 * pred_95_PCA2$se.fit,
            gam_pred_PCA2_mean_lat_upper_ci = pred_95_PCA2$fit + 1.96 * pred_95_PCA2$se.fit
      )

gam_plot_PCA2_mean_95 <- ggplot(data_e1_95_PCA2, aes(x = lat_centroid_m, y = PCA2_mean_95rich)) +
      geom_vline(xintercept = latitudes_proj_y, linetype = "solid", color = "gray30", linewidth = 1) +
      geom_point(aes(fill = PCA2_mean_95rich), shape = 21, color = "black", stroke = 0.1, size = 2, alpha = 1) +
      geom_ribbon(aes(ymin = gam_pred_PCA2_mean_lat_lower_ci, ymax = gam_pred_PCA2_mean_lat_upper_ci),
                  fill = "grey30", alpha = 0.45) +
      geom_line(aes(y = gam_pred_PCA2_mean_lat), color = "black", linewidth = 1) +
      coord_flip() +
      scale_fill_viridis_c(option = "B", direction = -1, limits = PCA2_range_95, name = "") +
      scale_x_continuous(limits = c(lat_min_95_PCA2, lat_max_95_PCA2), labels = scales::label_comma(),
                         sec.axis = dup_axis(name = "Latitude")) +
      scale_y_continuous(name = "") +
      theme_minimal() +
      labs(x = NULL) +
      theme(
            axis.title.y.left  = element_blank(),
            axis.text.y.left   = element_blank(),
            axis.text.y.right  = element_text(size = 9),
            axis.title.y.right = element_text(size = 12),
            axis.text.x        = element_text(size = 9),
            plot.title         = element_text(hjust = 0.5),
            aspect.ratio       = 1.2,
            legend.position    = "none"
      )

PCA2_mean_map_lat_95 <- map_plot_PCA2_mean_95 | gam_plot_PCA2_mean_95


### Final plot Figure 3 before visual enhancement in Inkscape
# Figure 3: Geographic patterns of mean ant colour variation across ecoregions.

final_plot_maps <- PCA1_mean_map_lat_95 / PCA2_mean_map_lat_95

final_plot_maps

##### Climates variables extraction ##### 
##### Climates variables extraction ##### 
##### Climates variables extraction ##### 

noms_fichiers <- paste0(".../Data/Rasters/wc2.1_30s_bio/wc2.1_30s_bio_", 1:19, ".tif")

bioclim_stack <- rast(noms_fichiers)

noms_fichiers_UV <- paste0(".../Data/Rasters/UV/UVB", 1:6, ".tif")

UV_stack <- terra::rast(noms_fichiers_UV)

ele_rast <- rast(".../Data/Raster/wc2.1_30s_elev/wc2.1_30s_elev.tif")

canopy_height <- rast("...Data/Raster/Canopy/canopy_height_30sc_noseas.tif")

l <- unique(data_e1$ECO_NAME)

dfs_list <- list()

j=0
set.seed(666)

Ecoreg$ECO_NAME <- gsub(" ", "_", Ecoreg$ECO_NAME)

list_eco_failed <- c()

for (eco_name in l) {
      
      tryCatch({
            
            rd_pts_current <- st_sample(Ecoreg[which(Ecoreg$ECO_NAME == eco_name),], size = 5000, type = "random")
            rd_pts_current_sp <- st_as_sf(rd_pts_current)
            
            bioclims_val <- terra::extract(bioclim_stack, rd_pts_current_sp)
            uvb_val <- terra::extract(UV_stack, rd_pts_current_sp)
            ele_val <- terra::extract(ele_rast,rd_pts_current_sp)
            canopy_val <- terra::extract(canopy_height,rd_pts_current_sp)
            
            val <- cbind(bioclims_val, uvb_val, ele_val, canopy_val)
            val <-  val[, !grepl("ID", colnames(val))]  # remoove ID coloumns
            
            means <- as.data.frame(t(sapply(val, mean, na.rm = TRUE)))
            colnames(means)[c(1:27)] <- c("mean_bio1",
                                          "mean_bio2",
                                          "mean_bio3",
                                          "mean_bio4",
                                          "mean_bio5",
                                          "mean_bio6",
                                          "mean_bio7",
                                          "mean_bio8",
                                          "mean_bio9",
                                          "mean_bio10",
                                          "mean_bio11",
                                          "mean_bio12",
                                          "mean_bio13",
                                          "mean_bio14",
                                          "mean_bio15",
                                          "mean_bio16",
                                          "mean_bio17",
                                          "mean_bio18",
                                          "mean_bio19",
                                          
                                          "mean_UVB1",
                                          "mean_UVB2",
                                          "mean_UVB3",
                                          "mean_UVB4",
                                          "mean_UVB5",
                                          "mean_UVB6",
                                          
                                          "mean_elev",
                                          "mean_canopy")
            
            medians <- as.data.frame(t(sapply(val, median, na.rm = TRUE)))
            colnames(medians)[c(1:27)] <- c("med_bio1",
                                            "med_bio2",
                                            "med_bio3",
                                            "med_bio4",
                                            "med_bio5",
                                            "med_bio6",
                                            "med_bio7",
                                            "med_bio8",
                                            "med_bio9",
                                            "med_bio10",
                                            "med_bio11",
                                            "med_bio12",
                                            "med_bio13",
                                            "med_bio14",
                                            "med_bio15",
                                            "med_bio16",
                                            "med_bio17",
                                            "med_bio18",
                                            "med_bio19",
                                            
                                            "med_UVB1",
                                            "med_UVB2",
                                            "med_UVB3",
                                            "med_UVB4",
                                            "med_UVB5",
                                            "med_UVB6",
                                            
                                            "med_elev",
                                            "med_canopy")
            
            
            
            eco_df <- cbind(means,medians)
            
            
            # add eco name 
            eco_df$ECO_NAME <- eco_name
            #add to the list
            dfs_list[[eco_name]] <- eco_df
            
            j=j+1
            print(j)
            
      }, error = function(e) {
            # Leave the region empty in case of error
            message(paste("Error for ecoreg:", eco_name))
            list_eco_failed <<- append(list_eco_failed,eco_name)
      })
}

length(list_eco_failed) # 0
dim(do.call(rbind, dfs_list)) # 633  55

values_5000_633eco <- do.call(rbind, dfs_list)

id_na <- which(is.na(Ecoreg_4_map$PCA1_mean_95rich))
data_95 <- st_drop_geometry(Ecoreg_4_map[-id_na,]) # the data set to use for the analyses

data_95_full <- data_95 %>%
      left_join(values_5000_633eco, by = "ECO_NAME")

#### phylogenetic dissimilarity ####
#### phylogenetic dissimilarity ####
#### phylogenetic dissimilarity ####

diss <- readRDS(".../Data/phy_diss.RDS") # Provided but from : https://github.com/EddiePerochon/Ant_Sociality_Env/blob/main/outputs/phy_diss.RDS
D_glob <- diss$glob 
class(D_glob) #"matrix" "array" 
head(D_glob)

colnames(D_glob)
rownames(D_glob)
rownames(D_glob) <- colnames(D_glob)

na_per_eco <- rowSums(is.na(D_glob)) # 2 per ecoregion, exept 2 that have 725 
sort(na_per_eco, decreasing = TRUE)[1:10]  #ID the two problematic ecoregions
bad_eco <- c("e335", "e522") # no phylo distances with any ecoregions..

keep <- !(rownames(D_glob) %in% bad_eco)
D_glob_clean <- D_glob[keep, keep]
any(is.na(D_glob_clean)) #FALSE

D_dist <- as.dist(D_glob_clean) #necessary for the NMDS

nmds <- metaMDS(
      D_dist,
      k = 5,
      trymax = 200
)

nmds$stress # 0.083 really good (k=5)

nmds_scores <- as.data.frame(vegan::scores(nmds, display = "sites"))

nmds_scores$OBJECTID <- rownames(nmds_scores)

head(nmds_scores)
head(data_95_full)

table(is.na(nmds_scores$NMDS1))

data_95_full$OBJECTID <- paste0("e",data_95_full$OBJECTID)

data_95_full <- data_95_full %>%
      left_join(nmds_scores, by = "OBJECTID")

head(data_95_full)

table(is.na(data_95_full$NMDS1)) #586    12 
table(is.na(data_95_full$NMDS2)) #586    12 

id_na_nmds <- which(is.na(data_95_full$NMDS1))
data_95_full[id_na_nmds,c("ECO_NAME","SHAPE_AREA")]
dim(data_95_full) #598  81

data_95_full <- data_95_full[-id_na_nmds,]
dim(data_95_full) #586  81 --> final dataset 

##### Modelisation #####
##### Modelisation #####
##### Modelisation #####

# Variables correlations 2 by 2

data_95_full$count_species_log <- log(data_95_full$count_species)
hist(data_95_full$count_species)
hist(data_95_full$count_species_log)

data_95_full$SHAPE_AREA_log <- log(data_95_full$SHAPE_AREA)
hist(data_95_full$SHAPE_AREA)
hist(data_95_full$SHAPE_AREA_log)


vars_to_scale <- c(
      "count_species_log",
      "SHAPE_AREA_log",
      "NMDS1", "NMDS2", "NMDS3", "NMDS4", "NMDS5",
      "mean_bio10", "mean_bio16", "mean_bio17",
      "mean_UVB1",
      "mean_elev",
      "mean_canopy"
)

data_95_full[paste0(vars_to_scale, "_scaled")] <- scale(data_95_full[vars_to_scale])

vars <- data_95_full[, c(
      "count_species_log_scaled",
      "SHAPE_AREA_log_scaled",
      
      "NMDS1_scaled"  ,
      "NMDS2_scaled"   ,
      
      "mean_bio10_scaled",
      "mean_bio16_scaled",
      "mean_bio17_scaled",
      "mean_UVB1_scaled",
      "mean_elev_scaled",
      "mean_canopy_scaled"
      
      
)]


var_labels <- c(
      "Species richness (log)",
      "Ecoregion area (log)",
      
      "NMDS1",
      "BIO10",
      "BIO16",
      "BIO17",
      "UVB1",
      "Elevation",
      "Canopy height"
)


cor_mat <- cor(vars, use = "complete.obs")

colnames(cor_mat) <- var_labels
rownames(cor_mat) <- var_labels

# Figure S1: Correlation structure among explanatory variables

corrplot(
      cor_mat,
      method = "color",
      type = "upper",
      addCoef.col = "black",
      tl.col = "black",
      tl.srt = 45,
      diag = FALSE
)

### PCA1 Pale-Dark ###
### PCA1 Pale-Dark ###
### PCA1 Pale-Dark ###

data_95_full$X_km <- data_95_full$lon_centroid_m / 1000 #meters to kilometers
data_95_full$Y_km <- data_95_full$lat_centroid_m / 1000

### data_95_full is provided directly as an .xlsx file 
### (.../Data/data_95_full.xlsx) at this stage to avoid having to rerun the entire pipeline,
### especially the environmental and climatic data extraction steps.

#  AIC selection for the best mesh

cutoffs_km <- seq(100, 2000, by = 10)

results <- list()

for (c_km in cutoffs_km) {
      cat("\n====================\nTesting cutoff =", c_km, "km\n====================\n")
      
      # Mesh
      mesh_c <- make_mesh(
            data_95_full,
            xy_cols = c("X_km", "Y_km"),
            cutoff = c_km
      )
      
      # Fit model
      mod_c <- try(
            sdmTMB(
                  PCA1_mean_95rich ~ 
                        count_species_log_scaled +
                        NMDS2_scaled +
                        SHAPE_AREA_log_scaled +
                        mean_bio10_scaled +
                        mean_bio16_scaled +
                        mean_bio17_scaled +
                        mean_UVB1_scaled +
                        mean_elev_scaled +
                        mean_canopy_scaled,
                  
                  data   = data_95_full,
                  mesh   = mesh_c,
                  spatial = "on",
                  family  = student(link = "identity")
            ),
            silent = TRUE
      )
      
      if (inherits(mod_c, "try-error")) {
            cat("!! cutoff", c_km, "error, AIC = NA\n")
            results[[as.character(c_km)]] <- list(
                  mesh  = mesh_c,
                  model = NA,
                  aic   = NA_real_
            )
      } else {
            results[[as.character(c_km)]] <- list(
                  mesh  = mesh_c,
                  model = mod_c,
                  aic   = AIC(mod_c)
            )
      }
}

res_cutoff <- data.frame(
      cutoff_km = cutoffs_km,
      AIC       = sapply(results, function(x) x$aic)
)

res_cutoff_ok <- na.omit(res_cutoff)

res_cutoff_ok
plot(res_cutoff_ok$cutoff_km, res_cutoff_ok$AIC, type = "b")

best_mesh <- res_cutoff_ok[order(res_cutoff_ok$AIC), ][1:10, ]
best_mesh



mesh <- make_mesh(
      data_95_full,
      type = "cutoff",
      xy_cols = c("X_km", "Y_km"),
      cutoff = 500 # <-- best_mesh
)

mod_sdm <- sdmTMB(
      
      PCA1_mean_95rich ~ 
            
            count_species_log_scaled +
            NMDS2_scaled +
            SHAPE_AREA_log_scaled +
            mean_bio10_scaled +
            mean_bio16_scaled +
            mean_bio17_scaled +
            mean_UVB1_scaled +
            mean_elev_scaled +
            mean_canopy_scaled,
      
      data = data_95_full,
      mesh = mesh,
      spatial = "on",
      family  = student(link = "identity"))

summary(mod_sdm)
sanity(mod_sdm)
r2(mod_sdm)

# p-values does not make sense with this type of model
# we look at the CIs instead

coefs_90   <- tidy(mod_sdm, effects = "fixed", conf.int = TRUE, conf.level = 0.90)  %>% 
      rename(conf.low.90 = conf.low, conf.high.90 = conf.high)


coefs_95   <- tidy(mod_sdm, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
      rename(conf.low.95 = conf.low, conf.high.95 = conf.high)

coefs_99   <- tidy(mod_sdm, effects = "fixed", conf.int = TRUE, conf.level = 0.99) %>% 
      rename(conf.low.99 = conf.low, conf.high.99 = conf.high)

coefs_999  <- tidy(mod_sdm, effects = "fixed", conf.int = TRUE, conf.level = 0.999) %>% 
      rename(conf.low.999 = conf.low, conf.high.999 = conf.high)

coefs_all <- coefs_95 %>%
      select(term, estimate, std.error) %>%
      left_join(coefs_90  %>% select(term, conf.low.90,  conf.high.90),  by = "term") %>%
      left_join(coefs_95  %>% select(term, conf.low.95,  conf.high.95),  by = "term") %>%
      left_join(coefs_99  %>% select(term, conf.low.99,  conf.high.99),  by = "term") %>%
      left_join(coefs_999 %>% select(term, conf.low.999, conf.high.999), by = "term") %>%
      mutate(
            significance = case_when(
                  conf.low.999 > 0 | conf.high.999 < 0 ~ "***",
                  conf.low.99  > 0 | conf.high.99  < 0 ~ "**",
                  conf.low.95  > 0 | conf.high.95  < 0 ~ "*",
                  conf.low.90  > 0 | conf.high.90  < 0 ~ ".",
                  TRUE ~ ""
            )
      )

coefs_all[,c("term","estimate","significance")]
coefs_all

var_labels <- c(
      "count_species_log_scaled"   = "Species richness (log)",
      "NMDS2_scaled"           = "Phylogenetic similarity",
      "SHAPE_AREA_log_scaled"  = "Ecoregion area (log)",
      "mean_bio10_scaled"      = "Temperature - warmest quarter (BIO10)",
      "mean_bio16_scaled"      = "Precipitation - wettest quarter (BIO16)",
      "mean_bio17_scaled"      = "Precipitation - driest quarter (BIO17)",
      "mean_UVB1_scaled"       = "UV-B radiation (UVB1)",
      "mean_elev_scaled"       = "Elevation",
      "mean_canopy_scaled"     = "Canopy height"
)

tab_df_PCA1 <- coefs_all %>%
      transmute(
            Predictor = recode(term, !!!var_labels, "(Intercept)" = "(Intercept)"),
            Estimate = estimate,
            SE = std.error,
            ci90_txt   = sprintf("[%.3f, %.3f]", conf.low.90,  conf.high.90),
            ci95_txt   = sprintf("[%.3f, %.3f]", conf.low.95,  conf.high.95),
            ci99_txt   = sprintf("[%.3f, %.3f]", conf.low.99,  conf.high.99),
            ci999_txt  = sprintf("[%.3f, %.3f]", conf.low.999, conf.high.999),
            Sig = significance,
            max_ci = case_when(
                  (conf.low.999 > 0 & conf.high.999 > 0) | (conf.low.999 < 0 & conf.high.999 < 0) ~ "ci999_txt",
                  (conf.low.99  > 0 & conf.high.99  > 0) | (conf.low.99  < 0 & conf.high.99  < 0) ~ "ci99_txt",
                  (conf.low.95  > 0 & conf.high.95  > 0) | (conf.low.95  < 0 & conf.high.95  < 0) ~ "ci95_txt",
                  (conf.low.90  > 0 & conf.high.90  > 0) | (conf.low.90  < 0 & conf.high.90  < 0) ~ "ci90_txt",
                  TRUE ~ NA_character_
            )
      )

sig_rows <- tab_df_PCA1$Sig %in% c("***", "**", "*", ".")

tab_PCA1 <- tab_df_PCA1 %>%
      gt() %>%
      cols_label(
            ci90_txt  = "CI 90%",
            ci95_txt  = "CI 95%",
            ci99_txt  = "CI 99%",
            ci999_txt = "CI 99.9%"
      ) %>%
      cols_hide(columns = max_ci) %>%
      fmt_number(columns = c(Estimate, SE), decimals = 3) %>%
      tab_header(title = "PCA1: Pale-Dark") %>%
      cols_align(align = "left", columns = Predictor) %>%
      cols_align(align = "center", columns = c(ci90_txt, ci95_txt, ci99_txt, ci999_txt, Sig)) %>%
      
      tab_style(
            style = cell_text(weight = "bold"),
            locations = cells_body(
                  columns = c(Predictor, Estimate),
                  rows = sig_rows
            )
      ) %>%
      
      tab_style(
            style = cell_text(size = "large", weight = "bold"),
            locations = cells_body(
                  columns = Sig,
                  rows = Sig %in% c("***", "**", ".")
            )
      ) %>%
      
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_body(columns = ci999_txt, rows = max_ci == "ci999_txt")) %>%
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_body(columns = ci99_txt,  rows = max_ci == "ci99_txt")) %>%
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_body(columns = ci95_txt,  rows = max_ci == "ci95_txt")) %>%
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_body(columns = ci90_txt,  rows = max_ci == "ci90_txt"))



N_points <- nrow(data_95_full)
R2_fixed_in_sample <- performance::r2(mod_sdm)$R2

tab_PCA1 <- tab_PCA1 %>%
      tab_source_note(
            source_note = paste0(
                  "R² (fixed effects) = ", sprintf("%.3f", R2_fixed_in_sample),
                  "; n = ", N_points, "."
            )
      ) %>%
      tab_source_note(
            source_note = "Moran's I = -0.0035 (p-value = 0.55)" #see below for the spatial autocor test
      )

# Table S1: Fixed‐effect estimates from the spatial SDM model explaining variation in PCA1 (Pale–Dark axis) using sdmTMB
tab_PCA1

set.seed(123)
samps <- sdmTMBextra::predict_mle_mcmc(mod_sdm, mcmc_iter = 1000, mcmc_warmup = 500)
mcmc_res <- residuals(mod_sdm, type = "mle-mcmc", mcmc_samples = samps)

#Q-Q plot

qq_df <- data.frame(
      sample = sort(mcmc_res),
      theo   = qnorm(ppoints(length(mcmc_res)))
)

p1_sdm <- ggplot(qq_df, aes(x = theo, y = sample)) +
      geom_point(color = "black", size = 1.5,alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color="darkred") +
      labs(
            title = "Q–Q plot of MCMC-based quantile residuals",
            x = "Theoretical Quantiles",
            y = "MCMC quantile residuals"
      ) +
      theme_minimal()


#Histogram
p2_sdm <- ggplot(data.frame(res = mcmc_res), aes(x = res)) +
      geom_histogram(bins = 30, fill = "black", color = "white", alpha = 0.5) +
      labs(
            title = "MCMC-based quantile residuals",
            x = "Residuals",
            y = "Count"
      ) +
      theme_minimal()


#Residuals vs fitted

# Fitted values 
fit_sdm <- predict(mod_sdm, re_form = NULL)$est  # conditionnal

p3_sdm <- ggplot(data.frame(fit = fit_sdm, res = mcmc_res),
                 aes(x = fit, y = res)) +
      geom_point(alpha = 0.5, color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
      labs(
            title = "Residuals vs fitted",
            x = "Fitted values",
            y = "MCMC quantile residuals"
      ) +
      theme_minimal()


# PCA1 Observed vs predicted
plot_df <- data.frame(
      observed = data_95_full$PCA1_mean_95rich,
      predicted = fit_sdm
)

p4_sdm <- ggplot(plot_df, aes(x = predicted, y = observed)) +
      geom_point(alpha = 0.5, color= "black") +
      geom_abline(intercept = 0, slope = 1, colour = "darkred", linetype = "dashed") +
      labs(
            title = "PCA1 (pale-dark) - Observed vs Predicted",
            x = "Predicted",
            y = "Observed"
      ) +
      theme_minimal()

# Figure S2: Diagnostic plots for the spatial model of the pale–dark colour axis (PCA1)
(p1_sdm | p2_sdm) / (p3_sdm | p4_sdm)


# Spatial autocor

coords_km <- cbind(data_95_full$X_km, data_95_full$Y_km)

# 3. Neighbours: points within 1446 km --> min too have at leat 1 neigboor per point
dmax_km <- 1446
nb_1446<- dnearneigh(coords_km, d1 = 0, d2 = dmax_km)

summary(nb_1446)  # check how many neighbours each point has
table(card(nb_1446))

listw_1446 <- nb2listw(nb_1446, style = "W")

# 5. Moran's I on sdmTMB residuals
set.seed(111)
moran_sdm_1446 <- spdep::moran.test(mcmc_res, listw_1446)
moran_sdm_1446 #### p-val=0.5487 (seed 111)

# Residuals map
# 
# ggplot(data_95_full,
#        aes(x = X_km, y = Y_km, colour = mcmc_res)) +
#       geom_point() +
#       scale_colour_gradient2(
#             midpoint = 0,
#             name = "Rsiduals"
#       ) +
#       coord_equal() +
#       theme_minimal() +
#       labs(
#             x = "Longitude",
#             y = "Latitude",
#             title = ""
#       )

### PCA2 Redness ### 
### PCA2 Redness ### 
### PCA2 Redness ### 

#  AIC selection for the best mesh

PCA2_cutoffs_km <- seq(100, 2000, by = 10)

PCA2_results <- list()

for (PCA2_c_km in PCA2_cutoffs_km) {
      cat("\n====================\nTesting cutoff =", PCA2_c_km, "km\n====================\n")
      
      # Mesh
      PCA2_mesh_c <- make_mesh(
            data_95_full,
            xy_cols = c("X_km", "Y_km"),
            cutoff = PCA2_c_km
      )
      
      # Fit model
      PCA2_mod_c <- try(
            sdmTMB(
                  PCA2_mean_95rich ~
                        count_species_log_scaled +
                        NMDS2_scaled +
                        SHAPE_AREA_log_scaled +
                        mean_bio10_scaled +
                        mean_bio16_scaled +
                        mean_bio17_scaled +
                        mean_UVB1_scaled +
                        mean_elev_scaled +
                        mean_canopy_scaled,
                  
                  data    = data_95_full,
                  mesh    = PCA2_mesh_c,
                  spatial = "on",
                  family  = student(link = "identity")
            ),
            silent = TRUE
      )
      
      if (inherits(PCA2_mod_c, "try-error")) {
            cat("!! cutoff", PCA2_c_km, "error, AIC = NA\n")
            PCA2_results[[as.character(PCA2_c_km)]] <- list(
                  mesh  = PCA2_mesh_c,
                  model = NA,
                  aic   = NA_real_
            )
      } else {
            PCA2_results[[as.character(PCA2_c_km)]] <- list(
                  mesh  = PCA2_mesh_c,
                  model = PCA2_mod_c,
                  aic   = AIC(PCA2_mod_c)
            )
      }
}

PCA2_res_cutoff <- data.frame(
      cutoff_km = PCA2_cutoffs_km,
      AIC       = sapply(PCA2_results, function(x) x$aic)
)

PCA2_res_cutoff_ok <- na.omit(PCA2_res_cutoff)

PCA2_res_cutoff_ok
plot(PCA2_res_cutoff_ok$cutoff_km, PCA2_res_cutoff_ok$AIC, type = "b")

PCA2_best_mesh <- PCA2_res_cutoff_ok[order(PCA2_res_cutoff_ok$AIC), ][1:10, ]
PCA2_best_mesh

PCA2_mesh <- make_mesh(
      data_95_full,
      type = "cutoff",
      xy_cols = c("X_km", "Y_km"),
      cutoff = 400 # <-- PCA2_best_mesh
)

mod_sdm2 <- sdmTMB(
      PCA2_mean_95rich ~
            
            count_species_log_scaled +
            NMDS2_scaled +
            SHAPE_AREA_log_scaled +
            mean_bio10_scaled +
            mean_bio16_scaled +
            mean_bio17_scaled +
            mean_UVB1_scaled +
            mean_elev_scaled +
            mean_canopy_scaled,
      
      data    = data_95_full,
      mesh    = PCA2_mesh,
      spatial = "on",
      family  = student(link = "identity")
)

summary(mod_sdm2)
sanity(mod_sdm2)
r2(mod_sdm2)

PCA2_coefs_90 <- tidy(mod_sdm2, effects = "fixed", conf.int = TRUE, conf.level = 0.90) %>%
      rename(conf.low.90 = conf.low, conf.high.90 = conf.high)

PCA2_coefs_95 <- tidy(mod_sdm2, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
      rename(conf.low.95 = conf.low, conf.high.95 = conf.high)

PCA2_coefs_99 <- tidy(mod_sdm2, effects = "fixed", conf.int = TRUE, conf.level = 0.99) %>%
      rename(conf.low.99 = conf.low, conf.high.99 = conf.high)

PCA2_coefs_999 <- tidy(mod_sdm2, effects = "fixed", conf.int = TRUE, conf.level = 0.999) %>%
      rename(conf.low.999 = conf.low, conf.high.999 = conf.high)

PCA2_coefs_all <- PCA2_coefs_95 %>%
      select(term, estimate, std.error) %>%
      left_join(PCA2_coefs_90  %>% select(term, conf.low.90,  conf.high.90),  by = "term") %>%
      left_join(PCA2_coefs_95  %>% select(term, conf.low.95,  conf.high.95),  by = "term") %>%
      left_join(PCA2_coefs_99  %>% select(term, conf.low.99,  conf.high.99),  by = "term") %>%
      left_join(PCA2_coefs_999 %>% select(term, conf.low.999, conf.high.999), by = "term") %>%
      mutate(
            significance = dplyr::case_when(
                  conf.low.999 > 0 | conf.high.999 < 0 ~ "***",
                  conf.low.99  > 0 | conf.high.99  < 0 ~ "**",
                  conf.low.95  > 0 | conf.high.95  < 0 ~ "*",
                  conf.low.90  > 0 | conf.high.90  < 0 ~ ".",
                  TRUE ~ ""
            )
      )

PCA2_coefs_all[, c("term", "estimate", "significance")]
PCA2_coefs_all

PCA2_var_labels <- c(
      "count_species_log_scaled"   = "Species richness (log)",
      "NMDS2_scaled"               = "Phylogenetic similarity",
      "SHAPE_AREA_log_scaled"      = "Ecoregion area (log)",
      "mean_bio10_scaled"          = "Temperature - warmest quarter (BIO10)",
      "mean_bio16_scaled"          = "Precipitation - wettest quarter (BIO16)",
      "mean_bio17_scaled"          = "Precipitation - driest quarter (BIO17)",
      "mean_UVB1_scaled"           = "UV-B radiation (UVB1)",
      "mean_elev_scaled"           = "Elevation",
      "mean_canopy_scaled"         = "Canopy height"
)

PCA2_tab_df <- PCA2_coefs_all %>%
      dplyr::transmute(
            Predictor = dplyr::recode(term, !!!PCA2_var_labels, "(Intercept)" = "(Intercept)"),
            Estimate  = estimate,
            SE        = std.error,
            ci90_txt  = sprintf("[%.3f, %.3f]", conf.low.90,  conf.high.90),
            ci95_txt  = sprintf("[%.3f, %.3f]", conf.low.95,  conf.high.95),
            ci99_txt  = sprintf("[%.3f, %.3f]", conf.low.99,  conf.high.99),
            ci999_txt = sprintf("[%.3f, %.3f]", conf.low.999, conf.high.999),
            Sig       = significance,
            
            max_ci    = dplyr::case_when(
                  Sig == "." ~ NA_character_,  
                  
                  (conf.low.999 > 0 & conf.high.999 > 0) | (conf.low.999 < 0 & conf.high.999 < 0) ~ "ci999_txt",
                  (conf.low.99  > 0 & conf.high.99  > 0) | (conf.low.99  < 0 & conf.high.99  < 0) ~ "ci99_txt",
                  (conf.low.95  > 0 & conf.high.95  > 0) | (conf.low.95  < 0 & conf.high.95  < 0) ~ "ci95_txt",
                  TRUE ~ NA_character_
            )
      )

PCA2_sig_rows <- PCA2_tab_df$Sig %in% c("***", "**", "*")

PCA2_tab <- PCA2_tab_df %>%
      gt::gt() %>%
      gt::cols_label(
            ci90_txt  = "CI 90%",
            ci95_txt  = "CI 95%",
            ci99_txt  = "CI 99%",
            ci999_txt = "CI 99.9%"
      ) %>%
      gt::cols_hide(columns = max_ci) %>%
      gt::fmt_number(columns = c(Estimate, SE), decimals = 3) %>%
      gt::tab_header(title = "PCA2: Redness") %>%
      gt::cols_align(align = "left",   columns = Predictor) %>%
      gt::cols_align(align = "center", columns = c(ci90_txt, ci95_txt, ci99_txt, ci999_txt, Sig)) %>%
      
      gt::tab_style(
            style = gt::cell_text(weight = "bold"),
            locations = gt::cells_body(
                  columns = c(Predictor, Estimate),
                  rows = PCA2_sig_rows
            )
      ) %>%
      
      gt::tab_style(
            style = gt::cell_text(size = "large", weight = "bold"),
            locations = gt::cells_body(
                  columns = Sig,
                  rows = Sig %in% c("***", "**", "*")
            )
      ) %>%
      
      gt::tab_style(
            style = gt::cell_text(weight = "bold"),
            locations = gt::cells_body(columns = ci999_txt, rows = max_ci == "ci999_txt")
      ) %>%
      gt::tab_style(
            style = gt::cell_text(weight = "bold"),
            locations = gt::cells_body(columns = ci99_txt, rows = max_ci == "ci99_txt")
      ) %>%
      gt::tab_style(
            style = gt::cell_text(weight = "bold"),
            locations = gt::cells_body(columns = ci95_txt, rows = max_ci == "ci95_txt")
      )

PCA2_N_points <- nrow(data_95_full)
PCA2_R2_fixed  <- performance::r2(mod_sdm2)$R2  

PCA2_tab <- PCA2_tab %>%
      tab_source_note(
            source_note = paste0(
                  "R² (fixed effects) = ", sprintf("%.3f", PCA2_R2_fixed),
                  "; n = ", PCA2_N_points, "."
            )
      ) %>%
      tab_source_note(
            source_note = "Moran's I = -0.0059 (p-value = 0.61)" #see below for the spatial autocor test
      )

# Table S2: Fixed‐effect estimates from the spatial SDM model explaining variation in PCA2 (redness axis) using sdmTMB

PCA2_tab

#MCMC-based quantile residuals

set.seed(123)
PCA2_samps <- sdmTMBextra::predict_mle_mcmc(mod_sdm2, mcmc_iter = 1000, mcmc_warmup = 500)
PCA2_mcmc_res <- residuals(mod_sdm2, type = "mle-mcmc", mcmc_samples = PCA2_samps)

PCA2_qq_df <- data.frame(
      sample = sort(PCA2_mcmc_res),
      theo   = qnorm(ppoints(length(PCA2_mcmc_res)))
)

PCA2_p1 <- ggplot(PCA2_qq_df, aes(x = theo, y = sample)) +
      geom_point(color = "black", size = 1.5, alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkred") +
      labs(
            title = "Q–Q plot of MCMC-based quantile residuals",
            x = "Theoretical Quantiles",
            y = "MCMC quantile residuals"
      ) +
      theme_minimal()

PCA2_p2 <- ggplot(data.frame(res = PCA2_mcmc_res), aes(x = res)) +
      geom_histogram(bins = 30, fill = "black", color = "white", alpha = 0.5) +
      labs(
            title = "MCMC-based quantile residuals",
            x = "Residuals",
            y = "Count"
      ) +
      theme_minimal()

PCA2_fit <- predict(mod_sdm2, re_form = NULL)$est

PCA2_p3 <- ggplot(data.frame(fit = PCA2_fit, res = PCA2_mcmc_res),
                  aes(x = fit, y = res)) +
      geom_point(alpha = 0.5, color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
      labs(
            title = "Residuals vs fitted",
            x = "Fitted values",
            y = "MCMC quantile residuals"
      ) +
      theme_minimal()


PCA2_plot_df <- data.frame(
      observed  = data_95_full$PCA2_mean_95rich,
      predicted = PCA2_fit
)

PCA2_p4 <- ggplot(PCA2_plot_df, aes(x = predicted, y = observed)) +
      geom_point(alpha = 0.5, color = "black") +
      geom_abline(intercept = 0, slope = 1, colour = "darkred", linetype = "dashed") +
      labs(
            title = "PCA2 (redness) - Observed vs Predicted",
            x = "Predicted",
            y = "Observed"
      ) +
      theme_minimal()

#Figure S3: Diagnostic plots for the spatial model of the redness axis (PCA2)

(PCA2_p1 | PCA2_p2) / (PCA2_p3 | PCA2_p4)

#Spatial autocor
coords_km <- cbind(data_95_full$X_km, data_95_full$Y_km)

# distance chosen so each point has >=1 neighbour
PCA2_dmax_km <- 1446
PCA2_nb <- spdep::dnearneigh(coords_km, d1 = 0, d2 = PCA2_dmax_km)

summary(PCA2_nb)
table(spdep::card(PCA2_nb))

PCA2_listw <- spdep::nb2listw(PCA2_nb, style = "W")

set.seed(111)
PCA2_moran <- spdep::moran.test(PCA2_mcmc_res, PCA2_listw)
PCA2_moran #pval=0.6117 (seed 111)

PCA2_map_df <- data_95_full
PCA2_map_df$PCA2_mcmc_res <- PCA2_mcmc_res

PCA2_res_map <- ggplot(PCA2_map_df,
                       aes(x = X_km, y = Y_km, colour = PCA2_mcmc_res)) +
      geom_point() +
      scale_colour_gradient2(midpoint = 0, name = "Residuals") +
      coord_equal() +
      theme_minimal() +
      labs(
            x = "Longitude",
            y = "Latitude",
            title = "PCA2 residuals (MCMC quantile residuals)"
      )

PCA2_res_map

############ bi-forest plot ############ 

# Forest plot from the model on darkness in X and redness in Y 
# Plot Figure 4 prior to visual enhancement in Inkscape

coef1 <- coefs_all %>%
      filter(term != "(Intercept)") %>%
      transmute(
            term,
            est1  = estimate,
            low1  = conf.low.95,
            high1 = conf.high.95,
            p1    = significance
      )

coef2 <- PCA2_coefs_all %>%
      filter(term != "(Intercept)") %>%
      transmute(
            term,
            est2  = estimate,
            low2  = conf.low.95,
            high2 = conf.high.95,
            p2    = significance
      )

coef_biv <- left_join(coef1, coef2, by = "term")

label_map <- c(
      "count_species_log_scaled"   = "Species richness (log)",
      "NMDS2_scaled"           = "Phylogenetic similarity",
      "SHAPE_AREA_log_scaled"  = "Ecoregion area (log)",
      "mean_bio10_scaled"      = "Temperature - warmest quarter (BIO10)",
      "mean_bio16_scaled"      = "Precipitation - wettest quarter (BIO16)",
      "mean_bio17_scaled"      = "Precipitation - driest quarter (BIO17)",
      "mean_UVB1_scaled"       = "UV-B radiation (UVB1)",
      "mean_elev_scaled"       = "Elevation",
      "mean_canopy_scaled"     = "Canopy height"
)

coef_biv <- coef_biv %>%
      mutate(
            label = label_map[term]
      )

coef_biv <- coef_biv %>%
      mutate(
            sig1 = p1 %in% c("*", "**", "***"),
            sig2 = p2 %in% c("*", "**", "***"),
            
            sig_cat = case_when(
                  sig1 & sig2 ~ "Both models",
                  sig1 & !sig2 ~ "Model 1 only",
                  !sig1 & sig2 ~ "Model 2 only",
                  TRUE         ~ "Non-significant"
            ),
            sig_cat = factor(sig_cat,
                             levels = c("Both models", "Model 1 only", "Model 2 only", "Non-significant"))
      )

#Bivariate forest plot

bi_forest_plot <- ggplot(coef_biv, aes(x = est1, y = est2)) +
      
      # IC horizontal (Model 1)
      geom_errorbarh(
            aes(xmin = low1, xmax = high1, alpha = sig_cat),
            height = 0.008,
            linewidth = 0.6
      ) +
      
      # IC vertical (Model 2)
      geom_errorbar(
            aes(ymin = low2, ymax = high2, alpha = sig_cat),
            width = 0.008,
            linewidth = 0.6
      ) +
      
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      
      geom_point(
            aes(colour = sig_cat, shape = sig_cat, alpha = sig_cat),
            size = 3
      ) +
      
      ggrepel::geom_text_repel(
            aes(label = label),
            size = 3.2,
            max.overlaps = Inf,
            box.padding = 0.4,
            point.padding = 0.3
      ) +
      
      scale_colour_manual(
            name = "Significance",
            values = c(
                  "Both models"     = "black",
                  "Model 1 only"    = "steelblue3",
                  "Model 2 only"    = "firebrick2",
                  "Non-significant" = "grey70"
            )
      ) +
      
      scale_shape_manual(
            name = "Significance",
            values = c(
                  "Both models"     = 16,
                  "Model 1 only"    = 17,
                  "Model 2 only"    = 15,
                  "Non-significant" = 1
            )
      ) +
      
      scale_alpha_manual(
            name = "Significance",
            values = c(
                  "Both models"     = 1,
                  "Model 1 only"    = 1,
                  "Model 2 only"    = 1,
                  "Non-significant" = 0.2
            )
      ) +
      
      scale_x_continuous(expand = expansion(mult = 0.05)) +
      scale_y_continuous(expand = expansion(mult = 0.05)) +
      
      labs(
            x = "Std estimate (Model 1: PCA1 Pale–Dark axis)",
            y = "Std estimate (Model 2: PCA2 Redness axis)"
      )


bi_forest_plot +
      theme_bw() +
      theme(
            panel.grid.major = element_line(
                  colour = "grey90",
                  linewidth = 0.4
            ),
            
            axis.text    = element_text(size = 12),
            axis.title   = element_text(size = 14),
            legend.title = element_text(size = 12),
            legend.text  = element_text(size = 12)
      )

