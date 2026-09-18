library(tidyverse)
library(ggplot2)
setwd("C:/Users/Henry/Desktop/WSU courses/Thesis/Lab protocols/RAD-seq/PCA")

pca <- read_table("pca_075_gam_tur.eigenvec")
eigenval <- scan("pca_075_gam_tur.eigenval")

# Remove nuisance column.
pca <- pca[,-1]

# Alter first column.
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# pca$ind <- substr(pca$ind, 1, 6)

# Load species and location info.
info <- read.csv("spp_loc_GAM_TUR.csv")
names(info)[1] <- "ind"

# Merge data frames.
combine <- merge(pca, info, by = "ind", all = TRUE)

# Convert eigenvalues to percentage variance explained/
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

# Plot variance explained by each principal component.
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# Plot PCA.
ggplot(combine, aes(PC1, PC2, col = Species, shape = Location)) + 
  geom_point(size = 3) +
  scale_colour_manual(values = c("#d95f02", "#1b9e77")) +
  scale_shape_manual(breaks = c('Payson',
                                'Globe',
                                'Gila',
                                'Lincoln',
                                'San Juan',
                                'Golden Valley'),
                     values = c(5, 1, 2, 0, 16, 15)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  coord_equal() +
  theme_light()
