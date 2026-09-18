library(tidyverse)
library(ggplot2)
setwd("C:/Users/Henry/Desktop/WSU courses/Thesis/Lab protocols/RAD-seq/PCA")

pca <- read.table("pca_075.eigenvec", header = FALSE)
eigenval <- scan("pca_075.eigenval")

# Remove nuisance column.
pca <- pca[,-1]

# Alter first column.
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# pca$ind <- substr(pca$ind, 1, 6)

# Load species and location info.
info <- read.csv("spp_loc_new.csv")
names(info)[1] <- "ind"

# Merge data frames.
combine <- merge(pca, info, by = "ind", all = TRUE)

# Convert eigenvalues to percentage variance explained.
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

# Plot variance explained by each principal component.
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()


# Subset to keep only allopatric points.
allopatry <- combine %>%
  filter(Location == 'San Juan' |
           Location == 'Santa Fe' |
           Location == 'Golden Valley')

# PCA with only allopatry (for slides).
ggplot(allopatry, aes(PC1, PC2, col = Species, shape = Location)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("#d95f02", 
                                 "#7570b3", 
                                 "#1b9e77")) +
  scale_shape_manual(breaks = c('San Juan',
                                'Santa Fe',
                                'Golden Valley'),
    values = c(16, 17, 15)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  coord_equal() +
  theme_light()


# PCA with all samples.
ggplot(combine, aes(PC1, PC2, col = Species, shape = Location)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("#d95f02", 
                                 "#7570b3", 
                                 "#1b9e77")) +
  scale_shape_manual(breaks = c('Payson',
                                'Globe',
                                'Gila',
                                'Lincoln',
                                'San Juan',
                                'Santa Fe',
                                'Golden Valley'),
    values = c(5, 1, 2, 0, 16, 17, 15)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  coord_equal() +
  theme_light()



#b <- ggplot(combine, aes(
#  PC1, PC2, col = Species, shape = Location)) + 
#  geom_point(size = 3)
#b <- b + scale_colour_manual(values = 
#                               c("#d95f02", 
#                                 "#7570b3", 
#                                 "#1b9e77"))
#b <- b + scale_shape_manual(
#  values = c(0, 1, 15, 2, 5, 16, 17))
#b <- b + coord_equal() + theme_light()
#b + xlab(
#  paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
#  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))