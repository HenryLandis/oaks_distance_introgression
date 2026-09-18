library(geosphere)
library(gtools)
library(vegan)
library(tidyverse)
library(ggplot2)


setwd("C:/Users/Henry/Desktop/WSU courses/Thesis/Analyses/RAD-seq/hierfstat")

# Get centroids of 12 taxon-site combos.
pol_PAY_GAM <- read.csv("./centroid/pay_gam.csv", header = FALSE)
centroid_PAY_GAM <- centroid(pol_PAY_GAM)
pol_PAY_GRI <- read.csv("./centroid/pay_gri.csv", header = FALSE)
centroid_PAY_GRI <- centroid(pol_PAY_GRI)
pol_PAY_TUR <- read.csv("./centroid/pay_tur.csv", header = FALSE)
centroid_PAY_TUR <- centroid(pol_PAY_TUR)
pol_GIL_GAM <- read.csv("./centroid/gil_gam.csv", header = FALSE)
centroid_GIL_GAM <- centroid(pol_GIL_GAM)
pol_GIL_GRI <- read.csv("./centroid/gil_gri.csv", header = FALSE)
centroid_GIL_GRI <- centroid(pol_GIL_GRI)
pol_LIN_GAM <- read.csv("./centroid/lin_gam.csv", header = FALSE)
centroid_LIN_GAM <- centroid(pol_LIN_GAM)
pol_LIN_GRI <- read.csv("./centroid/lin_gri.csv", header = FALSE)
centroid_LIN_GRI <- centroid(pol_LIN_GRI)
pol_GLO_GRI <- read.csv("./centroid/glo_gri.csv", header = FALSE)
centroid_GLO_GRI <- centroid(pol_GLO_GRI)
pol_GLO_TUR <- read.csv("./centroid/glo_tur.csv", header = FALSE)
centroid_GLO_TUR <- centroid(pol_GLO_TUR)
pol_SJN_GAM <- read.csv("./centroid/sjn.csv", header = FALSE)
centroid_SJN_GAM <- centroid(pol_SJN_GAM)
pol_SFE_GRI <- read.csv("./centroid/sfe.csv", header = FALSE)
centroid_SFE_GRI <- centroid(pol_SFE_GRI)
pol_GLV_TUR <- read.csv("./centroid/glv.csv", header = FALSE)
centroid_GLV_TUR <- centroid(pol_GLV_TUR)

# Get centroids for all sampled individuals in sympatric sites.
# Used for sampling sites figure.
pol_Payson<- rbind(pol_PAY_GAM, pol_PAY_GRI, pol_PAY_TUR)
centroid_Payson <- centroid(pol_Payson)
pol_Globe <- rbind(pol_GLO_GRI, pol_GLO_TUR)
centroid_Globe <- centroid(pol_Globe)
pol_Gila <- rbind(pol_GIL_GAM, pol_GIL_GRI)
centroid_Gila <- centroid(pol_Gila)
pol_Lincoln <- rbind(pol_LIN_GAM, pol_LIN_GRI)
centroid_Lincoln <- centroid(pol_Lincoln)

centroid_list <- c(centroid_PAY_GAM, centroid_PAY_GRI,
                   centroid_PAY_TUR, centroid_GIL_GAM,
                   centroid_GIL_GRI, centroid_LIN_GAM,
                   centroid_LIN_GRI, centroid_GLO_GRI,
                   centroid_GLO_TUR, centroid_SJN_GAM,
                   centroid_SFE_GRI, centroid_GLV_TUR)

# Use distHaversine() function for distances between centroids.
# Values recorded in CSV files.

tbl <- read_csv('pairwise_fst_km_txs.csv')
tbl <- tbl %>% add_column(ln_km = log(tbl$hav_dist_km))

# All samples.
label = "r = 0.03318
p = 0.0592"
ggplot(tbl, aes(ln_km, Fst, fill = factor(code),
       shape = factor(pairtype))) +
  geom_point(size = 2.5, 
             stroke = 1.25,
             position = 'jitter') +
  #geom_text(aes(label = txs1), 
  #          nudge_x = 0.6,
  #          nudge_y = 0.005) +
  xlab('ln (km)') +
  ylab('Fst / (1 - Fst)') +
  labs(fill = 'Taxon pairing', shape = 'Site type pairing') +
  geom_text(x = -0.5, y = 0.1, label = label) +
  scale_fill_manual(values = c("#bf812d", 
                               "#bababa", 
                               "#35978f")) + 
  scale_shape_manual(values = c(21, 22, 24)) +
  guides(fill = guide_legend(
    override.aes = list(shape=23))) +
  theme_bw()
# theme(legend.position = "None")
# hjust/vjust can be used to manually adjust labels if I'm crazy

# Only crosses with different sites and same taxa,
# meaning intraspecific crosses.
tbl_DS <- tbl[tbl$code == 'Different site, same taxa',]
label_DS = "r = 0.3176
p = 0.0114"
ggplot(tbl_DS, aes(ln_km, Fst)) +
  geom_point(size = 2.5, 
             stroke = 1.25,
             position = 'jitter') +
  geom_text(x = 4.8, y = 0.023, label = label_DS) +
  #geom_text(aes(label = txs1), 
  #          nudge_x = 0.6,
  #          nudge_y = 0.005) +
  xlab('ln (km)') +
  ylab('Fst / (1 - Fst)') +
  scale_shape_manual(values = c(22, 24)) +
  theme_bw()

# Only crosses with different taxa, meaning interspecific
# crosses.
tbl_inter <- subset(tbl, 
                    code == "Same site, different taxa" | 
                    code == "Different site, different taxa")
label_inter = "r = 0.03573
p = 0.2456"
ggplot(tbl_inter, aes(ln_km, Fst)) +
  geom_point(size = 2.5, 
             stroke = 1.25,
             position = 'jitter') +
  geom_text(x = 0, y = 0.12, label = label_inter) +
  #geom_text(aes(label = txs1), 
  #          nudge_x = 0.6,
  #          nudge_y = 0.005) +
  xlab('ln (km)') +
  ylab('Fst / (1 - Fst)') +
  scale_shape_manual(values = c(22, 24)) +
  theme_bw()


# Mantel test to show correlation between geodist & Fst
dist_geo = dist(tbl$ln_km)
dist_fst = dist(tbl$Fst)
mtest = mantel(x = dist_geo, 
               y = dist_fst,
               method = 'pearson',
               permutations = 9999)
mtest
# r = 0.03318, p = 0.0592

tbl_AA <- tbl[tbl$pairtype == 'Allopatric-allopatric',]
distAA_geo = dist(tbl_AA$ln_km)
distAA_fst = dist(tbl_AA$Fst)
mtestAA = mantel(x = distAA_geo, 
               y = distAA_fst,
               method = 'pearson',
               permutations = 9999)
mtestAA
# r = -0.69000, p = 0.8333

tbl_AS <- tbl[tbl$pairtype == 'Allopatric-sympatric',]
distAS_geo = dist(tbl_AS$ln_km)
distAS_fst = dist(tbl_AS$Fst)
mtestAS = mantel(x = distAS_geo, 
                 y = distAS_fst,
                 method = 'pearson',
                 permutations = 9999)
mtestAS
# r = -0.01972, p = 0.5721

tbl_SS <- tbl[tbl$pairtype == 'Sympatric-sympatric',]
distSS_geo = dist(tbl_SS$ln_km)
distSS_fst = dist(tbl_SS$Fst)
mtestSS = mantel(x = distSS_geo, 
                 y = distSS_fst,
                 method = 'pearson',
                 permutations = 9999)
mtestSS
# r = 0.01113, p = 0.2673

tbl_DD <- tbl[tbl$code == 'Different site, different taxa',]
distDD_geo = dist(tbl_DD$ln_km)
distDD_fst = dist(tbl_DD$Fst)
mtestDD = mantel(x = distDD_geo, 
                 y = distDD_fst,
                 method = 'pearson',
                 permutations = 9999)
mtestDD
# r = 0.08349, p = 0.0775

tbl_DS <- tbl[tbl$code == 'Different site, same taxa',]
distDS_geo = dist(tbl_DS$ln_km)
distDS_fst = dist(tbl_DS$Fst)
mtestDS = mantel(x = distDS_geo, 
                 y = distDS_fst,
                 method = 'pearson',
                 permutations = 9999)
mtestDS
# r = 0.3176, p = 0.0114 *****

tbl_SD <- tbl[tbl$code == 'Same site, different taxa',]
distSD_geo = dist(tbl_SD$ln_km)
distSD_fst = dist(tbl_SD$Fst)
mtestSD = mantel(x = distSD_geo, 
                 y = distSD_fst,
                 method = 'pearson',
                 permutations = 9999)
mtestSD
# r = -0.2072, p = 0.72917

# All interspecific crosses.
tbl_inter <- subset(tbl, 
                    code == "Same site, different taxa" | 
                    code == "Different site, different taxa")
dist_inter_geo = dist(tbl_inter$ln_km)
dist_inter_fst = dist(tbl_inter$Fst)
mtestSD = mantel(x = dist_inter_geo, 
                 y = dist_inter_fst,
                 method = 'pearson',
                 permutations = 9999)
mtestSD
# r = 0.03573, p = 0.2456






# Mantel test variation permuting only rows, not columns,
# on interspecific crosses only.
# Original Mantel result: r = 0.08349, p = 0.0775
distDD_geo = dist(tbl_DD$ln_km)
distDD_fst = dist(tbl_DD$Fst)

perms <- numeric(9999)
for (x in 1:9999) {
  temp <- as.matrix(distDD_geo) # is this the original matrix?
  rmat <- temp[sample(nrow(temp)), ]
  temp2 <- as.dist(rmat)
  value <- cor(temp2, distDD_fst)
  perms[x] <- value
}
pvalue <- sum(perms > 0.08349) / length(perms)
pvalue # 0.159616


