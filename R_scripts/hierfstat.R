library(gaston)
library(hierfstat)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(vcfR)
library(corrplot)

setwd("C:/Users/Henry/Desktop/WSU courses/Thesis/Analyses/RAD-seq/hierfstat")
x1 <- read.VCF("genotypes_final.vcf", convert.chr = FALSE,
               BiAllelic = TRUE)
dosage <- as.matrix(x1) # to matrix
df <- as.data.frame((dosage)) # to dataframe
names(df) <- 1:ncol(df)

# Get location info. 
info <- read.csv("spp_loc.csv")


# Extract vector, and insert into dataframe.
vec <- info$BIN
df <- add_column(df, new_col = vec, .before = 1,
                 .name_repair = 'minimal')

# Threshold to filter out loci with excess missing data.
df2 <- select(df, where(~ is.numeric(.x) && (mean(!is.na(.x))) >= 0.5))

# the long way
# df[136,] = NA
# for(i in 2:ncol(df)){
#   df[136,i] = sum(is.na(df[1:135,i])) / 135
#  }
# hist(df[136,])

bs <- basic.stats(df2)
summary(bs)
bs$overall

# Use perloc to observe distribution of heterozygosity.
fis <- bs$perloc$Fis
hist(fis, main="FIS (population assignment)")


# Variant with species assignment instead of population.
vec2 <- info$SPE
df3 <- add_column(df, new_col = vec2, .before = 1)
df4 <- select(df3, where(~ is.numeric(.x) && (mean(!is.na(.x))) >= 0.5))

bs_species <- basic.stats(df4)
summary(bs_species)
bs_species$overall

fis2 <- bs_species$perloc$Fis
hist(fis2, main="FIS (species assignment)")


# Variant with species-site groupings.
vec3 <- info$TxS
df5 <- add_column(df, new_col = vec3, .before = 1)
df6 <- select(df5, where(~ is.numeric(.x) && (mean(!is.na(.x))) >= 0.5))

bs_txs <- basic.stats(df6)
summary(bs_txs)
bs_txs$overall

fis3 <- bs_txs$perloc$Fis
hist(fis3, main="FIS (TxS assignment)")

# ggplot version of the histogram.
perloc <- bs_txs$perloc
ggplot() +
  geom_histogram(data = perloc, aes(x = Fis), binwidth = 0.01,
                 fill = 'lightgray', color = 'black') +
  xlab(bquote(F[IS])) +
  theme_bw()


# FST calculations.
pairwise_matrix <- pairwise.WCfst(df2)
popnames <- c("Payson_S", "Gila_S", "Lincoln_S", "Globe_S",
          "SantaFe_A", "SanJuan_A", "GValley_A")
rownames(pairwise_matrix) <- popnames
colnames(pairwise_matrix) <- popnames
write.csv(pairwise_matrix, "FST_pairwise_populations.csv")

pairwise_matrix_species <- pairwise.WCfst(df4)
speciesnames <- c("Quercus gambelli", "Quercus grisea",
                  "Quercus turbinella")
rownames(pairwise_matrix_species) <- speciesnames
colnames(pairwise_matrix_species) <- speciesnames
write.csv(pairwise_matrix_species, "FST_pairwise_species.csv")

pairwise_matrix_txs <- pairwise.WCfst(df6)
txs_names <- c('PAY_GAM', 'PAY_GRI', 'PAY_TUR',
               'GLO_GRI', 'GLO_TUR', 'GIL_GAM', 
               'GIL_GRI', 'LIN_GAM', 'LIN_GRI', 
               'SJN_GAM', 'SFE_GRI', 'GLV_TUR'
               )
rownames(pairwise_matrix_txs) <- txs_names
colnames(pairwise_matrix_txs) <- txs_names
write.csv(pairwise_matrix_txs, "FST_pairwise_txs.csv")


# Plot Fst as a heatmap.
df <- read.csv('FST_pairwise_txs_half.csv', row.names = 1)
df2 <- as.matrix(df, diag = FALSE)
corrplot(df2, method = 'number', 
         type = 'lower', 
         is.corr = FALSE,
         col = COL2('PiYG'),
         number.digits = 3,
         number.cex = 0.95,
         tl.cex = 0.8,
         tl.srt = 50,
         tl.col = 'black')

# Fst violin plot.
df3 <- read.csv('violin_fst.csv')
ggplot(df3, aes(x = factor(pair), y = Fst)) +
  geom_violin(fill = 'lightblue', trim = FALSE) +
  geom_jitter(width = 0.08) +
  labs(x = "Taxon pair", y = 'Fst') +
  theme_bw()
