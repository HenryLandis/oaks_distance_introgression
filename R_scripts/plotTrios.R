library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(corrplot)

setwd('C:/Users/Henry/Desktop/WSU courses/Thesis/Analyses/RAD-seq/Dsuite/')

PAY_trios <- read.csv('inds_payson_tree_FDR.csv')
GLO_trios <- read.csv('inds_globe_tree_FDR.csv')
GIL_trios <- read.csv('inds_gila_tree_FDR.csv')
LIN_trios <- read.csv('inds_lincoln_tree_FDR.csv')
SJN_trios <- read.csv('inds_sanjuan_FDR.csv')
SFE_trios <- read.csv('inds_santafe_FDR.csv')
GLV_trios <- read.csv('inds_goldenvalley_FDR.csv')
# The dashes in 'p-value' need to be removed.


# Payson
p01 <- ggplot() +
  geom_histogram(data = PAY_trios, aes(x = Dstatistic),
                 binwidth = 0.001, linewidth = 0.05,
                 fill = 'lightgray', color = 'black') +
  geom_vline(data = PAY_trios, aes(
    xintercept = mean(Dstatistic)),
    color = 'tan') +
  labs(title = "Payson (n=6545)") +
  xlim(0, 0.3) +
  theme_bw()
p02 <- ggplot() +
  geom_vline(aes(xintercept = 0.05)) +
  geom_histogram(data = PAY_trios, aes(x = pvalue),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'red', color = 'black',
                 alpha = 0.25) +
  geom_histogram(data = PAY_trios, aes(x = pvalue_FDR),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'blue', color = 'black',
                 alpha = 0.25) +
  xlim(0, 1) + ylim(0, 180) +
  theme_bw()

# Globe
p03 <- ggplot() +
  geom_histogram(data = GLO_trios, aes(x = Dstatistic),
                 binwidth = 0.001, linewidth = 0.05,
                 fill = 'lightgray', color = 'black') +
  geom_vline(aes(
    xintercept = mean(GLO_trios$Dstatistic)),
    color = 'tan') +
  labs(title = "Globe (n=1330)") +
  xlim(0, 0.3) + ylim(0, 80) + 
  theme_bw()
p04 <- ggplot() +
  geom_vline(aes(xintercept = 0.05)) +
  geom_histogram(data = GLO_trios, aes(x = pvalue),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'red', color = 'black',
                 alpha = 0.25) +
  geom_histogram(data = GLO_trios, aes(x = pvalue_FDR),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'blue', color = 'black',
                 alpha = 0.25) +
  xlim(0, 1) + ylim(0, 180) +
  theme_bw()

# Gila
p05 <- ggplot() +
  geom_histogram(data = GIL_trios, aes(x = Dstatistic),
                 binwidth = 0.001, linewidth = 0.05,
                 fill = 'lightgray', color = 'black') +
  geom_vline(aes(
    xintercept = mean(GIL_trios$Dstatistic)),
    color = 'tan') +
  labs(title = "Gila (n=1330)") +
  xlim(0, 0.3) + ylim(0, 80) + 
  theme_bw()
p06 <- ggplot() +
  geom_vline(aes(xintercept = 0.05)) +
  geom_histogram(data = GIL_trios, aes(x = pvalue),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'red', color = 'black',
                 alpha = 0.25) +
  geom_histogram(data = GIL_trios, aes(x = pvalue_FDR),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'blue', color = 'black',
                 alpha = 0.25) +
  xlim(0, 1) + ylim(0, 180) +
  theme_bw()

# Lincoln
p07 <- ggplot() +
  geom_histogram(data = LIN_trios, aes(x = Dstatistic),
                 binwidth = 0.001, linewidth = 0.05,
                 fill = 'lightgray', color = 'black') +
  geom_vline(aes(
    xintercept = mean(LIN_trios$Dstatistic)),
    color = 'tan') +
  labs(title = "Lincoln (n=1140)") +
  xlim(0, 0.3) + ylim(0, 80) + 
  theme_bw()
p08 <- ggplot() +
  geom_vline(aes(xintercept = 0.05)) +
  geom_histogram(data = LIN_trios, aes(x = pvalue),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'red', color = 'black',
                 alpha = 0.25) +
  geom_histogram(data = LIN_trios, aes(x = pvalue_FDR),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'blue', color = 'black',
                 alpha = 0.25) +
  xlim(0, 1) + ylim(0, 180) +
  theme_bw()

# San Juan
p09 <- ggplot() +
  geom_histogram(data = SJN_trios, aes(x = Dstatistic),
                 binwidth = 0.001, linewidth = 0.05,
                 fill = 'lightgray', color = 'black') +
  geom_vline(aes(
    xintercept = mean(SJN_trios$Dstatistic)),
    color = 'tan') +
  labs(title = "San Juan (n=57)") +
  xlim(0, 0.3) + ylim(0, 8) +
  theme_bw()
p10 <- ggplot() +
  geom_vline(aes(xintercept = 0.05)) +
  geom_histogram(data = SJN_trios, aes(x = pvalue),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'red', color = 'black',
                 alpha = 0.25) +
  geom_histogram(data = SJN_trios, aes(x = pvalue_FDR),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'blue', color = 'black',
                 alpha = 0.25) +
  xlim(0, 1) + ylim(0, 40) +
  theme_bw()

# Santa Fe
p11 <- ggplot() +
  geom_histogram(data = SFE_trios, aes(x = Dstatistic),
                 binwidth = 0.001, linewidth = 0.05,
                 fill = 'lightgray', color = 'black') +
  geom_vline(aes(
    xintercept = mean(SFE_trios$Dstatistic)),
    color = 'tan') +
  labs(title = "Santa Fe (n=220)") +
  xlim(0, 0.3) + ylim(0, 8) +
  theme_bw()
p12 <- ggplot() +
  geom_vline(aes(xintercept = 0.05)) +
  geom_histogram(data = SFE_trios, aes(x = pvalue),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'red', color = 'black',
                 alpha = 0.25) +
  geom_histogram(data = SFE_trios, aes(x = pvalue_FDR),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'blue', color = 'black',
                 alpha = 0.25) +
  xlim(0, 1) + ylim(0, 40) +
  theme_bw()

# Golden Valley
p13 <- ggplot() +
  geom_histogram(data = GLV_trios, aes(x = Dstatistic),
                 binwidth = 0.001, linewidth = 0.05,
                 fill = 'lightgray', color = 'black') +
  geom_vline(aes(
    xintercept = mean(PAY_trios$Dstatistic)),
    color = 'tan') +
  labs(title = "Golden Valley (n=120)") +
  xlim(0, 0.3) + ylim(0, 8) +
  theme_bw()
p14 <- ggplot() +
  geom_vline(aes(xintercept = 0.05)) +
  geom_histogram(data = GLV_trios, aes(x = pvalue),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'red', color = 'black',
                 alpha = 0.25) +
  geom_histogram(data = GLV_trios, aes(x = pvalue_FDR),
                 binwidth = 0.01, linewidth = 0.05,
                 fill = 'blue', color = 'black',
                 alpha = 0.25) +
  xlim(0, 1) + ylim(0, 40) +
  theme_bw()

# Patchwork multi panel arrangement.
(p01 + p02) / (p03 + p04) / (p05 + p06) / (p07 + p08)
(p09 + p10) / (p11 + p12) / (p13 + p14)


# KS tests for all site pairs.
# ks.test(PAY_trios$Dstatistic, GIL_trios$Dstatistic)
# etc

# Plot KS results.
ksmat <- as.matrix(read.csv('ks_tests.csv', row.names = 1))
ksmat_p <- as.matrix(read.csv('ks_tests_p.csv', row.names = 1))

trace(corrplot, edit=TRUE)
# Line 449 must be edited to the following:
# text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs]-0.25, 
# Must be repeated on every startup of R.

corrplot(ksmat, method = 'number', 
         type = 'lower', 
         is.corr = FALSE,
         col = COL2('PiYG'),
         number.digits = 3,
         number.cex = 0.95,
         tl.cex = 0.8,
         tl.col = 'black',
         tl.srt = 50,
         p.mat = ksmat_p,
         sig.level = c(0.001, 0.05),
         insig = 'label_sig',
         pch.cex = 0.9)
