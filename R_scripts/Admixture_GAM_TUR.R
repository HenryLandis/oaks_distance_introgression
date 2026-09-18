library(tidyverse)
setwd("C:/Users/Henry/Desktop/WSU courses/Thesis/Analyses/RAD-seq/Admixture")

# K/CV values.
K <- c(1:12)
CV <- c(0.600716,
        0.431540,
        0.454068,
        0.476568,
        0.504929,
        0.536269,
        0.577581,
        0.618423,
        0.664226,
        0.708610,
        0.757400,
        0.807461)
admix <- data.frame(K, CV)

# Plot CV chart.
ggplot(data = admix, aes(x = K, y = CV)) +
  geom_line() +
  labs(title = "Cross-validation error: GAM & TUR") +
  scale_x_continuous(breaks = seq(1, 12, by = 1))+
  theme_minimal()

# Load table.
tbl=read.table("pca_075_gam_tur.2_2.Q")

# Load sample set.
samples=read.csv('spp_loc_GAM_TUR.csv')
colnames(samples)=c('ID','taxon','site')

# Bind genetic data and labels.
tbl=cbind(tbl, samples)
colnames(tbl)[1:2]=c('Q1','Q2')

# Order by site and/or taxon.
# tbl=tbl[order(tbl$site,tbl$taxon), ]
# ranks=1:nrow(tbl)
# tbl$order=ranks

# Shorthand taxon labels.
tbl[tbl$taxon=='Quercus gambelii',]$taxon='GA'
tbl[tbl$taxon=='Quercus turbinella',]$taxon='TU'
# tbl$label=paste(tbl$site,tbl$order,sep=':')

# Subset by site or taxon.
# tbl=subset(tbl,tbl$site %in% c('Gila','Globe','Payson','Lincoln'))

# Make room for vertical label names.
# par(mar=c(4,4,4,0.01))

# Get the names as a standalone vector.
# namelist <- c(tbl[,3])

# Get "long" data.
tbl2 <- tbl %>% pivot_longer(cols = starts_with("Q"), 
                             names_to = "cluster", 
                             values_to = "q")

# Get cluster assignments by Q value.
tbl3 <- tbl2 %>% group_by(ID) %>% mutate(
  likely_assignment = which.min(q), # can only arrange by min for k=2
  assignment_prob = min(q)) %>%
  arrange(case_when(likely_assignment == 1 ~ assignment_prob,
                    likely_assignment == 2 ~ -assignment_prob)) %>%
  arrange(match(likely_assignment, c(1, 2))) %>%
  ungroup() %>% 
  mutate(ID = forcats::fct_inorder(factor(ID)))

# Arrange sympatric sites first, then allopatric.
tbl4 <- tbl3 %>% arrange(match(site, c('Payson',
                                       'Globe',
                                       'Gila',
                                       'Lincoln',
                                       'San Juan',
                                       'Golden Valley'
)))

# Build plot.
plot1 <- tbl4 %>% ggplot() +
  # Refactoring so the gambelii cluster is "first" for reading plot.
  geom_col(aes(x = ID, y = q, fill = factor(
    cluster, levels = c('Q2', 'Q1'))),
    linewidth = 0.05, color = 'black') +
  scale_fill_manual(values = c('#fddaec', '#ffff33'),
                    labels = c('Q1', 'Q2')) +
  facet_grid(~fct_inorder(site), 
             switch = "x", scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_x_discrete(expand = expansion(add = 1)) +
  labs(x = NULL, y = "Q-values", fill = "Cluster") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

tbl5 <- tbl4 %>% filter(row_number() %% 2 != 0) %>%
  mutate(dummy = 1)
plot2 <- tbl5 %>% ggplot() +
  geom_col(aes(x = ID, y = dummy, fill = taxon),
           linewidth = 0.05, color = 'black') +
  scale_fill_manual(values = c("#d95f02", 
                               "#1b9e77")) +
  facet_grid(~fct_inorder(site), 
             switch = "x", scales = "free", space = "free") +
  labs(x = NULL, y = NULL, fill = "Taxon") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank())

plot1 / plot2 + plot_layout(heights = (c(60, 1)),
                            guides = 'collect')
