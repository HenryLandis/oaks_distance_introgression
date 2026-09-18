library(tidyverse)
library(dplyr)
library(ggplot2)
setwd("C:/Users/Henry/Desktop/WSU courses/Thesis/Analyses/RAD-seq/pixy/pixy_txs_100kb_hudson_hardfilter")

# Function to load pixy files.
pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
      
    }
    
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}

# Load pixy files.
pixy_files <- list.files(getwd(), full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

chromosomes <- c('chr1',
                 'chr2',
                 'chr3',
                 'chr4',
                 'chr5',
                 'chr6',
                 'chr7',
                 'chr8',
                 'chr9',
                 'chr10',
                 'chr11',
                 'chr12')


# The following code is for calculating pi/dxy/Fst
# in 100kb genomic windows. For pi, there are 12 data points
# per window, one for each TxS grouping. For dxy and Fst,
# there are 66 data points per window, one for each TxS pair.

# Average dxy of pairwise TxS comparisons per genomic window.
for (k in chromosomes){
  print(k)
  check1 <- filter(pixy_df, chromosome == k)
  check2 <- pull(check1, window_pos_1)
  wp1 = 1
  repeat{
    val1 <- filter(pixy_df, 
                 pop2 != 'NA'
                 & chromosome == k 
                 & window_pos_1 == wp1 
                 & statistic == 'count_diffs') %>% 
    pull(value) %>% 
    sum()
  if (wp1 %in% check2 == FALSE){
    break
  } else if (is.null(val1) == TRUE) {
    wp1 <- wp1 + 100000
  } else {
  val2 <- filter(pixy_df, 
                 pop2 != 'NA'
                 & chromosome == k 
                 & window_pos_1 == wp1 
                 & statistic == 'count_comparisons') %>% 
    pull(value) %>% 
    sum()
  pixy_df <- pixy_df %>% add_row(chromosome = k,
                                 window_pos_1 = wp1,
                                 window_pos_2 = wp1 + 99999,
                                 statistic = 'posthoc_dxy',
                                 value = val1 / val2)
  wp1 <- wp1 + 100000
    }
  }
}

# Average pi of Txs groups per genomic window.
for (k in chromosomes){
  print(k)
  check1 <- filter(pixy_df, chromosome == k)
  check2 <- pull(check1, window_pos_1)
  wp1 = 1
  repeat{
    val1 <- filter(pixy_df, 
                   is.na(pop2) == TRUE
                   & chromosome == k 
                   & window_pos_1 == wp1 
                   & statistic == 'count_diffs') %>% 
      pull(value) %>% 
      sum()
    if (wp1 %in% check2 == FALSE){
      break
    } else if (is.null(val1) == TRUE) {
      wp1 <- wp1 + 100000
    } else {
      val2 <- filter(pixy_df, 
                     is.na(pop2) == TRUE
                     & chromosome == k 
                     & window_pos_1 == wp1 
                     & statistic == 'count_comparisons') %>% 
        pull(value) %>% 
        sum()
      pixy_df <- pixy_df %>% add_row(chromosome = k,
                                     window_pos_1 = wp1,
                                     window_pos_2 = wp1 + 99999,
                                     statistic = 'posthoc_pi',
                                     value = val1 / val2)
      wp1 <- wp1 + 100000
    }
  }
}

# Average Fst of pairwise TxS comparisons per genomic window.
for (k in chromosomes){
  print(k)
  check1 <- filter(pixy_df, chromosome == k)
  check2 <- pull(check1, window_pos_1)
  wp1 = 1
  repeat{
    val1 <- filter(pixy_df, 
                   pop2 != 'NA'
                   & chromosome == k 
                   & window_pos_1 == wp1 
                   & statistic == 'wc_fst_a') %>% 
      pull(value) %>% 
      sum()
    if (wp1 %in% check2 == FALSE){
      break
    } else if (is.null(val1) == TRUE) {
      wp1 <- wp1 + 100000
    } else {
      val2 <- filter(pixy_df, 
                     pop2 != 'NA'
                     & chromosome == k 
                     & window_pos_1 == wp1 
                     & statistic == 'wc_fst_b') %>% 
        pull(value) %>% 
        sum()
      val3 <- filter(pixy_df, 
                     pop2 != 'NA'
                     & chromosome == k 
                     & window_pos_1 == wp1 
                     & statistic == 'wc_fst_c') %>% 
        pull(value) %>% 
        sum()
      pixy_df <- pixy_df %>% add_row(chromosome = k,
                                     window_pos_1 = wp1,
                                     window_pos_2 = wp1 + 99999,
                                     statistic = 'posthoc_fst',
                                     value = val1 / (
                                       val1 + val2 + val3))
      wp1 <- wp1 + 100000
    }
  }
}

# Using Hudson Fst instead.
# Average Fst of pairwise TxS comparisons per genomic window.
for (k in chromosomes){
  print(k)
  check1 <- filter(pixy_df, chromosome == k)
  check2 <- pull(check1, window_pos_1)
  wp1 = 1
  repeat{
    val1 <- filter(pixy_df, 
                   pop2 != 'NA'
                   & chromosome == k 
                   & window_pos_1 == wp1 
                   & statistic == 'hudson_fst_num') %>% 
      pull(value) %>% 
      sum()
    if (wp1 %in% check2 == FALSE){
      break
    } else if (is.null(val1) == TRUE) {
      wp1 <- wp1 + 100000
    } else {
      val2 <- filter(pixy_df, 
                     pop2 != 'NA'
                     & chromosome == k 
                     & window_pos_1 == wp1 
                     & statistic == 'hudson_fst_den') %>% 
        pull(value) %>% 
        sum()
      pixy_df <- pixy_df %>% add_row(
        chromosome = k,
        window_pos_1 = wp1,
        window_pos_2 = wp1 + 99999,
        statistic = 'posthoc_fst_hudson',
        value = val1 / val2)
      wp1 <- wp1 + 100000
    }
  }
}

# Create a custom labeller for special characters 
# in pi/dxy/Fst.
pixy_labeller <- as_labeller(c(
  posthoc_pi = "pi",
  posthoc_dxy = "D[XY]",
  # posthoc_fst = "F[ST]"),
  posthoc_fst_hudson = "F[ST]"),
  default = label_parsed)

# Manhattan plot of genomic window summary statistics 
# across all chromosomes.
pixy_df %>%
  mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
                                       chromosome == "X" ~ "even",
                                       TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = c(
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
    "chr7", "chr8", "chr9", "chr10", "chr11", "chr12"))) %>%
  filter(statistic %in% c("posthoc_pi", 
                          "posthoc_dxy", 
                          "posthoc_fst_hudson")) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
  geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
  facet_grid(statistic ~ chromosome,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Chromsome")+
  ylab("Statistic Value")+
  scale_color_manual(values = c("grey50", "black"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))


# The following code is for calculating genome-wide
# pi/dxy/Fst values. At base they consider the full
# dataset. With additional arguments, they can consider
# specific TxS or species groupings.

# Genome-wide dxy for entire dataset.
val1_dxy <- filter(pixy_df, 
                   pop2 != 'NA' & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_dxy <- filter(pixy_df, 
                   pop2 != 'NA'
                   & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
dxy_genome <- val1_dxy / val2_dxy

# Genome-wide dxy for three species comparisons.
# Do string matching for pop1 and pop2.
val1_dxy <- filter(pixy_df, 
                   grepl("GAM", pop1)
                   & grepl("GRI", pop2)
                   & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_dxy <- filter(pixy_df, 
                   grepl("GAM", pop1)
                   & grepl("GRI", pop2)
                   & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
dxy_genome <- val1_dxy / val2_dxy
# 0.007133676

val1_dxy <- filter(pixy_df, 
                   grepl("GAM", pop1)
                   & grepl("TUR", pop2)
                   & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_dxy <- filter(pixy_df, 
                   grepl("GAM", pop1)
                   & grepl("TUR", pop2)
                   & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
dxy_genome <- val1_dxy / val2_dxy
# 0.007200053

val1_dxy <- filter(pixy_df, 
                   grepl("GRI", pop1)
                   & grepl("TUR", pop2)
                   & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_dxy <- filter(pixy_df, 
                   grepl("GRI", pop1)
                   & grepl("TUR", pop2)
                   & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
dxy_genome <- val1_dxy / val2_dxy
# 0.004983223

# Genome-wide pi for entire dataset.
val1_pi <- filter(pixy_df, 
                   is.na(pop2) == TRUE 
                   & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004701714

# Genome-wide pi for 12 TxS groupings.
val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'PAY_GAM'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'PAY_GAM'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004281458

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'PAY_GRI'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'PAY_GRI'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004800728

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'PAY_TUR'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'PAY_TUR'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004726104

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'GIL_GAM'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'GIL_GAM'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004316693

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'GIL_GRI'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'GIL_GRI'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.00458493

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'LIN_GAM'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'LIN_GAM'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004605702

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'LIN_GRI'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'LIN_GRI'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004865751

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'GLO_GRI'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'GLO_GRI'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004803413

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'GLO_TUR'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'GLO_TUR'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004634786

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'SJN_GAM'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'SJN_GAM'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004451797

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'SFE_GRI'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'SFE_GRI'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.005338269

val1_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE
                  & pop1 == 'GVL_TUR'
                  & statistic == 'count_diffs') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_pi <- filter(pixy_df, 
                  is.na(pop2) == TRUE 
                  & pop1 == 'GVL_TUR'
                  & statistic == 'count_comparisons') %>%
  pull(value) %>% sum(na.rm = TRUE)
pi_genome <- val1_pi / val2_pi
# 0.004564722

# Barplot of genome-wide pi.
taxon <- c('Q. gambelii',
             'Q. grisea',
             'Q. turbinella',
             'Q. gambelii',
             'Q. grisea',
             'Q. gambelii',
             'Q. grisea',
             'Q. grisea',
             'Q. turbinella',
             'Q. gambelii',
             'Q. grisea',
             'Q. turbinella')
site <- c('Payson', 'Payson', 'Payson', 'Gila', 'Gila', 'Lincoln',
          'Lincoln', 'Globe', 'Globe', 'San Juan', 'Santa Fe',
          'Golden Valley')
values <- c(0.004281458,
            0.004800728,
            0.004726104,
            0.004316693,
            0.00458493,
            0.004605702,
            0.004865751,
            0.004803413,
            0.004634786,
            0.004451797,
            0.005338269,
            0.004564722
)
data <- data.frame(taxon, site, values) %>%
  mutate(site = factor(
    site, levels = c('Payson', 'Globe', 'Gila', 'Lincoln',
                     'San Juan', 'Santa Fe', 'Golden Valley')
  ))
ggplot(data, aes(x = site, y = values, fill = taxon)) +
  geom_col(position = position_dodge2(preserve = 'single'), 
           stat = "identity") +
  scale_fill_manual(values = c("#d95f02", "#7570b3", "#1b9e77")) +
  labs(x = NULL, y = "pi") +
  lims(y = c(0, 0.006)) +
  theme_minimal()

# geom_col()# Genome-wide WC Fst for entire dataset.
val1_fst <- filter(pixy_df, 
                   pop2 != 'NA' & statistic == 'wc_fst_a') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_fst <- filter(pixy_df, 
                   pop2 != 'NA' & statistic == 'wc_fst_b') %>%
  pull(value) %>% sum(na.rm = TRUE)
val3_fst <- filter(pixy_df, 
                   pop2 != 'NA' & statistic == 'wc_fst_c') %>%
  pull(value) %>% sum(na.rm = TRUE)
fst_genome <- val1_fst / (val1_fst + val2_fst + val3_fst)

# Genome-wide Hudson Fst for entire dataset.
val1_fsth <- filter(pixy_df, 
                    pop2 != 'NA'
                  & statistic == 'hudson_fst_num') %>%
  pull(value) %>% sum(na.rm = TRUE)
val2_fsth <- filter(pixy_df, 
                    pop2 != 'NA'
                  & statistic == 'hudson_fst_den') %>%
  pull(value) %>% sum(na.rm = TRUE)
fsth_genome <- val1_fsth / val2_fsth
