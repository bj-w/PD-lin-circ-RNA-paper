library(tximport)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(data.table)
library(tidyverse)
library(janitor)
library(ragg)

source('combinedLinCirc/PD-RNA/functions.R')
theme_set(plot_theme)

# IMPORT OG METADATA ------------------------------------------------------
meta.icicle <- fread('icicle/raw_data/icicle_samples.csv') %>% 
  rename(id = ID) %>% 
  mutate(sex = factor(Gender, levels = c(1, 2)),
         sex = recode(sex, `1` = 'Male', `2` = 'Female'))
meta.ppmi <- fread('ppmi/output/BL_iPD_HC_metadata.csv') %>% 
  mutate(id = as.character(id))

# IMPORT DDS --------------------------------------------------------------
dds.ppmi <- read_rds('linear/data/ppmi_DDS_prefilter.rds')
dds.icicle <- read_rds('linear/data/icicle_DDS_prefilter.rds')

# CHECK CLINCALLY REPORTED SEX --------------------------------------------
# RPS4Y1, KDM5D, DDX3Y and USP9Y (XIST not there after filtering)
sex_genes <- c('ENSG00000129824', 'ENSG00000012817', 'ENSG00000067048', 'ENSG00000114374')

sexPCA <- function(dds, meta){
  vsd <- vst(dds, blind = TRUE)
  vsd<- assay(vsd)
  sex_vsd <- vsd[row.names(vsd) %in% sex_genes, ]
  sex_pca <- prcomp(t(sex_vsd), center = TRUE, scale. = FALSE)
  sex_pca.df <- sex_pca$x %>% 
    as.data.frame() %>% 
    rownames_to_column('id') %>% 
    left_join(meta[, c('id', 'sex')], by = 'id')
  return(sex_pca.df)
}
# run PCA
sex.ppmi <- sexPCA(dds.ppmi, meta.ppmi) %>% mutate(study = 'PPMI')
sex.icicle <- sexPCA(dds.icicle, meta.icicle) %>% mutate(study = 'ICICLE-PD')
# plot PC1 and PC2
bind_rows(sex.ppmi, sex.icicle) %>% 
  mutate(study = factor(study, levels = c('PPMI', 'ICICLE-PD'))) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = sex, shape = sex)) +
  geom_point(alpha = 0.6) +
  labs(colour = 'Clinically reported sex',
       shape = "Clinically reported sex") +
  scale_color_brewer(palette = 'Set2') +
  facet_wrap(~ study, scales = 'free') +
  theme(legend.position = 'bottom')
ggsave('combinedLinCirc/output/figures/supp/sexQC.png', device = agg_png,
       height = 3, width = 5, dpi = 600)
ggsave('combinedLinCirc/output/figures/supp/sexQC.pdf',
       height = 3, width = 5, dpi = 600)

# which is it
sex.icicle %>% filter(PC1 > 0 & sex == 'Female')  # IN008 