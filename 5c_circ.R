# set up ------------------------------------------------------------------
library(here)
library(tximport)
library(easystats)
library(data.table)
library(scales)
library(broom)
library(stringr)
library(ggtext)
library(ggrastr)
library(ragg)
library(ggVennDiagram)
library(ggrepel)
library(ggpmisc)
library(ggpubr)
library(ggside)
library(cowplot)
library(patchwork)
library(colorspace)
library(RColorBrewer)
library(MetBrewer)
library(DT)
library(clipr)
library(pROC)
library(fst)
library(DT)
library(gtsummary)
library(tidyverse)
library(furrr)

set.seed(1997)

# import common functions
source(here("combinedLinCirc/PD-RNA/functions.R"))
# set theme
theme_set(plot_theme)


# import ------------------------------------------------------------------
## metadata
meta.ppmi <- read_rds(here("circRNA/data/ppmi_metadata.rds"))
meta.icicle <- read_rds(here("circRNA/data/icicle_metadata.rds"))


interesting_columns <- c("id", "study", "condition", "age_at_consent", "sex", "batch", "pct_usable_bases", "pct_intronic_bases", "pct_coding_bases", "median_cv_coverage", "total_sequences", "salmon_mapped", "mapped_reads", "sum_bsj", "unique_bsj", "unique_hostBSJ", "sum_bsj_perMapped", "unique_bsj_perMapped", "unique_hostBSJ_perMapped")
meta.bound <- bind_rows(
  meta.ppmi %>% select(all_of(interesting_columns)),
  meta.icicle %>% select(all_of(interesting_columns))
) %>%
  mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")))

# gene annotations
gene_anno <- read_delim(here("geneAnnotations_ensembl_v101.txt"))

# BSJ annotation info
junc_info <- read_fst(here("circRNA/data/bound_juncInfo.fst"))

# BSJ results
bsj_results.ppmi <- read_csv(here("circRNA/output/ppmi_BSJresults.csv")) %>%
  mutate(study = "PPMI", type = "BSJ") %>% 
  rename(gene_symbol = gene_name)
bsj_results.icicle <- read_csv(here("circRNA/output/icicle_BSJresults.csv")) %>%
  mutate(study = "ICICLE-PD", type = "BSJ") %>% 
  rename(gene_symbol = gene_name)

bsj_results.merged <- full_join(bsj_results.ppmi, bsj_results.icicle,
                                by = c("coord_id", "gene_id", "gene_symbol", "type"),
                                suffix = c(".ppmi", ".icicle")
)
bsj_results.bind <- bind_rows(bsj_results.ppmi, bsj_results.icicle) %>%
  mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")))

# FSJ results
fsj_results.ppmi <- read_csv(here("circRNA/output/ppmi_FSJresults.csv")) %>%
  mutate(study = "PPMI", type = "FSJ") %>% 
  rename(gene_symbol = gene_name)
fsj_results.icicle <- read_csv(here("circRNA/output/icicle_FSJresults.csv")) %>%
  mutate(study = "ICICLE-PD", type = "FSJ") %>% 
  rename(gene_symbol = gene_name)

# import normalised counts
norm_linear_counts.ppmi <- read_csv(here("linear/data/ppmi_deseqFilteredNormalisedCounts.csv")) %>%
  rename(gene_id = ensembl) %>%
  pivot_longer(cols = !gene_id, names_to = "id", values_to = "norm_count") %>%
  mutate(study = "PPMI")
norm_linear_counts.icicle <- read_csv(here("linear/data/icicle_deseqFilteredNormalisedCounts.csv")) %>%
  rename(gene_id = ensembl) %>%
  pivot_longer(cols = !gene_id, names_to = "id", values_to = "norm_count") %>%
  mutate(study = "ICICLE-PD")

vst_linear_counts.ppmi <- read_csv(here("linear/data/ppmi_vstCounts.csv")) %>%
  rename(gene_id = ensembl) %>%
  pivot_longer(cols = !gene_id, names_to = "id", values_to = "vst_count") %>%
  mutate(study = "PPMI")
vst_linear_counts.icicle <- read_csv(here("linear/data/icicle_vstCounts.csv")) %>%
  rename(gene_id = ensembl) %>%
  pivot_longer(cols = !gene_id, names_to = "id", values_to = "vst_count") %>%
  mutate(study = "ICICLE-PD")

### vst counts
### add to junc info for ease
junc_info <- left_join(junc_info,
                       full_join(
                         bind_rows(
                           # PPMI BSJ VST
                           read_csv(here("circRNA/data/ppmi_vstBSJCounts.csv")) %>%
                             pivot_longer(cols = !coord_id, names_to = "id", values_to = "bsj_vst") %>%
                             mutate(study = "PPMI"),
                           # ICICLE-PD BSJ vst
                           read_csv(here("circRNA/data/icicle_vstBSJCounts.csv")) %>%
                             pivot_longer(cols = !coord_id, names_to = "id", values_to = "bsj_vst") %>%
                             mutate(study = "ICICLE-PD")
                         ),
                         bind_rows(
                           # PPMI FSJ vst
                           read_csv(here("circRNA/data/ppmi_vstFSJCounts.csv")) %>%
                             pivot_longer(cols = !coord_id, names_to = "id", values_to = "fsj_vst") %>%
                             mutate(study = "PPMI"),
                           # ICICLE-PD FSJ vst
                           read_csv(here("circRNA/data/icicle_vstFSJCounts.csv")) %>%
                             pivot_longer(cols = !coord_id, names_to = "id", values_to = "fsj_vst") %>%
                             mutate(study = "ICICLE-PD")
                         ),
                         by = c("coord_id", "id", "study")
                       ),
                       by = c("coord_id", "id", "study")
) %>% 
  rename(gene_symbol = gene_name)


# detection ---------------------------------------------------------------
# How many BSJs detected in each cohort?
map(cohorts, function(study) {
  junc_info %>%
    filter(study == {{ study }}) %>%
    pull(coord_id) %>%
    unique() %>%
    length()
})

# how many of the abundant BSJs overlap between cohorts?
round(length(intersect(bsj_results.ppmi$coord_id, bsj_results.icicle$coord_id)) / length(unique(c(bsj_results.ppmi$coord_id, bsj_results.icicle$coord_id)))*100, 2)

# abundant overenrichment -------------------------------------------------
# Import
abundant_bsj_enrich.ppmi <- read_rds(here("circRNA/output/ppmi_abundantBSJ_enrich.rds"))
abundant_bsj_enrich.ppmi <- abundant_bsj_enrich.ppmi@result

abundant_bsj_enrich.icicle <- read_rds(here("circRNA/output/icicle_abundantBSJ_enrich.rds"))
abundant_bsj_enrich.icicle <- abundant_bsj_enrich.icicle@result


# Any sig in PPMI?
filter(abundant_bsj_enrich.ppmi, p.adjust < 0.05)

# Any sig in ICICLE?
filter(abundant_bsj_enrich.icicle, p.adjust < 0.05)

# Output merged file for supp
abundant_bsj_enrich <- full_join(abundant_bsj_enrich.ppmi, abundant_bsj_enrich.icicle, by = c("ONTOLOGY", "ID", "Description"), suffix = c(".ppmi", ".icicle"))
# add a space infront of the ratios so excel doesn't change it to a date
abundant_bsj_enrich <- abundant_bsj_enrich %>% 
  mutate(GeneRatio.ppmi = paste0(" ", GeneRatio.ppmi),
         GeneRatio.icicle = paste0(" ", GeneRatio.icicle),
         BgRatio.ppmi = paste0(" ", BgRatio.ppmi),
         BgRatio.icicle = paste0(" ", BgRatio.icicle))
write_csv(abundant_bsj_enrich, "combinedLinCirc/output/tables/merged_abundantBSJ_go_enrich.csv")


# differential expression ------------------------------------------------
# How many significant in PPMI?
sig_bsj_results.ppmi <- bsj_results.ppmi %>%
  filter(
    padj < 0.05,
    log2FoldChange > 0.1 | log2FoldChange < -0.1
  ) %>%
  mutate(direction = case_when(
    log2FoldChange > 0.1 ~ "Increased",
    log2FoldChange < -0.1 ~ "Decreased",
    TRUE ~ "error"
  ))
nrow(sig_bsj_results.ppmi)

# Volcano plot
bsj_volcano.plot <- plotVolcano(bsj_results.bind, sig_bsj_results.ppmi) +
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  scale_colour_manual(values = c("Not DE" = "gray70", "DE" = "#486de8")) +
  guides(colour = guide_legend(title = "circRNA differential expression in PPMI")) +
  theme(legend.position = "top") +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5))
bsj_volcano.plot
ggsave(here("combinedLinCirc/output/figures/individual/circ_volcano.png"),
       height = 4, width = 6, dpi = 600, device = agg_png
)

### Replication in ICICLE-PD
# Do any replicate?
replicate_bsj <- sig_bsj_results.ppmi %>%
  left_join(bsj_results.icicle,
            by = c("coord_id", "gene_id", "gene_symbol"),
            suffix = c(".ppmi", ".icicle")
  ) %>%
  mutate(replicate_fdr = p.adjust(pvalue.icicle, method = "fdr")) %>% 
  select(-contains(c("study", "type", "Direction", "padj.icicle")))


# Plot fold changes of BSJs differentially expressed in PPMI
# log2 fold change and lfcSE to long format 
# do it separately and then bind rows together as I can't figure out how to do it in pivot
plot.replicate_bsj <- full_join(replicate_bsj %>% 
                                  select(contains(c("coord_id", "log2FoldChange"))) %>% 
                                  pivot_longer(cols = !coord_id, names_to = "study", values_to = "log2FoldChange") %>% 
                                  mutate(study = gsub(study, pattern = "log2FoldChange.", replacement = "")),
                                replicate_bsj %>% 
                                  select(contains(c("coord_id", "lfcSE"))) %>% 
                                  pivot_longer(cols = !coord_id, names_to = "study", values_to = "lfcSE") %>% 
                                  mutate(study = gsub(study, pattern = "lfcSE.", replacement = ""))
                                , by = c("coord_id", "study")) %>% 
  # add gene name from annoation
  left_join(unique(junc_info[, c("coord_id", "gene_symbol")]), by = "coord_id") %>% 
  # change cohort names to correct names + order for colours
  mutate(study = case_match(study, "ppmi" ~ "PPMI","icicle" ~ "ICICLE-PD"),
         study = factor(study, levels = c("ICICLE-PD", "PPMI")),
         coord_label = paste0(coord_id, "\n",
                              "(circ", gene_symbol, ")")) %>% 
  # calculate 95% CIs
  mutate(ci_lower = log2FoldChange - (qnorm(0.025) * lfcSE),
         ci_higher = log2FoldChange + (qnorm(0.025) * lfcSE)) %>% 
  ggplot(aes(x = log2FoldChange, y = coord_label, colour = study)) +
  geom_vline(xintercept = 0, alpha = 0.2, linetype = "dashed") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_higher), 
                  position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = study_colours, breaks = c("PPMI", "ICICLE-PD")) +
  xlim(-0.7, 0.7) +
  labs(x = "log<sub>2</sub>(Fold Change)",
       y = "BSJ position",
       colour = "Cohort") +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5))
plot.replicate_bsj


# Output combined BSJ results table for supp
bsj_results.merged %>% 
  left_join(replicate_bsj[, c("coord_id", "replicate_fdr")], by = "coord_id") %>% 
  select(-contains(c("study", "type"))) %>% 
  write_csv(here("combinedLinCirc/output/tables/circRNA_DE.csv"))

## GO enrichment
sig_bsj_enrich.ppmi <- read_rds(here("circRNA/output/ppmi_sigBSJ_enrich.rds"))
sig_bsj_enrich.ppmi@result %>% filter(p.adjust < 0.05)



# prev reported PD circRNAs ----------------------------------------------
# get previously identified PD BSJs
circatlas <- read_delim(here("circRNA/data/circAtlas_circRNA_bed.txt"))
circatlas$format <- paste0(circatlas$Chro, ":", circatlas$Start, "|", circatlas$End)
circatlas$my_format <- paste0(circatlas$Chro, ":", circatlas$Start, "-", circatlas$End, ":", circatlas$Stand)
circatlas <- circatlas[, c("format", "my_format")]
prev_circs <- read_csv(here("circRNA/data/raw/manuallyCollatedPDcircs.csv")) %>%
  rename(format = coords_GRCh38_circAtlas)
prev_circs <- left_join(prev_circs, circatlas, by = "format") %>%
  separate(my_format, into = c("chr", "coords", "strand"), sep = ":") %>%
  separate(coords, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start) - 1) %>%
  mutate(chr = gsub(chr, pattern = "chr", replacement = ""))
prev_circs$coord_id <- paste0(prev_circs$chr, ":", prev_circs$start, "-", prev_circs$end, ":", prev_circs$strand)
rm(circatlas)

# GET RESULTS OF PREV BSJS
prev_circs.table <- bsj_results.merged %>%
  filter(coord_id %in% prev_circs$coord_id) %>%
  full_join(prev_circs[, c("coord_id", "gene", "study", "direction")], by = "coord_id") %>%
  mutate(
    tissue = recode(study, "Ravanidis.etal2021" = "Blood", "Hanan.etal2020" = "Brain"),
    direction = recode(direction, "down" = "Decreased in PD", "up" = "Increased in PD")
  ) %>%
  arrange(pvalue.ppmi) %>%
  select(
    coord_id, gene, tissue, direction,
    baseMean.ppmi, log2FoldChange.ppmi, lfcSE.ppmi, stat.ppmi, pvalue.ppmi,
    baseMean.icicle, log2FoldChange.icicle, lfcSE.icicle, stat.icicle, pvalue.icicle
  ) %>%
  mutate(
    replicate_fdr.ppmi = p.adjust(pvalue.ppmi, method = "fdr"),
    relicate_fdr.icicle = p.adjust(pvalue.icicle, method = "fdr")
  )

# Attempting to replicate prev PD BSJs
prev_circs.table %>%
  write_csv(here("combinedLinCirc/output/tables/replicate_prev_PD_bsj.csv"))

df.prev_circs <- prev_circs.table[, c("coord_id", "direction")] %>%
  left_join(bsj_results.bind, by = "coord_id", multiple = "all") %>%
  mutate(
    direction = recode(direction,
                       "Decreased in PD" = "↓",
                       "Increased in PD" = "↑"
    ),
    coord_id = paste(direction, coord_id, sep = " ")
  )

# correlation of fold changes in both
prev_circs_cor <- df.prev_circs[, c("coord_id", "log2FoldChange", "study")] %>%
  pivot_wider(id_cols = coord_id, names_from = study, values_from = log2FoldChange)
prev_circs_cor <- cor.test(prev_circs_cor$PPMI, prev_circs_cor$`ICICLE-PD`, method = "pearson") %>% tidy()
prev_circs_cor

# plot
plot.prev_circs <- df.prev_circs %>%
  drop_na() %>% 
  mutate(study = factor(study, levels = c("ICICLE-PD", "PPMI")),
         # add on gene symbol undernath coord id
         coord_label = paste0(coord_id, "\n",
                              "(circ", gene_symbol, ")")) %>%
  ggplot(aes(x = log2FoldChange, y = coord_label, colour = study)) +
  #ggforestplot::geom_stripes(aes(colour = coord_id), odd = "white", even = "gray95") +
  geom_vline(xintercept = 0, alpha = 0.2, linetype = "dashed") +
  geom_pointrange(aes(xmin = log2FoldChange - (qnorm(0.025) * lfcSE), xmax = log2FoldChange + (qnorm(0.025) * lfcSE)), position = position_dodge(0.5), size = 0.3) +
  scale_colour_manual(values = study_colours, breaks = c("PPMI", "ICICLE-PD")) +
  xlim(-0.7, 0.7) +
  labs(
    x = "log<sub>2</sub>(Fold Change)",
    y = "BSJ position",
    colour = "Cohort"
  ) +
  theme(legend.position = "top") +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5))
plot.prev_circs
ggsave(here("combinedLinCirc/output/figures/individual/prev_pd_bsj.png"),
       height = 4, width = 5, dpi = 600, device = agg_png
)

# relate to previous work -------------------------------------------------
monogenicPD <- read_csv(here("combinedLinCirc/data/pdGenes_genomicsEnglandPanel.csv")) %>% 
  rename(gene_symbol = gene_name)

# monogenic
bsj_monogenicPD <- monogenicPD %>%
  inner_join(bsj_results.merged, by = "gene_symbol") %>%
  select(
    gene_symbol, gene_id, coord_id, baseMean.ppmi, baseMean.icicle,
    log2FoldChange.ppmi, log2FoldChange.icicle, pvalue.ppmi, pvalue.icicle
  ) %>%
  mutate(
    replicate_fdr.ppmi = p.adjust(pvalue.ppmi, method = "fdr"),
    replicate_fdr.icicle = p.adjust(pvalue.icicle, method = "fdr")
  )
bsj_monogenicPD %>%
  write_csv(here("combinedLinCirc/output/bsj_monogenicPD_overlap.csv"))

# How many genes host BSJs?
length(unique(bsj_monogenicPD$gene_symbol))

# GWAS
gwas_genes <- read_csv(here("combinedLinCirc/data/nalls_et_al2019_META5_genes.csv")) %>%
  janitor::clean_names() %>%
  rename(gene_symbol = nearest_gene)

length(unique(gwas_genes$gene_symbol))

bsj_gwas <- gwas_genes %>%
  inner_join(bsj_results.merged, by = "gene_symbol", multiple = "all") %>%
  distinct() %>%
  select(
    gene_symbol, coord_id, baseMean.ppmi, baseMean.icicle,
    log2FoldChange.ppmi, log2FoldChange.icicle, pvalue.ppmi, pvalue.icicle
  ) %>%
  mutate(
    new_pval.ppmi = p.adjust(pvalue.ppmi, method = "fdr"),
    new_pval.icicle = p.adjust(pvalue.icicle, method = "fdr"),
  ) %>%
  arrange(pvalue.ppmi)
bsj_gwas %>% datatable()
length(unique(bsj_gwas$gene_symbol)) # number of BSJ-hosting GWAS genes
bsj_gwas %>%
  write_csv(here("combinedLinCirc/output/bsj_GWAS_overlap.csv"))


# figure panel ------------------------------------------------------------
(bsj_volcano.plot | (plot.replicate_bsj / plot.prev_circs) + plot_layout(guides = 'collect', heights = c(1.5, 2))) + plot_layout(widths = c(2, 1)) & theme(legend.position = 'top')
ggsave(here("combinedLinCirc/output/figures/panels/circ_DE.svg"),
       height = 8, width = 9)
