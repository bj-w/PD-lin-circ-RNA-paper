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



### linear
linear_results.ppmi <- read_csv(here("linear/output/ppmi_deseqResults.csv")) %>%
  mutate(study = "PPMI", type = "linear") %>%
  rename(gene_id = ensembl) %>%
  arrange(pvalue) %>% 
  select(-gene_name) %>% 
  left_join(gene_anno, by = "gene_id")

linear_results.icicle <- read_csv(here("linear/output/icicle_deseqResults.csv")) %>%
  mutate(study = "ICICLE-PD", type = "linear") %>%
  rename(gene_id = ensembl) %>%
  arrange(pvalue) %>% 
  select(-gene_name) %>% 
  left_join(gene_anno, by = "gene_id")
linear_results.bind <- bind_rows(linear_results.ppmi, linear_results.icicle) %>%
  mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")))

linear_results.merged <- full_join(linear_results.ppmi, linear_results.icicle, by = c("gene_id", "gene_symbol"), suffix = c(".ppmi", ".icicle"))



# DE ----------------------------------------------------------------------

# how many of the tested genes overlap between PPMI and ICICLE?
round(length(intersect(linear_results.ppmi$gene_id, linear_results.icicle$gene_id)) / length(unique(c(linear_results.ppmi$gene_id, linear_results.icicle$gene_id)))*100, 2)


sig_linear_results.ppmi <- linear_results.ppmi %>%
  filter(
    padj < 0.05,
    log2FoldChange > 0.1 | log2FoldChange < -0.1
  )
sig_linear_results.ppmi %>% datatable()


# How many sig? How many increased + decreased? (TRUE = Increased)
nrow(sig_linear_results.ppmi)
table(sig_linear_results.ppmi$log2FoldChange > 0.1)

# previously reported PPMI DEGs
ppmi_degs <- read_csv("combinedLinCirc/data/PPMI_DEGs.csv") %>% 
  # remove gene version number from gene_id
  separate(gene_id, into = "gene_id", sep = "\\.") %>%
  filter(adj.P.Val < 0.05,
         logFC < -0.1 | logFC > 0.1)
sig_linear_results.ppmi %>% filter(gene_id %in% ppmi_degs$gene_id) %>% view()


# Volcano plot highlighting significant PPMI RNAs in both cohorts
linear_volcano.plot <- plotVolcano(linear_results.bind, sig_linear_results.ppmi) +
  scale_x_continuous(limits = c(-2.1, 2.1)) +
  scale_colour_manual(values = c("Not DE" = "gray70", "DE" = "#c55305")) +
  guides(colour = guide_legend(
    title = "Gene differential expression in PPMI")) +
  theme(legend.position = "top") +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5))
linear_volcano.plot
ggsave(here("combinedLinCirc/output/figures/individual/linear_volcano.png"),
       height = 4, width = 6, dpi = 600, device = agg_png
)


# Genes that replicate in ICICLE-PD
replicate_linear <- linear_results.merged %>%
  filter(
    padj.ppmi < 0.05,
    log2FoldChange.ppmi > 0.1 | log2FoldChange.ppmi < -0.1
  ) %>%
  mutate(replicate_fdr = p.adjust(pvalue.icicle, method = "fdr")) %>%
  arrange(replicate_fdr) %>%
  mutate(gene_symbol = as_factor(gene_symbol)) %>% 
  select(-contains(c('padj.icicle', 'entrez', 'study', 'type')))
replicate_linear %>% datatable()
write_csv(replicate_linear, here("combinedLinCirc/output/tables/replicateLinearDE.csv"))

# How many replicate?
sig_replicate_linear <- replicate_linear %>%
  filter(
    replicate_fdr < 0.05,
    log2FoldChange.icicle > 0.1 | log2FoldChange.icicle < -0.1
  )
sig_replicate_linear %>% datatable()

# Plot replicated genes
plot.linear_replicated <- linear_results.bind %>%
  filter(gene_id %in% sig_replicate_linear$gene_id) %>%
  mutate(study = factor(study, levels = c("ICICLE-PD", "PPMI"))) %>%
  ggplot(aes(log2FoldChange, gene_symbol,
             colour = study
  )) +
  geom_pointrange(aes(xmin = log2FoldChange - (qnorm(0.025) * lfcSE), xmax = log2FoldChange + (qnorm(0.025) * lfcSE)),
                  position = position_dodge(width = 0.3)
  ) +
  scale_colour_manual(values = study_colours, breaks = c("PPMI", "ICICLE-PD")) +
  labs(
    x = "log<sub>2</sub>(Fold Change)",
    y = "Gene symbol",
    colour = "Cohort"
  ) +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "top"
    # legend.justification = c(1, 0), legend.position = c(1, 0),
    # legend.background = element_blank()
  ) +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5))
plot.linear_replicated
ggsave(here("combinedLinCirc/output/figures/individual/replicated_linear.png"),
       height = 5, width = 3, dpi = 600, device = agg_png
)

# Export linear results tables for supp
linear_results.merged %>% 
  left_join(replicate_linear[,c("gene_id", "replicate_fdr")], by = "gene_id") %>% 
  select(-contains(c("entrez", "study", "type"))) %>% 
  write_csv(here("combinedLinCirc/output/tables/linearRNA_DE.csv"))


# previously reported differentially expressed genes in blood? (RN --------
# Import data
prev_DE_genes_blood_rnaseq <- read_csv(here("linear/data/prev_DE_genes_blood_rnaseq.csv"))

# How many unique genes?
length(unique(prev_DE_genes_blood_rnaseq$gene_symbol))

# Overlap with DE results
table.prev_DE_genes_blood_rnaseq <- prev_DE_genes_blood_rnaseq %>% 
  left_join(linear_results.merged, by = "gene_id") %>%
  mutate(
    fdr.ppmi = p.adjust(pvalue.ppmi, method = "fdr"),
    fdr.icicle = p.adjust(pvalue.icicle, method = 'fdr')) %>% 
  select(-contains(c('padj', 'entrez', 'study', 'type'))) %>% 
  rename(reported_gene_symbol = gene_symbol.x,
         gene_symbol = gene_symbol.y)
table.prev_DE_genes_blood_rnaseq %>% datatable()

write_csv(table.prev_DE_genes_blood_rnaseq, here("combinedLinCirc/output/tables/prev_de_genes_blood_rnaseq.csv"))


# Make a volcano plot for the above genes
plot.prev_DE_genes_blood_rnaseq <- table.prev_DE_genes_blood_rnaseq %>% 
  select(contains(c("gene_id", "log2FoldChange", "pvalue", "fdr", "gene_symbol"))) %>% 
  drop_na(gene_id) %>% 
  select(-reported_gene_symbol) %>% 
  pivot_longer(cols = -c(gene_id, gene_symbol), names_to = c("column", "study"), names_sep = "\\.", values_to = "values") %>%
  distinct() %>% 
  pivot_wider(id_cols = c("gene_id", "study", "gene_symbol"), names_from = column, values_from = values) %>% 
  mutate(study = case_match(study,
                            "ppmi" ~ "PPMI",
                            "icicle" ~ "ICICLE-PD"),
         study = factor(study, levels = c("PPMI", "ICICLE-PD")))
plot.prev_DE_genes_blood_rnaseq <- plot.prev_DE_genes_blood_rnaseq %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue), label = gene_symbol)) +
  geom_point(aes(colour = study), alpha = 0.5) +
  geom_text_repel(
    data = filter(plot.prev_DE_genes_blood_rnaseq, fdr < 0.05),
    size = 3, min.segment.length = 0
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(limits = c(-1.1, 1.1)) +
  scale_colour_manual(values = study_colours) +
  labs(
    x = "log<sub>2</sub>(Fold Change)",
    y = "-log<sub>10</sub>(<i>P</i>)",
    colour = "Cohort"
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", alpha = 0.5) +
  facet_wrap(~study, scales = "free")
plot.prev_DE_genes_blood_rnaseq

# GWAS risk loci ----------------------------------------------------------
# Taken from [Nalls et al 2019](https://www.sciencedirect.com/science/article/pii/S1474442219303205?via%3Dihub)
gwas_genes <- read_csv(here("combinedLinCirc/data/nalls_et_al2019_META5_genes.csv")) %>%
  janitor::clean_names() %>%
  rename(gene_symbol = nearest_gene)

length(unique(gwas_genes$gene_symbol))


# Overlap
# DE in PPMI
gwas_linear.ppmi <- gwas_genes %>% 
  select(gene_id) %>% 
  distinct() %>% 
  left_join(linear_results.ppmi, by = "gene_id") %>% 
  mutate(fdr.ppmi = p.adjust(pvalue, method = "fdr"))

# replicate significant ones in ICICLE-PD
gwas_linear.icicle <- gwas_linear.ppmi %>% 
  filter(fdr.ppmi < 0.05,
         log2FoldChange > 0.1 | log2FoldChange < -0.1) %>% 
  select(gene_id) %>% 
  left_join(linear_results.icicle, by = "gene_id") %>% 
  mutate(fdr.icicle = p.adjust(pvalue, method = "fdr"))

# create combined table
table.linear_gwas <- gwas_genes %>% 
  select(snp, gene_symbol, gene_id) %>% 
  left_join(linear_results.merged, by = "gene_id") %>% 
  # add on FDR from each cohort
  left_join(gwas_linear.ppmi[, c("gene_id", "fdr.ppmi")], by = "gene_id") %>% 
  left_join(gwas_linear.icicle[, c("gene_id", "fdr.icicle")], by = "gene_id") %>% 
  select(-contains(c("padj", "study", "type", "gene_biotype", "entrez"))) %>% 
  select(-"gene_symbol.y") %>% 
  rename(gene_symbol = gene_symbol.x)
table.linear_gwas %>% datatable()
table.linear_gwas %>% 
  write_csv(here("combinedLinCirc/output/tables/linear_gwas_overlap.csv"))


# Plot
table.linear_gwas %>% glimpse()
plot.linear_gwas <- table.linear_gwas %>% 
  drop_na(gene_id) %>% 
  select(contains(c("gene_id", "log2FoldChange", "pvalue", "fdr"))) %>% 
  pivot_longer(cols = !gene_id, names_to = c("column", "study"), names_sep = "\\.", values_to = "values") %>% 
  distinct() %>% 
  pivot_wider(id_cols = c(gene_id, study), names_from = column, values_from = values) %>% 
  mutate(study = case_match(study,
                            "ppmi" ~ "PPMI",
                            "icicle" ~ "ICICLE-PD"),
         study = factor(study, levels = c("PPMI", "ICICLE-PD"))) %>% 
  # add on gene symbols
  left_join(unique(linear_results.bind[, c("gene_id", "gene_symbol")]), by = "gene_id")

# plot
plot.linear_gwas <- plot.linear_gwas %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue), label = gene_symbol)) +
  geom_point(aes(colour = study), alpha = 0.5) +
  geom_text_repel(data = filter(plot.linear_gwas,
                                fdr < 0.05,
                                log2FoldChange > 0.1 | log2FoldChange < -0.1), size = 3, min.segment.length = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_colour_manual(values = study_colours) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", alpha = 0.5) +
  labs(
    x = "log<sub>2</sub>(Fold change)",
    y = "-log<sub>10</sub>(<i>P</i>)", colour = "Cohort"
  ) +
  facet_wrap(~study)
plot.linear_gwas


# monogenic pd genes ------------------------------------------------------
# Import
monogenicPD <- read_csv(here("combinedLinCirc/data/pdGenes_genomicsEnglandPanel.csv")) %>% 
  rename(gene_symbol = gene_name)

# Overlap

# PPMI DE
linear_monogenic.ppmi <- linear_results.ppmi %>% 
  filter(gene_symbol %in% monogenicPD$gene_symbol) %>% 
  mutate(fdr.ppmi = p.adjust(pvalue, method = "fdr"))

# replicate in ICICLE-PD
linear_monogenic.icicle <- linear_monogenic.ppmi %>% 
  filter(fdr.ppmi < 0.05,
         log2FoldChange > 0.1 | log2FoldChange < -0.1) %>% 
  select(gene_id) %>% 
  left_join(linear_results.icicle, by = "gene_id") %>% 
  mutate(fdr.icicle = p.adjust(pvalue, method = "fdr"))

# combine into one table
table.linear_monogenic <- linear_results.merged %>% 
  filter(gene_symbol %in% monogenicPD$gene_symbol) %>% 
  # add on cohort FDR values
  left_join(linear_monogenic.ppmi[, c("gene_id", "fdr.ppmi")], by = "gene_id") %>% 
  left_join(linear_monogenic.icicle[, c("gene_id", "fdr.icicle")], by = "gene_id") %>%
  select(-contains(c("entrez", "study", "type", "gene_biotype", "padj"))) %>% 
  relocate(gene_symbol)
write_csv(table.linear_monogenic, here("combinedLinCirc/output/tables/linear_monogenicPD.csv"))


linear_monogenicPD <- monogenicPD %>%
  inner_join(linear_results.merged, by = "gene_symbol") %>%
  select(
    gene_symbol, gene_id, baseMean.ppmi, baseMean.icicle,
    log2FoldChange.ppmi, log2FoldChange.icicle, pvalue.ppmi, pvalue.icicle
  ) %>%
  mutate(
    replicate_fdr.ppmi = p.adjust(pvalue.ppmi, method = "fdr"),
    replicate_fdr.icicle = p.adjust(pvalue.icicle, method = "fdr")
  )
linear_monogenicPD %>% datatable()
linear_monogenicPD %>% write_csv(here("combinedLinCirc/output/linear_monogenicPD_overlap.csv"))

# GSEA --------------------------------------------------------------------

### Gene Ontologies

# Import GSEA results
linear_gsea_go.ppmi <- readRDS(here("linear/output/ppmi_GO_gsea.rds"))
linear_gsea_go.icicle <- readRDS(here("linear/output/icicle_GO_gsea.rds"))

# Overlap results between cohorts
# PPMI results
linear_gsea_go_fdr.ppmi <- linear_gsea_go.ppmi@result %>% 
  group_by(ONTOLOGY) %>% 
  mutate(fdr.ppmi = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

# replicate in ICICLE-PD
linear_gsea_go_fdr.icicle <- linear_gsea_go_fdr.ppmi %>% 
  filter(fdr.ppmi < 0.05) %>% 
  select(ID) %>% 
  left_join(linear_gsea_go.icicle@result, by = "ID") %>% 
  group_by(ONTOLOGY) %>% 
  mutate(fdr.icicle = p.adjust(pvalue, method = "fdr")) %>% 
  ungroup()

# combine into one table
table.linear_gsea_go <- full_join(linear_gsea_go.ppmi@result,
                                  linear_gsea_go.icicle@result,
                                  by = c("ONTOLOGY", "ID", "Description"),
                                  suffix = c(".ppmi", ".icicle")) %>% 
  select(contains(c("ONTOLOGY", "ID", "Description", "NES", "pvalue"))) %>% 
  # add on fdr columns for each cohort
  left_join(linear_gsea_go_fdr.ppmi[, c("ID", "fdr.ppmi")], by = "ID") %>% 
  left_join(linear_gsea_go_fdr.icicle[, c("ID", "fdr.icicle")], by = "ID") %>% 
  # add on column saying whether the NES direction agrees
  mutate(direction = case_when(
    NES.ppmi > 0 & NES.icicle > 0 ~ "Agree",
    NES.ppmi < 0 & NES.icicle < 0 ~ "Agree",
    TRUE ~ "Disagree"
  ))
table.linear_gsea_go %>% datatable()
write_csv(table.linear_gsea_go, here("combinedLinCirc/output/tables/linear_GSEA_GO.csv"))


# How many sig in each ontology (and agree on direction)
table.linear_gsea_go %>% 
  filter(direction == "Agree",
         fdr.ppmi < 0.05,
         fdr.icicle < 0.05) %>% 
  group_by(ONTOLOGY) %>% 
  count()

# Plot top ranked ontologies
plotGSEA <- function(df) {
  df %>%
    pivot_longer(cols = c("NES.ppmi", "NES.icicle"), names_to = "study", values_to = "NES") %>%
    mutate(study = recode_factor(study,
                                 "NES.ppmi" = "PPMI",
                                 "NES.icicle" = "ICICLE-PD"
    )) %>%
    ggplot(aes(x = study, y = fct_reorder(Description, pvalue.ppmi, .desc = TRUE), fill = NES)) +
    geom_tile(colour = "black") +
    scale_fill_gradient2() +
    scale_y_discrete(expand = c(0, 0), labels = wrap_format(40)) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(x = "", y = "", fill = "NES") +
    facet_wrap(~ONTOLOGY, scales = "free", ncol = 1) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
}


plot.linear_gsea_go <- table.linear_gsea_go %>% 
  as_tibble() %>% 
  filter(direction == "Agree",
         fdr.ppmi < 0.05,
         fdr.icicle < 0.05) %>% 
  mutate(ONTOLOGY = case_match(ONTOLOGY,
                               "BP" ~ "Biological Process",
                               "CC" ~ "Cellular Component",
                               "MF" ~ "Molecular Function")) %>% 
  group_by(ONTOLOGY) %>%
  slice_min(order_by = pvalue.ppmi, n = 10) %>%
  plotGSEA()
plot.linear_gsea_go



### KEGG pathways
# Import
linear_gsea_kegg.ppmi <- read_rds(here("linear/output/ppmi_KEGG_gsea.rds"))
linear_gsea_kegg.icicle <- read_rds(here("linear/output/icicle_KEGG_gsea.rds"))

# Overlap results
# PPMI sig 
linear_gsea_kegg_fdr.ppmi <- linear_gsea_kegg.ppmi@result %>% 
  as_tibble() %>% 
  mutate(fdr.ppmi = p.adjust(pvalue, method = "fdr"))
# replicate sig PPMI in ICICLE-PD
linear_gsea_kegg_fdr.icicle <- linear_gsea_kegg.ppmi@result %>% 
  as_tibble() %>% 
  filter(p.adjust < 0.05) %>% 
  select(ID, Description) %>% 
  left_join(linear_gsea_kegg.icicle@result, by = c("ID", "Description")) %>% 
  mutate(fdr.icicle = p.adjust(pvalue, method = "fdr"))

# combined table
table.linear_gsea_kegg <- full_join(linear_gsea_kegg.ppmi@result,
                                    linear_gsea_kegg.icicle@result,
                                    by = c("ID", "Description"),
                                    suffix = c(".ppmi", ".icicle")) %>% 
  as_tibble() %>% 
  select(contains(c("ID", "Description", "NES", "pvalue"))) %>% 
  left_join(linear_gsea_kegg_fdr.ppmi[, c("ID", "Description", "fdr.ppmi")], by = c("ID", "Description")) %>% 
  left_join(linear_gsea_kegg_fdr.icicle[, c("ID", "Description", "fdr.icicle")], by = c("ID", "Description")) %>%
  mutate(direction = case_when(
    NES.ppmi > 0 & NES.icicle > 0 ~ "Agree",
    NES.ppmi < 0 & NES.icicle < 0 ~ "Agree",
    TRUE ~ "Disagree"
  ))
write_csv(table.linear_gsea_kegg, here("combinedLinCirc/output/tables/linear_GSEA_KEGG.csv"))



# # Plot just ribosome GSEA plot
# ribo.ppmi <- enrichplot::gsearank(linear_gsea_kegg.ppmi, geneSetID = "hsa03010", title = "PPMI - Ribosome (hsa03010)") + plot_theme
# ribo.icicle <- enrichplot::gsearank(linear_gsea_kegg.icicle, geneSetID = "hsa03010", title = "ICICLE-PD - Ribosome (hsa03010)") + plot_theme
# plot.ribo_kegg <- plot_grid(ribo.ppmi, ribo.icicle)



## Figure panel
((((linear_volcano.plot | plot.linear_replicated) + plot_layout(widths = c(3, 1))) /
    ((plot.prev_DE_genes_blood_rnaseq / plot.linear_gwas) + plot_layout(guides = 'collect'))) | plot.linear_gsea_go) + plot_layout(widths = c(6, 1))
ggsave(here("combinedLinCirc/output/figures/panels/linear_DE.svg"),
       height = 12, width = 12)
