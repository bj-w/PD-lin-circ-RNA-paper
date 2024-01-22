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

# linear
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

# import linear normalised counts
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




# imbalance binomial test ---------------------------------------------------------------------
# function
imbalanceTest <- function(results) {
  filt_results <- results %>%
    filter(log2FoldChange > 0.1 | log2FoldChange < -0.1) %>%
    mutate(direction = case_when(
      log2FoldChange > 0.1 ~ "increased",
      log2FoldChange < -0.1 ~ "decreased",
      TRUE ~ "error"
    ),
    direction = factor(direction, levels = c("increased", "decreased")))
  numDirection <- table(filt_results$direction)
  print(numDirection)
  tidy(binom.test(x = numDirection))
}


# BSJ imbalance
imbalanceTest(bsj_results.ppmi)
imbalanceTest(bsj_results.icicle)


# FSJ imbalance
imbalanceTest(fsj_results.ppmi)
imbalanceTest(fsj_results.icicle)


# Linear RNA imbalance
imbalanceTest(linear_results.ppmi)
imbalanceTest(linear_results.icicle)


# Linear RNA restricted to those that host BSJs
linear_hostBSJ_results.ppmi <- linear_results.ppmi %>%
  filter(gene_id %in% bsj_results.ppmi$gene_id) %>%
  mutate(type = "linearHostBSJ")
linear_hostBSJ_results.icicle <- linear_results.icicle %>%
  filter(gene_id %in% bsj_results.icicle$gene_id) %>%
  mutate(type = "linearHostBSJ")
imbalanceTest(linear_hostBSJ_results.ppmi)
imbalanceTest(linear_hostBSJ_results.icicle)


# # imbalance volcano plots ---------------------------------------------------------------------

# format results df - rename coord_id/gene_id, select only the relevant columns
imbalance_volcanos <- map(
  list(
    bsj_results.ppmi, bsj_results.icicle,
    fsj_results.ppmi, fsj_results.icicle,
    linear_results.ppmi, linear_results.icicle,
    linear_hostBSJ_results.ppmi, linear_hostBSJ_results.icicle
  ),
  ~ rename(.x, "feature" = 1) %>%
    select(feature, log2FoldChange, pvalue, type, study)
) %>%
  list_rbind() %>%
  # add study levels for plotting
  mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")))
# calculate binomial test stats for each RNA in each cohort
binom_stats <- map(cohorts, function(study) {
  map2(
    .x = unique(imbalance_volcanos$type),
    .y = study,
    .f = function(type, study) {
      filt_df <- imbalance_volcanos %>%
        filter(
          type == {{ type }},
          study == {{ study }}
        )
      output <- imbalanceTest(filt_df) %>%
        mutate(
          type = {{ type }},
          study = {{ study }}
        )
      return(output)
    }
  ) %>% list_rbind()
}) %>%
  list_rbind() %>%
  mutate(
    p.value = if_else(p.value < 0.001,
                      expSup(p.value),
                      format.pval(p.value, digits = 1)
    ),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2)
  ) %>%
  mutate(
    type = recode(type,
                  "linear" = "Gene",
                  "linearHostBSJ" = "BSJ hosts (Gene)"
    )
  )
imbalance_volcanos <- imbalance_volcanos %>%
  mutate(
    type = recode(type,
                  "linear" = "Gene",
                  "linearHostBSJ" = "BSJ hosts (Gene)"
    ),
    fc_colour = case_when(
      type == "BSJ" & (log2FoldChange > 0.1 | log2FoldChange < -0.1) ~ "BSJ",
      type == "FSJ" & (log2FoldChange > 0.1 | log2FoldChange < -0.1) ~ "FSJ",
      type == "Gene" & (log2FoldChange > 0.1 | log2FoldChange < -0.1) ~ "Gene",
      type == "BSJ hosts (Gene)" & (log2FoldChange > 0.1 | log2FoldChange < -0.1) ~ "BSJ hosts (Gene)",
      TRUE ~ "noColour"
    ),
  )
plot.imbalance_volcanos <- imbalance_volcanos %>%
  mutate(type = factor(type, levels = c("BSJ", "FSJ", "BSJ hosts (Gene)", "Gene"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(colour = fc_colour), alpha = 0.4) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", alpha = 0.5) +
  geom_richtext(
    data = binom_stats, aes(x = Inf, y = Inf, label = paste0(
      "<i>P</i> = ", p.value, "  <br>",
      "N = ", parameter, "  "
    )),
    vjust = 1, hjust = "right", size = 3, fill = NA, label.color = NA
  ) +
  geom_xsidedensity(data = filter(imbalance_volcanos, log2FoldChange > 0.1 | log2FoldChange < -0.1), aes(y = stat(density), fill = fc_colour), colour = NA, alpha = 0.7) +
  theme_ggside_void() +
  scale_colour_manual(values = c(rna_colours, "noColor" = "gray90")) +
  scale_fill_manual(values = c(rna_colours)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_grid(factor(type, levels = c("BSJ", "FSJ", "BSJ hosts (Gene)", "Gene")) ~ factor(study, levels = c("PPMI", "ICICLE-PD")), scales = "free_y", labeller = labeller(type = label_wrap_gen(10))) +
  labs(
    x = "log<sub>2</sub>(Fold Change)",
    y = "-log<sub>10</sub>(<i>P</i>)",
    colour = "RNA type",
    fill = "RNA type"
  ) +
  theme(legend.position = "top")
plot.imbalance_volcanos
ggsave(here("combinedLinCirc/output/figures/individual/imbalance_volcanoes.png"),
       height = 7, width = 6, dpi = 600, device = agg_png)

# export df for supp
imbalance_volcanos %>%
  select(-fc_colour) %>%
  write_csv(here("combinedLinCirc/output/tables/imbalanceVolcanoes.csv"))


# imbalance forest plot ---------------------------------------------------
imbalance_forest <- map(
    list(
      bsj_results.ppmi, bsj_results.icicle,
      fsj_results.ppmi, fsj_results.icicle,
      linear_results.ppmi, linear_results.icicle,
      linear_hostBSJ_results.ppmi, linear_hostBSJ_results.icicle
    ),
    ~ rename(.x, "feature" = 1) %>%
      select(feature, log2FoldChange, pvalue, type, study)
  ) %>%
    list_rbind() %>%
    # add study levels for plotting
    mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")))
  # calculate binomial test stats for each RNA in each cohort
  binom_stats <- map(cohorts, function(study) {
    map2(
      .x = unique(imbalance_forest$type),
      .y = study,
      .f = function(type, study) {
        filt_df <- imbalance_forest %>%
          filter(
            type == {{ type }},
            study == {{ study }}
          )
        output <- imbalanceTest(filt_df) %>%
          mutate(
            type = {{ type }},
            study = {{ study }}
          )
        return(output)
      }
    ) %>% 
      list_rbind() %>% 
      mutate(bonf = p.adjust(p.value, method = "bonferroni"))
  }) %>%
    list_rbind() %>%
    mutate(
      bonf_label = if_else(bonf < 0.001,
                        expSup(bonf),
                        format.pval(bonf, digits = 1)
      ),
      estimate = round(estimate, 2),
    ) %>%
    mutate(
      type = recode(type,
                    "linear" = "Gene",
                    "linearHostBSJ" = "BSJ hosts (Gene)"
      )
    )
  
plot.imbalance_forest <- binom_stats %>% 
    mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")),
           type = factor(type, levels = c("BSJ", "FSJ", "BSJ hosts (Gene)", "Gene"))) %>% 
    ggplot(aes(x = estimate, y = type, colour = type, label = bonf_label)) +
    geom_vline(xintercept = 0.5, colour = "gray80", linetype = "dashed") +
    # geom_point() +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
    geom_richtext(nudge_y = 0.3, size = 3, fill = NA, label.color = NA) +
    #geom_text_repel(colour = "black", nudge_y = 0.2, size = 3) +
    scale_y_discrete(limits=rev, labels = label_wrap(width = 10)) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = rna_colours) +
    labs(x = "Proportion of loci increased in PD vs Controls", y = "RNA Type", colour = "RNA type") +
    facet_wrap(~ study , ncol = 2) +
    theme(legend.position = "bottom",
          strip.background = element_blank())
plot.imbalance_forest


# differences in mean expression --------------------------------------------------------------
# Function for calculating mean difference and plotting
meanNormCountDiff <- function(ppmiCountMatrix, icicleCountMatrix, ppmiMeta, icicleMeta) {
  ppmi <- ppmiCountMatrix %>%
    rename(rna = 1) %>%
    pivot_longer(cols = !rna, names_to = "id", values_to = "norm_count") %>%
    left_join(ppmiMeta[, c("id", "condition")], by = "id") %>%
    mutate(study = "PPMI")
  icicle <- icicleCountMatrix %>%
    rename(rna = 1) %>%
    pivot_longer(cols = !rna, names_to = "id", values_to = "norm_count") %>%
    left_join(icicleMeta[, c("id", "condition")], by = "id") %>%
    mutate(study = "ICICLE-PD")
  # bind study dfs together
  bound <- bind_rows(ppmi, icicle)
  bound <- bound %>%
    # calculate mean counts for each gene/rna for each study
    group_by(study, rna, condition) %>%
    summarise(mean_norm_count = mean(norm_count)) %>%
    # set factor levels for plotting
    mutate(
      condition = factor(condition, levels = c("PD", "Control")),
      study = factor(study, levels = c("PPMI", "ICICLE-PD"))
    )
  # calculate p values for plotting
  pvals <- map(cohorts, function(study) {
    wilcox_df <- tidy(wilcox.test(mean_norm_count ~ condition, data = filter(bound, study == {{ study }}, conf.int = TRUE))) %>%
      mutate(
        study = {{ study }},
        p.value = if_else(p.value < 0.001,
                          expSup(p.value),
                          format.pval(p.value, digits = 1)
        )
      )
  }) %>% list_rbind()
  # plot
  bound %>%
    ggplot(aes(x = condition, y = mean_norm_count + 1)) +
    geom_violin(aes(fill = condition), colour = NA) +
    geom_boxplot(outlier.shape = NA, width = 0.2) +
    geom_richtext(
      data = pvals, aes(x = 1.5, y = Inf, label = paste0(
        "<i>P</i> = ", p.value
      )),
      vjust = 1, size = 3, fill = NA, label.color = NA
    ) +
    scale_y_log10(expand = expansion(mult = c(0.1, 0.2))) +
    scale_fill_manual(values = condition_colours) +
    labs(x = "Study group", y = "Mean normalised count +1", fill = "Study group") +
    facet_wrap(~ factor(study, levels = c("PPMI", "ICICLE-PD")))
}

# BSJ
plot.bsj_counts_condition <- meanNormCountDiff(
  ppmiCountMatrix = read_csv(here("circRNA/data/ppmi_deseqBSJFilteredNormalisedCounts.csv")),
  icicleCountMatrix = read_csv(here("circRNA/data/icicle_deseqBSJFilteredNormalisedCounts.csv")),
  ppmiMeta = meta.ppmi, icicleMeta = meta.icicle
) +
  ylab("Mean BSJ expression +1")
plot.bsj_counts_condition

# FSJ
plot.fsj_counts_condition <- meanNormCountDiff(
  ppmiCountMatrix = read_csv(here("circRNA/data/ppmi_deseqFSJFilteredNormalisedCounts.csv")),
  icicleCountMatrix = read_csv(here("circRNA/data/icicle_deseqFSJFilteredNormalisedCounts.csv")),
  ppmiMeta = meta.ppmi, icicleMeta = meta.icicle
) +
  ggtitle("FSJ expression")
plot.fsj_counts_condition

# Gene BSJ hosts
plot.bsjHosts_counts_condition <- meanNormCountDiff(
  ppmiCountMatrix = norm_linear_counts.ppmi %>%
    select(-study) %>%
    filter(gene_id %in% bsj_results.ppmi$gene_id) %>%
    pivot_wider(id_cols = "gene_id", names_from = id, values_from = norm_count),
  icicleCountMatrix = norm_linear_counts.icicle %>%
    select(-study) %>%
    filter(gene_id %in% bsj_results.icicle$gene_id) %>%
    pivot_wider(id_cols = "gene_id", names_from = id, values_from = norm_count),
  ppmiMeta = meta.ppmi, icicleMeta = meta.icicle
) +
  ggtitle("circRNA host gene expression")
plot.bsjHosts_counts_condition

# All genes
plot.gene_counts_condition <- meanNormCountDiff(
  ppmiCountMatrix = norm_linear_counts.ppmi %>%
    select(-study) %>%
    pivot_wider(id_cols = "gene_id", names_from = id, values_from = norm_count),
  icicleCountMatrix = norm_linear_counts.icicle %>%
    select(-study) %>%
    pivot_wider(id_cols = "gene_id", names_from = id, values_from = norm_count),
  ppmiMeta = meta.ppmi, icicleMeta = meta.icicle
) +
  ggtitle("Gene expression")
plot.gene_counts_condition

# Merge into one figure for supp
plot.fsj_counts_condition + plot.bsjHosts_counts_condition + plot.gene_counts_condition + 
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
ggsave(here("combinedLinCirc/output/figures/supp/otherRNAtypeNormCounts.svg"),
       height = 4, width = 10)


# relationship between BSJ + FSJ fold changes -------------------------------------------------
BSJ_vs_FSJ <- bind_rows(
  full_join(bsj_results.ppmi, fsj_results.ppmi,
            by = c("coord_id", "study"),
            suffix = c(".bsj", ".fsj")
  ),
  full_join(bsj_results.icicle, fsj_results.icicle,
            by = c("coord_id", "study"),
            suffix = c(".bsj", ".fsj")
  )
)

# Calculate correlation
bsj_fsj_fc_cor <- map_dfr(cohorts, function(cohort) {
  BSJ_vs_FSJ %>%
    filter(study == {{ cohort }}) %>%
    lm(log2FoldChange.bsj ~ log2FoldChange.fsj, data = .) %>%
    glance() %>%
    mutate(
      study = {{ cohort }},
      r.squared = round(r.squared, digits = 2),
      p.value = if_else(p.value < 0.001,
                        expSup(p.value),
                        format.pval(p.value, digits = 1)
      )
    ) %>%
    select(r.squared, p.value, study)
})
bsj_fsj_fc_cor


# How many highly expressed BSJs in each cohort had no corresponding FSJs detected?
BSJ_vs_FSJ %>%
  filter(is.na(log2FoldChange.fsj)) %>%
  group_by(study) %>%
  count()

# plot
plot.junc_fc <- BSJ_vs_FSJ %>%
  ggplot(aes(x = log2FoldChange.bsj, y = log2FoldChange.fsj)) +
  geom_point(aes(colour = study), alpha = 0.7) +
  geom_xsidedensity(aes(y = stat(density)), fill = rna_colours["BSJ"], colour = NA) +
  geom_ysidedensity(aes(x = stat(density)), fill = rna_colours["FSJ"], colour = NA) +
  theme_ggside_void() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_richtext(
    data = bsj_fsj_fc_cor, aes(x = -Inf, y = Inf, label = paste0(
      " <i>R<sup>2</sup></i> = ", r.squared, "<br>",
      " <i>P</i> = ", p.value
    )),
    vjust = 1, hjust = "left", size = 3, fill = NA, label.color = NA
  ) +
  scale_colour_manual(values = study_colours) +
  labs(
    x = "log<sub>2</sub>(BSJ Fold Change)",
    y = "log<sub>2</sub>(FSJ Fold Change)"
  ) +
  facet_wrap(~ factor(study, levels = c("PPMI", "ICICLE-PD"))) +
  theme(legend.position = "none")
plot.junc_fc
ggsave(here("combinedLinCirc/output/figures/individual/bsj_vs_fsj_fc.png"),
       height = 5, width = 8, dpi = 600, device = agg_png
)

# compare fold changes between BSJ and FSJs
plot.bsj_fsj_fc <- BSJ_vs_FSJ %>%
  pivot_longer(cols = c(log2FoldChange.fsj, log2FoldChange.bsj), names_to = "type", values_to = "log2FoldChange") %>%
  mutate(
    type = gsub(type, pattern = "log2FoldChange.", replacement = ""),
    type = case_match(type, "bsj" ~ "BSJ", "fsj" ~ "FSJ")
  ) %>%
  ggplot(aes(x = type, y = log2FoldChange)) +
  geom_violin(aes(fill = type), colour = NA) +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  # add on wilcoxon test pvals
  geom_richtext(
    # iterate over cohort differences
    data = map(cohorts, function(study) {
      BSJ_vs_FSJ %>%
        select(coord_id, study, log2FoldChange.bsj, log2FoldChange.fsj) %>%
        filter(study == {{ study }}) %>%
        # convert to tidy format
        pivot_longer(cols = c(log2FoldChange.fsj, log2FoldChange.bsj), names_to = "type", values_to = "log2FoldChange") %>%
        # perform test
        wilcox.test(.$log2FoldChange ~ .$type, data = .) %>%
        tidy() %>%
        # add on study and format pvals for plotting
        mutate(
          study = {{ study }},
          p.value = if_else(p.value < 0.001,
                            expSup(p.value),
                            format.pval(p.value, digits = 1)
          )
        )
    }) %>%
      # combine into one df
      list_rbind(),
    aes(x = 1.5, y = Inf, label = paste0("<i>P</i> = ", p.value)), vjust = 1, size = 3, fill = NA, label.color = NA
  ) +
  labs(x = "Junction type", y = "log<sub>2</sub>(Fold Change)") +
  scale_fill_manual(values = rna_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ factor(study, levels = c("PPMI", "ICICLE-PD"))) +
  theme(legend.position = "none")
plot.bsj_fsj_fc




# circrna regulators --------------------------------------------------------------------------
# import regulators
circ_reg_lit <- read_csv(here("circRNA/data/circ_regulators_lit_evidence.csv"))
length(unique(circ_reg_lit$gene_id))

# DE in PPMI
bsj_reg.ppmi <- circ_reg_lit[, c("gene_id", "evidence")] %>% 
  left_join(linear_results.ppmi, by = "gene_id") %>% 
  mutate(fdr.ppmi = p.adjust(pvalue, method= 'fdr'))

# replicate in ICICLE-PD
bsj_reg.icicle <- bsj_reg.ppmi %>% 
  filter(fdr.ppmi < 0.05,
         log2FoldChange > 0.1 | log2FoldChange < -0.1) %>% 
  left_join(linear_results.icicle, by = "gene_id", suffix = c(".ppmi", ".icicle")) %>% 
  mutate(fdr.icicle = p.adjust(pvalue.icicle, method = "fdr"))

table.bsj_reg <- circ_reg_lit[, c("gene_id", "evidence")] %>% 
  left_join(linear_results.merged, by = "gene_id") %>% 
  left_join(bsj_reg.ppmi[, c("gene_id", "fdr.ppmi")], by = "gene_id") %>% 
  left_join(bsj_reg.icicle[, c("gene_id", "fdr.icicle")], by = "gene_id") %>% 
  select(-contains(c("study", "entrez", "type", "gene_biotype", "padj"))) %>% 
  relocate(gene_symbol, .before = gene_id) %>% 
  arrange(pvalue.ppmi)
# export as supp
write_csv(table.bsj_reg, here("combinedLinCirc/output/tables/circ_regulator_DE.csv"))



# Function to plot counts easily
plotNormCounts <- function(gene_id) {
  de_stats <- linear_results.bind %>%
    filter(gene_id == {{ gene_id }}) %>%
    mutate(
      pvalue = if_else(pvalue < 0.001,
                       expSup(pvalue),
                       format.pval(pvalue, digits = 1)
      ),
      log2FoldChange = round(log2FoldChange, 2)
    )
  bind_rows(norm_linear_counts.ppmi, norm_linear_counts.icicle) %>%
    filter(gene_id == {{ gene_id }}) %>%
    left_join(meta.bound[, c("id", "condition", "study")], by = c("id", "study")) %>%
    mutate(
      study = factor(study, levels = c("PPMI", "ICICLE-PD")),
      condition = factor(condition, levels = c("PD", "Control"))
    ) %>%
    ggplot(aes(x = condition, y = norm_count)) +
    geom_violin(aes(fill = condition), colour = NA) +
    geom_boxplot(outlier.shape = NA, width = 0.2) +
    # geom_jitter(aes(colour = condition), width = 0.1, alpha = 0.5) +
    # stat_summary(fun = "mean", size= 0.5, geom = "crossbar", width = 0.2) +
    scale_fill_manual(values = condition_colours) +
    scale_colour_manual(values = condition_colours) +
    scale_y_log10(expand = expansion(mult = c(0.1, 0.2))) +
    labs(
      x = "Study group",
      y = "Normalised counts",
      fill = "Study group"
    ) +
    facet_wrap(~study) +
    geom_richtext(
      data = de_stats, aes(x = 1.5, y = Inf, label = paste0(
        "log<sub>2</sub>FC = ", log2FoldChange, "<br>",
        "<i>P</i> = ", pvalue
      )),
      vjust = 1, size = 3, fill = NA, label.color = NA
    )
}


# RNASEL is significant - plot counts
plot.rnasel <- plotNormCounts("ENSG00000135828") +
  ggtitle("<i>RNASEL</i>")
plot.rnasel + theme(legend.position = "none")
ggsave(here("combinedLinCirc/output/figures/individual/RNASEL_DE.png"),
       height = 5, width = 4, dpi = 600, device = agg_png
)

# Is PKR differentially expressed? (ENSG00000055332)
plot.pkr <- plotNormCounts("ENSG00000055332") +
  ggtitle("<i>EIF2AK2</i> (PKR)")
plot.pkr


# figure panel --------------------------------------------------------------------------------
# (((plot.imbalance_volcanos / plot.bsj_counts_condition) + plot_layout(heights = c(3, 1))) | 
#     (plot.junc_fc / plot.bsj_fsj_fc / (plot.rnasel | plot.pkr)) + plot_layout(heights = c(2, 1, 2))) +
#   plot_layout(guides = "collect")
(plot.imbalance_forest / (plot.bsj_counts_condition | plot.junc_fc) /
(plot.bsj_fsj_fc | plot.rnasel | plot.pkr)) + plot_layout(heights = c(0.75, 1, 1), guides = "collect") + plot_annotation(tag_levels = "a")
ggsave(here("combinedLinCirc/output/figures/panels/imbalance.svg"),
       height = 10, width = 10)
