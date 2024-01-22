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

expSup <- function(w, digits=1) {
  sprintf(paste0("%.", digits, "fx10<sup>%d</sup>"), w/10^floor(log10(abs(w))), floor(log10(abs(w))))
}

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

## junc info
junc_info <- read_fst(here("circRNA/data/bound_juncInfo.fst"))
# add vst counts to junc info
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


# auc individual junction ---------------------------------------------------------------------

# Function to calculate indiviudal AUCs
coordAUC <- function(coord_id, column) {
  map_dfr(c("PPMI", "ICICLE-PD"), function(study) {
    cpm <- junc_info %>%
      filter(study == {{ study }} & coord_id == {{ coord_id }}) %>%
      select(id, coord_id, feature = {{ column }}, condition)
    perf <- roc(cpm$condition ~ cpm$feature, levels = c("Control", "PD"), ci = TRUE, quiet = TRUE)
    auc_ci <- as.numeric(perf$ci)
    data.frame(
      coord_id = coord_id,
      lower = auc_ci[1],
      auc = auc_ci[2],
      upper = auc_ci[3],
      study = study
    )
  })
}

# tpmAUC <- function(gene_id) {
#   map_dfr(c("PPMI", "ICICLE-PD"), function(study) {
#     tpm <- tpm_counts %>%
#       left_join(meta.bound[, c("id", "study", "condition")], by = c("id", "study")) %>%
#       filter(study == {{ study }} & gene_id == {{ gene_id }}) %>%
#       select(id, gene_id, TPM, condition)
#     perf <- roc(tpm$condition ~ tpm$TPM, levels = c("Control", "PD"), ci = FALSE, quiet = TRUE)
#     data.frame(
#       gene_id = gene_id,
#       auc = as.numeric(auc(perf)),
#       study = study
#     )
#   })
# }


# Run on all gene counts (TPMs) (doesn't work - predicted to take 12 hours)
# system.time(
# tpm_auc <- map(.x = linear_results.merged %>%
#                  filter(pvalue.ppmi < 0.05) %>% drop_na() %>% pull(gene_id),
#     .f = tpmAUC,
#     .progress = TRUE) %>%
#   list_rbind()
# )

# Run on all BSJs in both cohorts
#plan(multisession, workers = 2)
bsj_auc <- map2_dfr(
  .x = bsj_results.merged %>% drop_na() %>% pull(coord_id),
  .y = "bsj_vst",
  .f = coordAUC,
  .progress = TRUE
  #.options = furrr_options(seed = TRUE)
) %>% mutate(type = "BSJ")
#plan(sequential)

# Run on all FSJs in both cohorts
#plan(multisession, workers = 2)
# still use bsj results df as fsj coords are the same
fsj_auc <- map2_dfr(
  .x = bsj_results.merged %>% drop_na() %>% pull(coord_id),
  .y = "fsj_vst",
  .f = coordAUC,
  .progress = TRUE
  #.options = furrr_options(seed = TRUE)
) %>% mutate(type = "FSJ")
#plan(sequential)

# Run all ratios in both cohorts
#plan(multisession, workers = 2)
# still use bsj results df as fsj coords are the same
ratio_auc <- map2_dfr(
  .x = bsj_results.merged %>% drop_na() %>% pull(coord_id),
  .y = "junc_ratio",
  .f = coordAUC,
  .progress = TRUE
  #.options = furrr_options(seed = TRUE)
) %>% mutate(type = "ratio")
#plan(sequential)

# Combine all individual junction AUCs into one df
all_junc_auc <- bind_rows(bsj_auc, fsj_auc, ratio_auc)
all_junc_auc %>% 
  filter(coord_id %in% c("3:63912587-63913225:+", "3:182884752-182887713:+"))
all_junc_auc %>% 
  group_by(study) %>% 
  summarise(med_auc = median(auc))
all_junc_auc %>% 
  wilcox.test(auc ~ study, data = .) %>% 
  tidy()
all_junc_auc %>% 
  filter(auc > 0.6) %>% 
  group_by(type, study) %>% tally()


# Plot AUCs (BSJ, FSJ, BSJ:FSJ ratio) in each cohort
df.junc_auc <- bind_rows(bsj_auc, fsj_auc, ratio_auc) %>%
  pivot_wider(id_cols = c(coord_id, type), names_from = study, values_from = auc) %>%
  mutate(
    type = recode(type, "ratio" = "BSJ:FSJ ratio"),
    type = factor(type, levels = c("BSJ", "FSJ", "BSJ:FSJ ratio"))
  )
df.junc_auc %>% 
  filter(PPMI > 0.6 & `ICICLE-PD` > 0.6)
junc_auc_cor <- map(c("BSJ", "FSJ", "BSJ:FSJ ratio"), function(type) {
  cor.test(~ PPMI + `ICICLE-PD`, data = filter(df.junc_auc, type == {{ type }}), method = "spearman") %>%
    tidy() %>%
    mutate(type = {{ type }})
}) %>%
  list_rbind() %>%
  mutate(
    rho = round(estimate, digits = 2),
    p.value = if_else(p.value < 0.001,
      expSup(p.value, 1),
      format.pval(p.value, digits = 2)
    )
  )
plot.junc_auc <- df.junc_auc %>%
  ggplot(aes(x = PPMI, y = `ICICLE-PD`)) +
  geom_point(aes(colour = type), alpha = 0.4) +
  scale_colour_manual(values = rna_colours) +
  xlim(0.45, 0.65) +
  ylim(0.45, 0.65) +
  geom_xsidedensity(aes(y = stat(density)), fill = study_colours["PPMI"], colour = NA, alpha = 0.7) +
  geom_ysidedensity(aes(x = stat(density)), fill = study_colours["ICICLE-PD"], colour = NA, alpha = 0.7) +
  theme_ggside_void() +
  geom_vhlines(xintercept = 0.5, yintercept = 0.5, linetype = "dashed", alpha = 0.2) +
  geom_richtext(
    data = junc_auc_cor, aes(x = Inf, y = -Inf, label = paste0(
      " rho = ", rho, "<br>",
      " <i>P</i> = ", p.value
    )),
    vjust = 0, hjust = "right", size = 3, fill = NA, label.color = NA 
  ) +
  labs(x = "PPMI AUC", y = "ICICLE-PD AUC", colour = "Feature type") +
  facet_wrap(~ factor(type, levels = c("BSJ", "FSJ", "BSJ:FSJ ratio")), nrow = 3)
plot.junc_auc
ggsave(here("combinedLinCirc/output/figures/individual/individual_aucs.png"),
  device = agg_png,
  height = 3, width = 9, dpi = 600
)


# regularised logistic regression -------------------------------------------------------------
# Import models
models <- read_rds(here("combinedLinCirc/output/models.rds"))

# function to tidy ROC objects
tidyROC <- function(rocObject) {
  data.frame(
    threshold = rocObject[["thresholds"]],
    specificity = rocObject[["specificities"]],
    sensitivity = rocObject[["sensitivities"]]
  )
}

# function to extract performance
modelPerformance <- function(modelOutput, coefs) {
  ##### ROC AUC
  # get tidy roc object for the training set
  train_roc <- modelOutput[["train"]][["roc"]]
  train_auc <- ci.auc(train_roc)
  train_auc <- data.frame(ci_lower = train_auc[1], auc = train_auc[2], ci_higher = train_auc[3])
  train_auc <- train_auc %>% mutate(study = 'PPMI')
  train_tidy <- tidyROC(train_roc)
  train_tidy <- train_tidy %>% mutate(study = "PPMI")
  # get tidy roc object for test set
  test_df <- modelOutput[["test"]]
  test_roc <- roc(test_df$condition ~ test_df$s1, levels = c("Control", "PD"))
  test_auc <- ci.auc(test_roc)
  test_auc <- data.frame(ci_lower = test_auc[1], auc = test_auc[2], ci_higher = test_auc[3])
  test_auc <- test_auc %>% mutate(study = 'ICICLE-PD')
  test_tidy <- tidyROC(test_roc)
  test_tidy <- test_tidy %>% mutate(study = "ICICLE-PD")
  # merge tidy train + test dfs for plotting ROC curves
  roc_curve_df <- bind_rows(train_tidy, test_tidy) %>%
    mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")))
  # merge auc data frames for plotting AUCs with 95% CIs
  auc_df <- bind_rows(train_auc, test_auc) %>% 
    mutate(study = factor(study, levels = c('PPMI', 'ICICLE-PD')))
  # get coefficients from training data
  if (coefs == TRUE) {
    coef_df <- modelOutput[["train"]][["final_coef"]] %>%
      rownames_to_column("feature") %>%
      # remove intercept
      filter(feature != "(Intercept)") %>%
      arrange(desc(coef))
  } else if (coefs == FALSE) {
    coef_df <- data.frame(coef = "NONE")
  } else {
    stop("Set coefs to TRUE/FALSE depending on whether a logistic regression model has been made")
  }
  ##### Export plots as list
  output <- list(coefs = coef_df, rocs = roc_curve_df, aucs = auc_df)
  return(output)
}

### Linear RNA model performance
linear_model <- modelPerformance(modelOutput = models[["lin_model"]], coefs = TRUE)
linear_model$aucs %>% 
  mutate(across(where(is.numeric), \(x) round(x, 2)))

### BSJ model performance
bsj_model <- modelPerformance(modelOutput = models[["bsj_model"]], coefs = TRUE)
bsj_model$aucs %>% 
  mutate(across(where(is.numeric), \(x) round(x, 2)))

### FSJ model performance
fsj_model <- modelPerformance(modelOutput = models[["fsj_model"]], coefs = TRUE)
fsj_model$aucs %>% 
  mutate(across(where(is.numeric), \(x) round(x, 2)))

### BSJ:FSJ ratio model
ratio_model <- modelPerformance(modelOutput = models[["ratio_model"]], coefs = TRUE)
ratio_model$aucs %>% 
  mutate(across(where(is.numeric), \(x) round(x, 2)))

### Combined BSJ + linear RNA model
combined_lin_bsj_model <- modelPerformance(modelOutput = models[["combined_lin_bsj_model"]], coefs = TRUE)
combined_lin_bsj_model$aucs %>% 
  mutate(across(where(is.numeric), \(x) round(x, 2)))

### Imbalance model
imbalance_model <- models[["imbalance_model"]]
# do ROC analysis separate as it's not a logistic regression model
# roc objects
imbalance_roc.ppmi <- roc(condition ~ total_score, data = imbalance_model[["train"]])
imbalance_roc.icicle <- roc(condition ~ total_score, data = imbalance_model[["test"]])
# add AUCs from roc objects to modellist object
imbalance_model[["aucs"]] <- bind_rows(data.frame(ci_lower = ci.auc(imbalance_roc.ppmi)[1],
                                                  auc = ci.auc(imbalance_roc.ppmi)[2],
                                                  ci_higher = ci.auc(imbalance_roc.ppmi)[3],
                                                  study = "PPMI"),
                                       data.frame(ci_lower = ci.auc(imbalance_roc.icicle)[1],
                                                  auc = ci.auc(imbalance_roc.icicle)[2],
                                                  ci_higher = ci.auc(imbalance_roc.icicle)[3],
                                                  study = "ICICLE-PD"))
# add tidied roc objects to model list object
imbalance_model[["rocs"]] <- bind_rows(
  imbalance_roc.ppmi %>% tidyROC() %>% mutate(study = "PPMI"),
  imbalance_roc.icicle %>% tidyROC() %>% mutate(study = "ICICLE-PD")
) %>% 
  data.frame()
imbalance_model$aucs %>% 
  mutate(across(where(is.numeric), \(x) round(x, 2)))



## Number of predictors in each model
no_predictors_df <- data.frame(
  totalRNA = nrow(linear_model[["coefs"]]),
  BSJ = nrow(bsj_model[["coefs"]]),
  FSJ = nrow(fsj_model[["coefs"]]),
  ratios = nrow(ratio_model[["coefs"]]),
  combinedLinBSJ = nrow(combined_lin_bsj_model[["coefs"]])
) %>%
  pivot_longer(cols = everything(), names_to = "Predictor", values_to = "n") %>%
  mutate(Predictor = recode(Predictor,
                            totalRNA = "Gene",
                            ratios = "BSJ:FSJ ratio",
                            combinedLinBSJ = "Gene and BSJ"
  ))
no_predictors_df
no_predictors_plot <- no_predictors_df %>%
  ggplot(aes(x = Predictor, y = n)) +
  geom_col(aes(fill = Predictor)) +
  geom_text(aes(label = paste0("n = ", n)), vjust = -1, size = 3) +
  scale_fill_manual(values = rna_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Number of predictors in model") +
  theme(legend.position = "none")
no_predictors_plot



# roc curves ----------------------------------------------------------------------------------
# Plot ROC curves
# model_stats <- model_table %>%
#   pivot_wider(names_from = model, values_from = AUC)
roc_curves <- bind_rows(
  linear_model[["rocs"]] %>% mutate(type = "Gene"),
  bsj_model[["rocs"]] %>% mutate(type = "BSJ"),
  fsj_model[["rocs"]] %>% mutate(type = "FSJ"),
  ratio_model[["rocs"]] %>% mutate(type = "BSJ:FSJ ratio"),
  combined_lin_bsj_model[["rocs"]] %>% mutate(type = "Gene and BSJ"),
  imbalance_model[["rocs"]] %>% mutate(type = "Imbalance")
) %>%
  mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD"))) %>%
  ggplot(aes(x = specificity, y = sensitivity)) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed") +
  geom_path(aes(colour = type), size = 0.7) +
  # geom_text(
  #   data = model_stats, aes(x = -Inf, y = -Inf, label = paste0(
  #     "AUC \n",
  #     "Gene: ", formatC(`Linear RNA`, 2), "\n",
  #     "BSJ: ", formatC(BSJ, 2), "\n",
  #     "FSJ: ", formatC(FSJ, 2), "\n",
  #     "BSJ:FSJ ratio: ", formatC(ratio, 2), "\n",
  #     "Gene and BSJ: ", formatC(`Gene and BSJ`, 2), "\n",
  #     "Imbalance: ", formatC(Imbalance, 2), "\n"
  #   )),
  #   vjust = 0, hjust = 1, size = 3
# ) +
scale_y_continuous(expand = c(0, 0)) +
  scale_x_reverse() +
  scale_colour_manual(values = rna_colours) +
  labs(
    x = "Specificity",
    y = "Sensitivity",
    colour = "Predictor"
  ) +
  facet_wrap(~ factor(study, levels = c("PPMI", "ICICLE-PD")), scales = "free") +
  theme(legend.position = "none")
roc_curves
ggsave(here("combinedLinCirc/output/figures/individual/plr_ROCcurves.png"),
       height = 5, width = 8, dpi = 600, device = agg_png
)


# auc plot ------------------------------------------------------------------------------------
plot.model_aucs <- bind_rows(
  linear_model$aucs %>% mutate(type = "Gene"),
  bsj_model$aucs %>% mutate(type = "BSJ"),
  fsj_model$aucs %>% mutate(type = "FSJ"),
  ratio_model$aucs %>% mutate(type = "BSJ:FSJ ratio"),
  combined_lin_bsj_model$aucs %>% mutate(type = "Gene and BSJ"),
  imbalance_model$aucs %>% mutate(type = "Imbalance")
) %>% 
  mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD"))) %>% 
  ggplot(aes(y = type, x = auc, colour = type)) +
  facet_wrap(~ study, scales = "free_y") +
  geom_point(show.legend = FALSE) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_higher)) +
  scale_colour_manual(values = rna_colours) +
  labs(x = "AUC", y = "Predictor", colour = "Predictor") +
  theme(legend.position = 'bottom')
plot.model_aucs


# figure panel --------------------------------------------------------------------------------
right <- roc_curves / plot.model_aucs + plot_layout(heights = c(1, 1))
left <- plot.junc_auc + theme(legend.position = "none")
left + right + plot_layout(widths = c(1, 2))
ggsave(here('combinedLinCirc/output/figures/panels/classification.svg'),
       height = 9, width = 10)
