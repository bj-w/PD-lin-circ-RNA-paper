library(here)
library(data.table)
library(fst)
library(tidyverse)
library(glmnet)
library(nestedcv)
library(caret)
library(pROC)

set.seed(1997, "L'Ecuyer-CMRG")

# model function ----------------------------------------------------------
formatFeatureMatrix <- function(trainCounts, testCounts, trainMeta, testMeta) {
  # genes not present in both cohorts
  absent_lin <- setdiff(trainCounts$feature, testCounts$feature)
  # remove these genes from counts
  trainCounts <- trainCounts %>% filter(!feature %in% absent_lin)
  testCounts <- testCounts %>% filter(!feature %in% absent_lin)
  # make sure same features in both
  trainCounts <- trainCounts %>% filter(feature %in% testCounts$feature)
  testCounts <- testCounts %>% filter(feature %in% trainCounts$feature)
  # format count matrices
  trainCounts <- trainCounts %>% column_to_rownames("feature")
  testCounts <- testCounts %>% column_to_rownames("feature")
  # transpose
  trainCounts <- trainCounts %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    as_tibble()
  testCounts <- testCounts %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    as_tibble()
  # add condition onto counts
  trainCounts <- trainCounts %>%
    left_join(trainMeta[, c("id", "condition")], by = "id") %>%
    relocate(condition, .after = id) %>%
    column_to_rownames("id")
  testCounts <- testCounts %>%
    left_join(testMeta[, c("id", "condition")], by = "id") %>%
    relocate(condition, .after = id) %>%
    column_to_rownames("id")
  # export as list object
  output <- list(train = trainCounts, test = testCounts)
  return(output)
}

trainLR <- function(counts, filter) {
  set.seed(1997, "L'Ecuyer-CMRG")
  counts$condition <- as.factor(ifelse(counts$condition == "PD", 1, 0))
  x <- as.matrix(counts[, -1])
  y <- counts$condition
  if (filter == TRUE) {
    nested_results <- nestcv.glmnet(
      x = x, y = y,
      filterFUN = wilcoxon_filter,
      filter_options = list(p_cutoff = 0.05, rsq_cutoff = 0.9),
      n_outer_folds = 10, n_inner_folds = 10, outer_train_predict = TRUE,
      family = "binomial", alphaSet = seq(0, 1, 0.05),
      finalCV = TRUE, cv.cores = 2
    )
  } else {
    nested_results <- nestcv.glmnet(
      x = x, y = y,
      n_outer_folds = 10, n_inner_folds = 10, outer_train_predict = TRUE,
      family = "binomial", alphaSet = seq(0, 1, 0.05),
      finalCV = TRUE, cv.cores = 2
    )
  }
  return(nested_results)
}

modelFunction <- function(listModelCounts) {
  # train model on PPMI
  train_model <- trainLR(counts = listModelCounts[["train"]], filter = FALSE)
  print(train_model[["roc"]])
  # predict on ICICLE-PD
  test_predict <- predict(train_model,
    newdata = as.matrix(listModelCounts[["test"]][, -1]),
    type = "response"
  )
  test_predict <- test_predict %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    left_join(meta.icicle[, c("id", "condition")], by = "id")
  print(roc(test_predict$condition ~ test_predict$s1, levels = c("Control", "PD"), direction = "<"))
  # export
  return(list(
    train = train_model,
    test = test_predict
  ))
}

# IMPORT GENERAL DATA -----------------------------------------------------
# metadata
meta.ppmi <- read_rds("circRNA/data/ppmi_metadata.rds") %>%
  mutate(id = as.character(id))
meta.icicle <- read_rds("circRNA/data/icicle_metadata.rds")

# junction counts
jc <- read_fst("circRNA/data/bound_juncInfo.fst") %>%
  rename(feature = coord_id)

# CLASSIFICATION WITH TOTAL RNA -------------------------------------------

# import significant features from total rna results
sig_linear_features <- read_csv("linear/output/ppmi_deseqResults.csv") %>%
  filter(pvalue < 0.05) %>%
  pull(ensembl)

# import counts
counts_lin.ppmi <- fread("linear/data/ppmi_deseqFilteredNormalisedCounts.csv", header = TRUE) %>%
  rename(feature = ensembl) %>%
  filter(feature %in% sig_linear_features)
counts_lin.icicle <- fread("linear/data/icicle_deseqFilteredNormalisedCounts.csv", header = TRUE) %>%
  rename(feature = ensembl) %>%
  filter(feature %in% sig_linear_features)

# getTPM <- function(){
#   formatTPM <- function(txi, study){
#     txi[["abundance"]] %>%
#       as.data.frame() %>%
#       rownames_to_column('gene_id') %>%
#       pivot_longer(cols = !gene_id, names_to = 'path', values_to = 'TPM') %>%
#       mutate(study = {{ study }})
#   }
#   cat("Importing PPMI TPMs...\n")
#   linear_meta.ppmi <- read_csv(here('ppmi/output/BL_iPD_HC_metadata.csv')) %>%
#     mutate(id = as.character(id))
#   ppmi <- readRDS(here("linear/data/ppmi_txiForDeseq.rds")) %>%
#     formatTPM("PPMI") %>%
#     # add on sample ID that corresponds to fastq full name
#     left_join(linear_meta.ppmi[, c('id', 'path')], by = 'path') %>%
#     # remove fastq name path by selecting only the specific columns we need
#     select(gene_id, id, TPM, study)
#   cat("Importing ICICLE-PD TPMs...\n")
#   icicle <- readRDS(here("linear/data/icicle_txiForDeseq.rds")) %>%
#     formatTPM("ICICLE-PD") %>%
#     rename('id' = 'path')
#   tpm <- bind_rows(ppmi, icicle)
#   return(tpm)
# }
# # import TPM counts
# tpm_counts <- getTPM()
# # filter for PPMI and log2 + 1
# counts_lin.ppmi <- tpm_counts %>%
#   mutate(TPM = log2(TPM + 1)) %>%
#   filter(study == 'PPMI') %>%
#   pivot_wider(id_cols = 'gene_id', names_from = 'id', values_from = 'TPM', values_fill = 0) %>%
#   rename(feature = gene_id)
# # filter for ICICLE and log2 +1
# counts_lin.icicle <- tpm_counts %>%
#   mutate(TPM = log2(TPM + 1)) %>%
#   filter(study == 'ICICLE-PD') %>%
#   pivot_wider(id_cols = 'gene_id', names_from = 'id', values_from = 'TPM', values_fill = 0) %>%
#   rename(feature = gene_id)
# # clean up
# rm(tpm_counts)

# only use DEGs from PPMI
counts_lin.ppmi <- counts_lin.ppmi %>% filter(feature %in% sig_linear_features)
counts_lin.icicle <- counts_lin.icicle %>% filter(feature %in% sig_linear_features)

lin_model_counts <- formatFeatureMatrix(
  trainCounts = counts_lin.ppmi,
  testCounts = counts_lin.icicle,
  trainMeta = meta.ppmi,
  testMeta = meta.icicle
)

# model + prediction
lin_model <- modelFunction(lin_model_counts)


# CLASSIFICATION WITH BSJ COUNTS ------------------------------------------
# import sig BSJs
bsj_results <- read_csv("circRNA/output/ppmi_BSJresults.csv")
sig_bsj_features <- bsj_results %>%
  filter(pvalue < 0.05) %>%
  pull(coord_id)

# import counts
bsj_counts.ppmi <- read_csv("circRNA/data/ppmi_vstBSJCounts.csv") %>%
  rename(feature = coord_id) %>%
  filter(feature %in% sig_bsj_features)
bsj_counts.icicle <- read_csv("circRNA/data/icicle_vstBSJCounts.csv") %>%
  rename(feature = coord_id) %>%
  filter(feature %in% sig_bsj_features)

# bsj_counts.ppmi <- jc %>% 
#   filter(study == "PPMI") %>% 
#   select(feature, id, bsj_perMapped) %>% 
#   pivot_wider(id_cols = feature, names_from = id, values_from = bsj_perMapped, values_fill = 0) %>% 
#   mutate(across(where(is.numeric), ~log2(.x + 1)))
# 
# bsj_counts.icicle <- jc %>% 
#   filter(study == "ICICLE-PD") %>% 
#   select(feature, id, bsj_perMapped) %>% 
#   pivot_wider(id_cols = feature, names_from = id, values_from = bsj_perMapped, values_fill = 0) %>% 
#   mutate(across(where(is.numeric), ~log2(.x + 1)))

# format for model building
bsj_model_counts <- formatFeatureMatrix(
  trainCounts = bsj_counts.ppmi, testCounts = bsj_counts.icicle,
  trainMeta = meta.ppmi, testMeta = meta.icicle
)
# model + prediction
bsj_model <- modelFunction(bsj_model_counts)


# CLASSIFICATION WITH FSJ COUNTS ------------------------------------------

# import counts
fsj_counts.ppmi <- read_csv("circRNA/data/ppmi_vstFSJCounts.csv") %>%
  rename(feature = coord_id) %>%
  filter(feature %in% sig_bsj_features)
fsj_counts.icicle <- read_csv("circRNA/data/icicle_vstFSJCounts.csv") %>%
  rename(feature = coord_id) %>%
  filter(feature %in% sig_bsj_features)
# format for model building
fsj_model_counts <- formatFeatureMatrix(
  trainCounts = fsj_counts.ppmi, testCounts = fsj_counts.icicle,
  trainMeta = meta.ppmi, testMeta = meta.icicle
)
# model + prediction
fsj_model <- modelFunction(fsj_model_counts)


# BSJ:FSJ RATIO MODEL -----------------------------------------------------
ratio_counts.ppmi <- jc %>%
  filter(
    study == "PPMI",
    feature %in% sig_bsj_features
  ) %>%
  select(id, feature, junc_ratio) %>%
  pivot_wider(id_cols = feature, names_from = id, values_from = junc_ratio, values_fill = 0)
ratio_counts.icicle <- jc %>%
  filter(
    study == "ICICLE-PD",
    feature %in% sig_bsj_features
  ) %>%
  select(id, feature, junc_ratio) %>%
  pivot_wider(id_cols = feature, names_from = id, values_from = junc_ratio, values_fill = 0)
# format
ratio_model_counts <- formatFeatureMatrix(
  trainCounts = ratio_counts.ppmi, testCounts = ratio_counts.icicle,
  trainMeta = meta.ppmi, testMeta = meta.icicle
)
# model + prediction
ratio_model <- modelFunction(ratio_model_counts)


# COMBINED LINEAR + BSJ ---------------------------------------------------
# create a matrix from linear + BSJ counts (just bind the cols together)
combined_lin_bsj_model_counts <- list(
  train = bind_cols(lin_model_counts$train, bsj_model_counts$train[, -1]),
  test = bind_cols(lin_model_counts$test, bsj_model_counts$test[, -1])
)

# model + prediction
combined_lin_bsj_model <- modelFunction(combined_lin_bsj_model_counts)


# CIRCRNA IMBALANCE MODEL ---------------------------------------------------------------------

control_vst <- bsj_counts.ppmi %>%
  pivot_longer(cols = !feature, names_to = "id", values_to = "vst") %>%
  left_join(meta.ppmi[, c("id", "condition")], by = "id") %>%
  filter(condition == "Control") %>%
  group_by(feature) %>%
  summarise(control_vst = mean(vst))

fold_changes <- bsj_results %>%
  # encode whether bsj is up or down in PD
  mutate(direction = case_when(
    log2FoldChange < 0 ~ "down",
    log2FoldChange > 0 ~ "up",
    TRUE ~ "noDiff"
  )) %>%
  # rename coord id to feature for joining
  rename(feature = coord_id)
# convert counts to tidy format
imbalance_model <- bsj_counts.ppmi %>%
  pivot_longer(cols = !feature, names_to = "id", values_to = "vst")
# add on direction from DEA
imbalance_model <- imbalance_model %>% left_join(fold_changes[, c("feature", "direction")], by = "feature")
# add on control mean vst counts
imbalance_model <- imbalance_model %>% left_join(control_vst, by = "feature")
# add on study group
imbalance_model <- imbalance_model %>% left_join(meta.ppmi[, c("id", "condition")], by = "id")
# create model score - add 1 if expression agrees
imbalance_model <- imbalance_model %>% mutate(score = case_when(
  direction == "down" & vst < control_vst ~ 1,
  direction == "up" & vst > control_vst ~ 1,
  TRUE ~ 0
))
# calculate total score
imbalance_model <- imbalance_model %>%
  group_by(id) %>%
  mutate(total_score = sum(score)) %>%
  select(id, total_score, condition) %>%
  distinct()

# repeat in ICICLE-PD
imbalance_model_replicate <- bsj_counts.icicle %>% pivot_longer(cols = !feature, names_to = "id", values_to = "vst")
# add direction from PPMI
imbalance_model_replicate <- imbalance_model_replicate %>% left_join(fold_changes[, c("feature", "direction")], by = "feature")
# add on control vst from PPMI
imbalance_model_replicate <- imbalance_model_replicate %>% left_join(control_vst, by = "feature")
# calculate score
imbalance_model_replicate <- imbalance_model_replicate %>% mutate(score = case_when(
  direction == "down" & vst < control_vst ~ 1,
  direction == "up" & vst > control_vst ~ 1,
  TRUE ~ 0
))
# add on ICICLE study groups
imbalance_model_replicate <- imbalance_model_replicate %>% left_join(meta.icicle[, c("id", "condition")], by = "id")
# total score
imbalance_model_replicate <- imbalance_model_replicate %>%
  group_by(id) %>%
  mutate(total_score = sum(score)) %>%
  select(id, total_score, condition) %>%
  distinct()
# export model predictions
imbalance_model <- list(train = imbalance_model, test = imbalance_model_replicate)

# EXPORT MODELS -----------------------------------------------------------
models <- list(
  lin_model = lin_model,
  bsj_model = bsj_model,
  fsj_model = fsj_model,
  ratio_model = ratio_model,
  combined_lin_bsj_model = combined_lin_bsj_model,
  imbalance_model = imbalance_model
)
write_rds(models, "combinedLinCirc/output/models.rds")