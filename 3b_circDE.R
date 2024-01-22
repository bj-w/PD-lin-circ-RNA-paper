library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(fst)
library(DOSE)
library(clusterProfiler)
library(tidyverse)

# import common functions
source("combinedLinCirc/PD-RNA/functions.R")


# IMPORT METADATA ---------------------------------------------------------
meta.ppmi <- read_rds("circRNA/data/ppmi_metadata.rds") %>%
  mutate(
    pct_intronic_bases = scale(pct_intronic_bases),
  )

meta.icicle <- read_rds("circRNA/data/icicle_metadata.rds") %>%
  mutate(
    pct_intronic_bases = scale(pct_intronic_bases),
    median_cv_coverage = scale(median_cv_coverage)
  )

# IMPORT UNFILTERED AND DE DDS --------------------------------------------
bsj_dds.ppmi <- read_rds('circRNA/data/ppmi_bsj_DDS_preFilter.rds')
bsj_dds.icicle <- read_rds('circRNA/data/icicle_bsj_DDS_preFilter.rds')
fsj_dds.ppmi <- read_rds('circRNA/data/ppmi_fsj_DDS_preFilter.rds')
fsj_dds.icicle <- read_rds('circRNA/data/icicle_fsj_DDS_preFilter.rds')


# OVERWRITE DESIGN --------------------------------------------------------
# OG DESIGN JUST HAS CONDITION FOR QC AND SETUP
# specify design formula
design.ppmi <- ~ sex + batch + age_bin + pct_intronic_bases + condition
design.icicle <- ~ sex + batch + age_bin + pct_intronic_bases + median_cv_coverage + condition

bsj_dds.ppmi@design <- design.ppmi
bsj_dds.icicle@design <- design.ppmi
fsj_dds.ppmi@design <- design.ppmi
fsj_dds.icicle@design <- design.icicle

# PREFILTERING ------------------------------------------------------------
# filter BSJs based on a filtering strategy
filterDDS <- function(dds, metadata) {
  # apply raw count filtering
  keep <- rowSums(counts(dds, normalized = FALSE) > 10) > min(table(metadata$condition))
  # how many genes removed
  genesRemoved <- function(keep, metadata) {
    total <- nrow(data.frame(keep))
    ntrue <- data.frame(logical = keep) %>%
      filter(logical == TRUE) %>%
      nrow()
    nfalse <- data.frame(logical = keep) %>%
      filter(logical == FALSE) %>%
      nrow()
    print(paste0("Total number of genes tested: ", total))
    print(paste0("Number of genes kept: ", ntrue, " (", round(ntrue / total * 100), "%)"))
    print(paste0("Number of genes removed: ", nfalse, " (", round(nfalse / total * 100), "%)"))
  }
  genesRemoved(keep, meta_ppmi)
  dds <- dds[keep, ]
  return(dds)
}

bsj_dds.ppmi <- filterDDS(dds = bsj_dds.ppmi, metadata = meta.ppmi)
bsj_dds.icicle <- filterDDS(dds = bsj_dds.icicle, metadata = meta.icicle)

# now filter the FSJ counts for the same junctions
fsj_dds.ppmi <- fsj_dds.ppmi[rownames(fsj_dds.ppmi) %in% rownames(bsj_dds.ppmi), ]
fsj_dds.icicle <- fsj_dds.icicle[rownames(fsj_dds.icicle) %in% rownames(bsj_dds.icicle), ]


# DIFFERENTIAL EXPRESSION -------------------------------------------------
getDDS <- function(dds, runOrImport, cohort, type) {
  if (runOrImport == "run") {
    cat("Running DESeq2 (might take a while)\n")
    dds <- DESeq(dds)
    write_rds(dds, paste0("circRNA/data/", cohort, "_dds_", type, ".rds"))
    return(dds)
  } else if (runOrImport == "import") {
    cat("Importing previously run DESeq2 analysis DeseqDataSet\n")
    dds <- read_rds(paste0("circRNA/data/", cohort, "_dds_", type, ".rds"))
    return(dds)
  } else {
    stop("runOrLoad must be either run or import")
  }
}

# DE - RUN OR IMPORT ------------------------------------------------------

# BSJ
bsj_dds.ppmi <- getDDS(dds = bsj_dds.ppmi, runOrImport = "run", cohort = "ppmi", type = "bsj")
bsj_dds.icicle <- getDDS(dds = bsj_dds.icicle, runOrImport = "run", cohort = "icicle", type = "bsj")

fsj_dds.ppmi <- getDDS(dds = fsj_dds.ppmi, runOrImport = "run", cohort = "ppmi", type = "fsj")
fsj_dds.icicle <- getDDS(dds = fsj_dds.icicle, runOrImport = "run", cohort = "icicle", type = "fsj")

# RESULTS -----------------------------------------------------------------

# import junc file
junc_info <- read_fst("circRNA/data/bound_juncInfo.fst")

# import junc info file to get annotations
coord_anno <- junc_info %>%
  select(coord_id, gene_id, gene_name) %>%
  distinct()

# set contrast level
contrast <- c("condition", "PD", "Control")

getDeseqResults <- function(dds, rowname) {
  results <- results(dds, contrast = contrast, independentFiltering = FALSE, alpha = 0.05)
  print(summary(results))
  results <- results %>%
    as.data.frame() %>%
    rownames_to_column(rowname)
  return(results)
}

# BSJ
bsj_results.ppmi <- getDeseqResults(dds = bsj_dds.ppmi, rowname = "coord_id") %>%
  left_join(coord_anno, by = "coord_id")
write_csv(bsj_results.ppmi, "circRNA/output/ppmi_BSJresults.csv")

bsj_results.icicle <- getDeseqResults(dds = bsj_dds.icicle, rowname = "coord_id") %>%
  left_join(coord_anno, by = "coord_id")
write_csv(bsj_results.icicle, "circRNA/output/icicle_BSJresults.csv")

fsj_results.ppmi <- getDeseqResults(dds = fsj_dds.ppmi, rowname = "coord_id") %>%
  left_join(coord_anno, by = "coord_id")
write_csv(fsj_results.ppmi, "circRNA/output/ppmi_FSJresults.csv")

fsj_results.icicle <- getDeseqResults(dds = fsj_dds.icicle, rowname = "coord_id") %>%
  left_join(coord_anno, by = "coord_id")
write_csv(fsj_results.icicle, "circRNA/output/icicle_FSJresults.csv")


exportFilteredCounts <- function(dds, cohort, type) {
  # raw counts
  counts(dds, normalized = FALSE) %>%
    as.data.frame() %>%
    rownames_to_column("coord_id") %>%
    write_csv(paste0("circRNA/data/", cohort, "_deseq", type, "FilteredRawCounts.csv"))
  cat(paste0(cohort, ": exported filtered raw counts"), "\n")
  # normalised counts
  counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("coord_id") %>%
    write_csv(paste0("circRNA/data/", cohort, "_deseq", type, "FilteredNormalisedCounts.csv"))
  cat(paste0(cohort, ": exported filtered normalised counts"), "\n")
  # VST counts
  assay(varianceStabilizingTransformation(dds, blind = FALSE)) %>%
    as.data.frame() %>%
    rownames_to_column("coord_id") %>%
    write_csv(paste0("circRNA/data/", cohort, "_vst", type, "Counts.csv"))
  cat(paste0(cohort, ": exported VST counts"), "\n")
}

exportFilteredCounts(bsj_dds.ppmi, "ppmi", "BSJ")
exportFilteredCounts(bsj_dds.icicle, "icicle", "BSJ")

exportFilteredCounts(fsj_dds.ppmi, "ppmi", "FSJ")
exportFilteredCounts(fsj_dds.icicle, "icicle", "FSJ")



# # TEST (don't run) ----------------------------------------------------------------------------
# # BSJ limited covariates as test
# testDDS <- function(dds) {
#   test_dds <- dds
#   test_dds@design <- design <- ~ sex + batch + age_bin + condition
#   test_dds <- DESeq(test_dds)
#   test_results <- getDeseqResults(test_dds, "coord_id")
#   return(test_results)
# }
# test_dds.ppmi <- testDDS(bsj_dds.ppmi) %>% left_join(coord_anno, by = "coord_id")
# write_csv(test_dds.ppmi, 'circRNA/output/ppmi_BSJ_results_limitedModel.csv')
# test_dds.icicle <- testDDS(bsj_dds.icicle) %>% left_join(coord_anno, by = "coord_id")
# write_csv(test_dds.icicle, 'circRNA/output/icicle_BSJ_results_limitedModel.csv')


# FUNCTIONAL ANALYSES -----------------------------------------------------

# significant BSJs
sig_bsj_enrich.ppmi <- goEnrich(
  geneSet = bsj_results.ppmi %>% filter(
    padj < 0.05,
    log2FoldChange > 0.1 | log2FoldChange < -0.1
  ) %>% pull(gene_id),
  background = bsj_results.ppmi %>% pull(gene_id),
  ontology = "ALL"
)
saveRDS(sig_bsj_enrich.ppmi, "circRNA/output/ppmi_sigBSJ_enrich.rds")


# ABUNDANT BSJs
abundant_bsj_enrich.ppmi <- goEnrich(
  geneSet = bsj_results.ppmi %>% pull(gene_id),
  background = junc_info %>% filter(study == "PPMI") %>% pull(gene_id) %>% unique(),
  ontology = "ALL"
)
saveRDS(abundant_bsj_enrich.ppmi, "circRNA/output/ppmi_abundantBSJ_enrich.rds")

abundant_bsj_enrich.icicle <- goEnrich(
  geneSet = bsj_results.icicle %>% pull(gene_id),
  background = junc_info %>% filter(study == "ICICLE-PD") %>% pull(gene_id) %>% unique(),
  ontology = "ALL"
)
saveRDS(abundant_bsj_enrich.icicle, "circRNA/output/icicle_abundantBSJ_enrich.rds")
