library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)


# IMPORT ------------------------------------------------------------------

# import metadata
meta_ppmi <- read_csv("ppmi/output/BL_iPD_HC_metadata.csv") %>%
  mutate(
    id = as.character(id),
    condition = factor(condition, levels = c("PD", "HC")),
    sex = factor(sex, levels = c("Male", "Female")),
    age_bin = factor(age_bin, levels = c("under_55", "55_to_65", "over_65")),
    plate = factor(plate),
    pct_usable_bases = scale(pct_usable_bases)
  )

meta_icicle <- read_csv("icicle/output/ICICLE_METADATA.csv") %>%
  mutate(
    condition = factor(condition, levels = c("PD", "HC")),
    sex = factor(sex, levels = c("Male", "Female")),
    age_bin = factor(age_bin, levels = c("under_55", "55_to_65", "over_65")),
    run = factor(run),
    pct_intronic_bases = scale(pct_usable_bases),
    pct_coding_bases = scale(pct_coding_bases)
  )

# IMPORT UNFILTERED AND DE DDS --------------------------------------------
dds_ppmi <- read_rds('linear/data/ppmi_DDS_prefilter.rds')
dds_icicle <- read_rds('linear/data/icicle_DDS_prefilter.rds')


# DESIGN --------------------------------------------------------------------------------------
# OG DESIGN JUST HAS CONDITION FOR QC AND SETUP
# specify design formula
design.ppmi <- ~ sex + plate + age_bin + pct_usable_bases + condition
design.icicle <- ~ sex + run + age_bin + pct_intronic_bases + pct_coding_bases + condition
# overwrite
dds_ppmi@design <- design.ppmi
dds_icicle@design <- design.icicle

colData(dds_ppmi)

# FILTERING ---------------------------------------------------------------
filterDDS <- function(dds, metadata){
  # apply raw count filtering
  keep <- rowSums(counts(dds, normalized = FALSE) >  10) > min(table(metadata$condition))
  # how many genes removed
  genesRemoved <- function(keep){
    total <- nrow(data.frame(keep))
    ntrue <- data.frame(logical = keep) %>% filter(logical == TRUE) %>% nrow()
    nfalse <- data.frame(logical = keep) %>% filter(logical == FALSE) %>% nrow()
    print(paste0("Total number of genes tested: ", total))
    print(paste0("Number of genes kept: ", ntrue, " (", round(ntrue/total*100), "%)"))
    print(paste0("Number of genes removed: ", nfalse, " (", round(nfalse/total*100), "%)"))
  }
  genesRemoved(keep)
  dds <- dds[keep, ]
  return(dds)
}

dds_ppmi <- filterDDS(dds_ppmi, meta_ppmi)
dds_icicle <- filterDDS(dds_icicle, meta_icicle)


# RUN OR LOAD DESEQ2 ------------------------------------------------------
getDDS <- function(dds, runOrImport, cohort){
  if (runOrImport == "run"){
    cat("Running DESeq2 (might take a while)\n")
    dds <- DESeq(dds)
    write_rds(dds, paste0("linear/data/dds_", cohort, ".rds"))
    return(dds)
  } else if (runOrImport == "import"){
    cat("Importing previously run DESeq2 analysis DeseqDataSet\n")
   dds <- read_rds(paste0("linear/data/dds_", cohort, ".rds"))
   return(dds)
  } else {
    stop("runOrLoad must be either run or import")
  }
}
dds_ppmi <- getDDS(dds_ppmi, "import", "ppmi")
dds_icicle <- getDDS(dds_icicle, "import", "icicle")


# RESULTS -----------------------------------------------------------------
contrast <- c("condition", "PD", "HC")

getResults <- function(dds){
  results <- results(dds, contrast = contrast, alpha = 0.05, independentFiltering = FALSE)
  print(summary(results))
  return(results)
}
addSymbolEntrez <- function(results){
  results$gene_name <- mapIds(org.Hs.eg.db,
                           keys = results$ensembl,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
  results$entrez <- mapIds(org.Hs.eg.db,
                           keys = results$ensembl,
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")
  return(results)
}

results_ppmi <- getResults(dds_ppmi) %>% 
  as.data.frame() %>% 
  rownames_to_column("ensembl")
results_ppmi <- addSymbolEntrez(results_ppmi)
write_csv(results_ppmi, "linear/output/ppmi_deseqResults.csv")

results_icicle <- getResults(dds_icicle) %>% 
  as.data.frame() %>% 
  rownames_to_column("ensembl")
results_icicle <- addSymbolEntrez(results_icicle)
write_csv(results_icicle, "linear/output/icicle_deseqResults.csv")

full_join(results_ppmi, results_icicle, by = c('ensembl', 'gene_name'), suffix = c('.ppmi', '.icicle')) %>%
  filter(padj.ppmi < 0.05,
         log2FoldChange.ppmi < -0.1 | log2FoldChange.ppmi > 0.1) %>%
  mutate(replicate_fdr = p.adjust(pvalue.icicle, method = 'fdr')) %>% view()

# EXPORT COUNTS -----------------------------------------------------------
exportFilteredCounts <- function(dds, cohort) {
  # raw counts
  counts(dds, normalized = FALSE) %>%
    as.data.frame() %>%
    rownames_to_column("ensembl") %>%
    write_csv(paste0("linear/data/", cohort, "_deseqFilteredRawCounts.csv"))
  cat(paste0(cohort, ": exported filtered raw counts"), "\n")
  # normalised counts
  counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("ensembl") %>%
    write_csv(paste0("linear/data/", cohort, "_deseqFilteredNormalisedCounts.csv"))
  cat(paste0(cohort, ": exported filtered normalised counts"), "\n")
  # VST counts
}
exportFilteredCounts(dds_ppmi, "ppmi")
exportFilteredCounts(dds_icicle, "icicle")


# FUNCTIONAL ANALYSIS -----------------------------------------------------
goEnrich <- function(geneSet, background, ontology){
  set.seed(1997)
  enrich <- enrichGO(gene = geneSet,
                     universe = background,
                     keyType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                     ont = ontology,
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
  return(enrich)
}

keggGSEA <- function(results){
  set.seed(1997)
  # remove genes with no entrez ID
  results <- filter(results, entrez != "NA")
  # remove duplicate entrez IDs
  results <- results[which(duplicated(results$entrez) == FALSE), ]
  # get fold changes
  fc <- results$log2FoldChange
  # name the fold changes with entrez id
  names(fc) <- results$entrez
  # sort fold changes in decreasing order
  fc <- sort(fc, decreasing = TRUE)
  # run GSEA
  kegg <- gseKEGG(geneList = fc, 
                  organism = "hsa",
                  #nPerm = 1000,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 1,
                  eps = 0,
                  seed = TRUE,
                  verbose = FALSE)
  return(kegg)
}

goGSEA <- function(results, ontology){
  set.seed(1997)
  # remove genes with no entrez ID
  results <- filter(results, entrez != "NA")
  # remove duplicate entrez IDs
  results <- results[which(duplicated(results$entrez) == FALSE), ]
  # get fold changes
  fc <- results$log2FoldChange
  # name the fold changes with entrez id
  names(fc) <- results$entrez
  # sort fold changes in decreasing order
  fc <- sort(fc, decreasing = TRUE)
  go <- gseGO(geneList = fc, 
              ont = ontology, 
              OrgDb = org.Hs.eg.db::org.Hs.eg.db,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 1,
              eps = 0,
              seed = TRUE,
              verbose = FALSE)
  return(go)
}

# GO enrichment
# go_enrich.ppmi <- list(
  # mf = goEnrich(
  #   geneSet = results_ppmi %>% filter(padj < 0.05) %>% pull(ensembl),
  #   background = results_ppmi %>% pull(ensembl),
  #   ontology = "MF"
  # ),
#   cc = goEnrich(
#     geneSet = results_ppmi %>% filter(padj < 0.05) %>% pull(ensembl),
#     background = results_ppmi %>% pull(ensembl),
#     ontology = "CC"
#   ),
#   bp = goEnrich(
#     geneSet = results_ppmi %>% filter(padj < 0.05) %>% pull(ensembl),
#     background = results_ppmi %>% pull(ensembl),
#     ontology = "BP"
#   )
# )
  go_enrich.ppmi = goEnrich(
    geneSet = results_ppmi %>% filter(padj < 0.05) %>% pull(ensembl),
    background = results_ppmi %>% pull(ensembl),
    ontology = "ALL"
  )
saveRDS(go_enrich.ppmi, "linear/output/ppmi_GO_enrichment.rds")

# go_enrich.icicle <- list(
#   mf = goEnrich(
#     geneSet = results_icicle %>% filter(padj < 0.05) %>% pull(ensembl),
#     background = results_icicle %>% pull(ensembl),
#     ontology = "MF"
#   ),
#   cc = goEnrich(
#     geneSet = results_icicle %>% filter(padj < 0.05) %>% pull(ensembl),
#     background = results_icicle %>% pull(ensembl),
#     ontology = "CC"
#   ),
#   bp = goEnrich(
#     geneSet = results_icicle %>% filter(padj < 0.05) %>% pull(ensembl),
#     background = results_icicle %>% pull(ensembl),
#     ontology = "BP"
#   )
# )
go_enrich.icicle <- goEnrich(
      geneSet = results_icicle %>% filter(padj < 0.05) %>% pull(ensembl),
      background = results_icicle %>% pull(ensembl),
      ontology = "ALL"
)
saveRDS(go_enrich.icicle, "linear/output/icicle_GO_enrichment.rds")


# Gene set enrichment analysis
# go_gsea.ppmi <- list(
#   mf = goGSEA(results_ppmi, "BP"),
#   cc = goGSEA(results_ppmi, "CC"),
#   bp = goGSEA(results_ppmi, "BP")
# )
go_gsea.ppmi <- goGSEA(results_ppmi, "ALL")
saveRDS(go_gsea.ppmi, "linear/output/ppmi_GO_gsea.rds")

# go_gsea.icicle <- list(
#   mf = goGSEA(results_icicle, "BP"),
#   cc = goGSEA(results_icicle, "CC"),
#   bp = goGSEA(results_icicle, "BP")
# )
go_gsea.icicle <- goGSEA(results_icicle, 'ALL')
saveRDS(go_gsea.icicle, "linear/output/icicle_GO_gsea.rds")

kegg_gsea.ppmi <- keggGSEA(results_ppmi)
saveRDS(kegg_gsea.ppmi, "linear/output/ppmi_KEGG_gsea.rds")
kegg_gsea.icicle <- keggGSEA(results_icicle)
saveRDS(kegg_gsea.icicle, "linear/output/icicle_KEGG_gsea.rds")
