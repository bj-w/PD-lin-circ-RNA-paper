# consistent theme for plotting
plot_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = ggtext::element_markdown(face = 'plain', hjust = 0.5, size = 12),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.tag = element_text(size = 16, face = 'bold'),
    #strip.background = ggfun::element_roundrect(fill = 'gray93', colour = NA)
    strip.background = element_blank(),
    legend.position = "bottom")



# COLOURS -------------------------------------------------------------------------------------
condition_colours <- c("PD" = "#F7931E", "Control" = "#04545E")
study_colours <- c("PPMI" = "#CA5435", "ICICLE-PD" = "#35ABCA")


# FORMATTING ----------------------------------------------------------------------------------


# Converts long sample IDs in the column 'id' to shorter ones:
# PPMI-Phase1-IR1_3432_BL_PP0015-9806_5104-SL-1841_R1.fastq.gz
# is converted to
# 3432
cleanPPMIsampleID <- function(dfToClean){
  metadata <- data.table::fread(here::here('ppmi/output/BL_iPD_HC_metadata.csv'))
  ids <- metadata %>%
    select(id, fastq_name)
  clean <- dfToClean %>%
    rename(fastq_name = id) %>%
    left_join(ids, by = 'fastq_name') %>%
    relocate(id, .before = fastq_name) %>%
    mutate(id = as.character(id)) %>%
    select(!fastq_name)
  return(clean)
}

# CONVERTS EXPONENTION TO SUPERSCRIPT
expSup <- function(w, digits=1) {
  sprintf(paste0("%.", digits, "fx10<sup>%d</sup>"), w/10^floor(log10(abs(w))), floor(log10(abs(w))))
}

# FUNCTIONAL ANALYSES -------------------------------------------------------------------------

goEnrich <- function(geneSet, background, ontology) {
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

keggGSEA <- function(results) {
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

goGSEA <- function(results, ontology) {
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

addEntrezID <- function(results) {
  results$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                           keys = results$gene_id,
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first"
  )
  return(results)
}

# script specific functions
conditionBoxplot <- function(df, y) {
  df %>%
    ggplot(aes(x = condition, y = {{ y }})) +
    geom_violin(aes(fill = condition), alpha = 0.1, colour = NA) +
    geom_boxplot(aes(colour = condition), width = 0.2, lwd = 1, fill = "white") +
    scale_colour_manual(values = condition_colours) +
    scale_fill_manual(values = condition_colours) +
    xlab("Study group") +
    labs(
      fill = "Study group",
      colour = "Study group"
    ) +
    # stat_compare_means(comparisons = list(c("PD", "Control")), method = "wilcox.test",
    #                    aes(label = "..p.format"), label.x = 0.6, size = 3) +
    facet_grid(~study)
}

plotVolcano <- function(results, sigResults) {
  colnames(results)[1] <- "feature"
  colnames(sigResults)[1] <- "feature"
  results %>%
    mutate(colour = if_else(feature %in% sigResults$feature, "DE", "Not DE")) %>%
    arrange(desc(colour)) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue), colour = colour)) +
    geom_point(alpha = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = "log<sub>2</sub>(Fold Change)",
      y = "-log<sub>10</sub>(<i>P</i>)",
      colour = "Differential expression in PPMI"
    ) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", alpha = 0.2) +
    facet_grid(~study)
}

getTPM <- function() {
  formatTPM <- function(txi, study) {
    txi[["abundance"]] %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      pivot_longer(cols = !gene_id, names_to = "fastq_name", values_to = "TPM") %>%
      mutate(study = {{ study }})
  }
  cat("Importing PPMI TPMs...\n")
  linear_meta.ppmi <- read_csv(here("PPMI_RAW_DATA/metadata/output/BL_iPD_HC_metadata.csv")) %>%
    mutate(id = as.character(id))
  ppmi <- readRDS(here("linear/data/ppmi_txiForDeseq.rds")) %>%
    formatTPM("PPMI") %>%
    # add on sample ID that corresponds to fastq full name
    left_join(linear_meta.ppmi[, c("id", "fastq_name")], by = "fastq_name") %>%
    # remove fastq name path by selecting only the specific columns we need
    select(gene_id, id, TPM, study)
  cat("Importing ICICLE-PD TPMs...\n")
  icicle <- readRDS(here("linear/data/icicle_txiForDeseq.rds")) %>%
    formatTPM("ICICLE-PD") %>%
    rename("id" = "fastq_name")
  tpm <- bind_rows(ppmi, icicle)
  return(tpm)
}

# set helpful variables
cohorts <- c("PPMI", "ICICLE-PD")
rna_colours <- c(
  "BSJ" = "#1b9e77",
  "FSJ" = "#d95f02",
  "BSJ:FSJ ratio" = "#a6761d",
  "BSJ hosts (Gene)" = "#e7298a",
  "Gene" = "#7570b3",
  "Gene and BSJ" = "#e6ab02",
  "Imbalance" = "#66a61e"
)
