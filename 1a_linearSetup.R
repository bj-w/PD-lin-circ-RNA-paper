library(tximport)
library(AnnotationHub)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(DESeq2)
library(tidyverse)


# FUNCTIONS ---------------------------------------------------------------
# imports salmon quant files using tximport
createTxi <- function(fileList, software) {
  # gene ID from transcript ID
  txdb <- makeTxDbFromGFF(file = "Homo_sapiens.GRCh38.101.gtf.gz")
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  tx2gene$TXNAME <- sub("\\.\\d+", "", tx2gene$TXNAME)
  tx2gene$GENEID <- sub("\\.\\d+", "", tx2gene$GENEID)
  # run tximport
  if (software == "deseq") {
    cat("Software is DESeq2 - raw counts \n")
    txi <- tximport(fileList,
      type = "salmon",
      tx2gene = tx2gene,
      ignoreAfterBar = TRUE,
      ignoreTxVersion = TRUE
    )
  } else if (software == "limma") {
    cat("Software is Limma - setting countsFromAbundance to lengthScaledTPM \n")
    txi <- tximport(fileList,
      type = "salmon",
      tx2gene = tx2gene,
      countsFromAbundance = "lengthScaledTPM",
      ignoreAfterBar = TRUE,
      ignoreTxVersion = TRUE
    )
  }
  return(txi)
}


# IMPORT METADATA AND FILELIST --------------------------------------------
# import metadata for each cohort
ppmi_meta <- read_csv("ppmi/output/BL_iPD_HC_metadata.csv") %>% 
  mutate(id = as.character(id),
  condition = factor(condition, levels = c("PD", "HC")),
  sex = factor(sex, levels = c("Male", "Female")),
  age_bin = factor(age_bin, levels = c("under_55", "55_to_65", "over_65")),
  plate = factor(plate))

icicle_meta <- read_csv("icicle/output/ICICLE_METADATA.csv") %>% 
  mutate(
    condition = factor(condition, levels = c("PD", "HC")),
    sex = factor(sex, levels = c("Male", "Female")),
    age_bin = factor(age_bin, levels = c("under_55", "55_to_65", "over_65")),
    run = factor(run))

# create a list of salmon files for each cohort to be imported by tximport
ppmi_files <- file.path("ppmi/data/raw/BL/salmon", ppmi_meta$path, "quant.sf")
names(ppmi_files) <- ppmi_meta$path
icicle_files <- file.path("icicle/data/raw/salmon", icicle_meta$id, "quant.sf")
names(icicle_files) <- icicle_meta$id


# CREATE TXI OBJECT FOR EACH COHORT ---------------------------------------

# limma
ppmi_txi_limma <- createTxi(ppmi_files, "limma")
write_rds(ppmi_txi_limma, "linear/data/ppmi_txiForLimma.rds")
icicle_txi_limma <- createTxi(icicle_files, "limma")
write_rds(icicle_txi_limma, "linear/data/icicle_txiForLimma.rds")

# deseq
ppmi_txi_deseq <- createTxi(ppmi_files, "deseq")
write_rds(ppmi_txi_deseq, "linear/data/ppmi_txiForDeseq.rds")
icicle_txi_deseq <- createTxi(icicle_files, "deseq")
write_rds(icicle_txi_deseq, "linear/data/icicle_txiForDeseq.rds")

# match orders of metadata and IDs in txi counts
ppmi_txi_deseq$counts <- ppmi_txi_deseq$counts[, match(ppmi_meta$path, colnames(ppmi_txi_deseq$counts))]
icicle_txi_deseq$counts <- icicle_txi_deseq$counts[, match(icicle_meta$id, colnames(icicle_txi_deseq$counts))]
# check they're the same
identical(colnames(ppmi_txi_deseq$counts), ppmi_meta$path)
identical(colnames(icicle_txi_deseq$counts), icicle_meta$id)

# rename PPMI sample names to just be the ID (ICICLE are fine)
colnames(ppmi_txi_deseq$counts) <- ppmi_meta$id


# CREATE DESEQ DATASET ----------------------------------------------------
# just set design to condition for now (before covariate selection)
dds_ppmi <- DESeqDataSetFromTximport(
  ppmi_txi_deseq,
  colData = ppmi_meta,
  design = ~ condition
)

dds_icicle <- DESeqDataSetFromTximport(
  txi = icicle_txi_deseq,
  colData = icicle_meta,
  design = ~ condition
)

# NORMALISE COUNTS --------------------------------------------------------
dds_ppmi <- estimateSizeFactors(dds_ppmi)
dds_icicle <- estimateSizeFactors(dds_icicle)

exportSizeFactors <- function(dds, study){
  nm <- assays(dds)[["avgTxLength"]]
  sf <- estimateSizeFactorsForMatrix(counts(dds) / nm)
  path <- paste0("linear/data/", study, "_linearSizeFactors.rds")
  cat(paste0('Writing ', study, ' size factors to: ', path, '\n'))
  write_rds(sf, path)
}

exportSizeFactors(dds_ppmi, 'ppmi')
exportSizeFactors(dds_icicle, 'icicle')

 
# export VST counts
# vst normalised counts
vst_ppmi <- assay(vst(dds_ppmi, blind = FALSE)) %>%
  as.data.frame() %>%
  rownames_to_column("ensembl")
write_csv(vst_ppmi, "linear/data/ppmi_vstCounts.csv")

vst_icicle <- assay(vst(dds_icicle, blind = FALSE)) %>%
  as.data.frame() %>%
  rownames_to_column("ensembl")
write_csv(vst_icicle, "linear/data/icicle_vstCounts.csv")

# export DDS objects (before filtering and carrrying out DE)
write_rds(dds_ppmi, 'linear/data/ppmi_DDS_prefilter.rds')
write_rds(dds_icicle, 'linear/data/icicle_DDS_prefilter.rds')