library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(data.table)
library(fst)
library(tidyverse)


# import common functions
source("combinedLinCirc/PD-RNA/functions.R")

# IMPORT METADATA ---------------------------------------------------------

# Converts long sample IDs in the column 'id' to shorter ones:
## PPMI-Phase1-IR1_3432_BL_PP0015-9806_5104-SL-1841_R1.fastq.gz
# is converted to
## 3432
cleanPPMIsampleID <- function(dfToClean){
  metadata <- data.table::fread('ppmi/output/BL_iPD_HC_metadata.csv')
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

meta.ppmi <- fread("ppmi/output/BL_iPD_HC_metadata.csv") %>%
  rename(
    batch = plate,
    updrs3 = updrs3_score,
  ) %>%
  mutate(
    id = as.character(id),
    sex = factor(sex),
    age_bin = factor(age_bin),
    condition = recode(condition, "HC" = "Control"),
    condition = factor(condition, levels = c("Control", "PD")),
    batch = factor(batch),
    study = "PPMI"
  )

meta.icicle <- fread("icicle/output/ICICLE_METADATA.csv") %>%
  rename(
    batch = run,
    age_at_consent = age,
    updrs3 = mds_updrs_3_total,
    moca = mo_ca_total,
  ) %>%
  mutate(
    id = as.character(id),
    sex = factor(sex),
    age_bin = factor(age_bin),
    condition = recode(condition, "HC" = "Control"),
    condition = factor(condition, levels = c("Control", "PD")),
    batch = factor(batch),
    study = "ICICLE-PD"
  )

# COLLECT ALIGNMENT INFO --------------------------------------------------
# function to collect ciriquant alignment info
collectAlignInfo <- function(gtfDir) {
  # read in gtf file headers
  align_info <- map(list.files(gtfDir, full.names = TRUE),
    fread,
    nrows = 3
  )
  # tidy up alignment metric names
  align_info <- map(align_info, ~ .x %>%
    mutate(
      `##Sample:` = gsub(`##Sample:`, pattern = "##", replacement = ""),
      `##Sample:` = gsub(`##Sample:`, pattern = ":", replacement = "")
    ))
  # format metrics in table and combine all samples together
  align_info <- map_dfr(align_info, function(df) {
    data.frame(
      id = names(df)[2],
      total_reads = as.numeric(df[1, 2]),
      mapped_reads = as.numeric(df[2, 2]),
      circular_reads = as.numeric(df[3, 2])
    )
  })
  return(align_info)
}

# collect ciriquant alignment info
align.ppmi <- collectAlignInfo("circRNA/data/raw/rerun/ppmi/ciriquant/gtf/") %>%
  cleanPPMIsampleID()
align.icicle <- collectAlignInfo("circRNA/data/raw/rerun/icicle/ciriquant/gtf/")


# add alignment metrics onto metadata
meta.ppmi <- left_join(meta.ppmi, align.ppmi, by = "id")
meta.icicle <- left_join(meta.icicle, align.icicle, by = "id")


# JUNCTION INFO -----------------------------------------------------------

# function to create the table
createJunctionInfo <- function(gtfDir) {
  # use rtracklayer to import CIRIquant output as it's in gtf format
  cat("Importing GTFs...\n")
  juncInfo <- map(list.files(gtfDir, full.names = TRUE), rtracklayer::import)
  # convert to df
  juncInfo <- map(juncInfo, as.data.frame)
  # get ID names
  sample_names <- gsub(
    pattern = ".gtf",
    replacement = "",
    x = list.files(gtfDir)
  )
  # add ID as a name to each list df
  names(juncInfo) <- sample_names
  # Merge into one big df adding sample ID as a column
  cat("Merging into one DF...\n")
  output <- map_df(juncInfo, data.frame, .id = "id")
  # format output
  output <- output %>%
    # rename columns to make more consistent
    rename(chr = seqnames, genomic_span = width, CPM = score) %>%
    # drop unneeded columns
    select(-c(source, type, phase)) %>%
    # reformat circ_id into coord_id (chr:start-end:strand)
    unite("coord_id", c(circ_id, strand), sep = ":", remove = FALSE) %>% # add strand to the end
    mutate(coord_id = gsub(coord_id, pattern = "\\|", replacement = "-")) %>% # replace | with -
    # drop the now redundant circ_id column
    select(!circ_id)
  return(output)
}


# create the junction info df
junc_info.ppmi <- createJunctionInfo("circRNA/data/raw/rerun/ppmi/ciriquant/gtf/")
junc_info.ppmi <- junc_info.ppmi %>%
  cleanPPMIsampleID() %>%
  mutate(
    study = "PPMI",
    bsj = as.numeric(bsj),
    fsj = as.numeric(fsj),
    junc_ratio = as.numeric(junc_ratio)
  )

junc_info.icicle <- createJunctionInfo("circRNA/data/raw/rerun/icicle/ciriquant/gtf/") %>%
  mutate(
    study = "ICICLE-PD",
    bsj = as.numeric(bsj),
    fsj = as.numeric(fsj),
    junc_ratio = as.numeric(junc_ratio)
  )

# add on exon counts
## import exon counts
exons.ppmi <- read_delim("circRNA/data/exonCount_ppmi.bed", col_names = FALSE) %>%
  select(coord_id = X4, exonCount = X7)

exons.icicle <- read_delim("circRNA/data/exonCount_icicle.bed", col_names = FALSE) %>%
  select(coord_id = X4, exonCount = X7)

## add onto junc_info files
junc_info.ppmi <- left_join(junc_info.ppmi, exons.ppmi, by = "coord_id")
junc_info.icicle <- left_join(junc_info.icicle, exons.icicle, by = "coord_id")

# add on condition
junc_info.ppmi <- left_join(junc_info.ppmi, meta.ppmi[, c("id", "condition")], by = "id")
junc_info.icicle <- left_join(junc_info.icicle, meta.icicle[, c("id", "condition")], by = "id")


# BSJ-SAMPLE METRICS ------------------------------------------------------
# function to calculate sample specific BSJ metrics
getBSJsampleMetrics <- function(juncInfo, metadata) {
  # unique BSJ and summed backsplice CPM
  uniqueBSJ_and_sumBSJ <- juncInfo %>%
    group_by(id) %>%
    summarise(
      sum_bsj = sum(bsj),
      unique_bsj = n_distinct(coord_id)
    )
  # number of unique BSJ host genes (MULTIMAPPED BSJs ARE REMOVED!!!)
  unique_hostBSJ <- juncInfo %>%
    filter(!grepl(",", gene_name)) %>%
    group_by(id) %>%
    summarise(unique_hostBSJ = n_distinct(gene_id))
  # add onto metadata
  metadata <- left_join(metadata, uniqueBSJ_and_sumBSJ, by = "id")
  metadata <- left_join(metadata, unique_hostBSJ, by = "id")
  # normalise unique BSJ + unique BSJ host genes by mapped reads
  metadata <- metadata %>%
    mutate(
      sum_bsj_perMapped = (sum_bsj / salmon_mapped) * 1e+06,
      unique_bsj_perMapped = (unique_bsj / salmon_mapped) * 1e+06,
      unique_hostBSJ_perMapped = (unique_hostBSJ / salmon_mapped) * 1e+06
    )
  return(metadata)
}

# run function adding BSJ metrics to sample metadata
meta.ppmi <- getBSJsampleMetrics(juncInfo = junc_info.ppmi, metadata = meta.ppmi)
meta.icicle <- getBSJsampleMetrics(juncInfo = junc_info.icicle, metadata = meta.icicle)

# remove samples not ran through BSJ detection pipeline to make it easier to
# use the metadata file downstream
meta.ppmi <- meta.ppmi %>% filter(id %in% unique(junc_info.ppmi$id))
meta.icicle <- meta.icicle %>% filter(id %in% unique(junc_info.icicle$id))


# ADD COUNTS NORMALISED TO SALMON MAPPED TO JUNC INFO ---------------------
junc_info.ppmi <- junc_info.ppmi %>% 
  left_join(meta.ppmi[, c('id', 'salmon_mapped')], by = 'id') %>% 
  mutate(bsj_perMapped = (bsj / salmon_mapped) * 1e+06,
         fsj_perMapped = (fsj / salmon_mapped) * 1e+06) %>% 
  # drop salmon mapped reads to make junc file a bit smaller
  select(-salmon_mapped)

junc_info.icicle <- junc_info.icicle %>% 
  left_join(meta.icicle[, c('id', 'salmon_mapped')], by = 'id') %>% 
  mutate(bsj_perMapped = (bsj / salmon_mapped) * 1e+06,
         fsj_perMapped = (fsj / salmon_mapped) * 1e+06) %>% 
  # drop salmon mapped reads to make junc file a bit smaller
  select(-salmon_mapped)

# EXPORT METADATA ---------------------------------------------------------
write_rds(meta.ppmi, "circRNA/data/ppmi_metadata.rds")
write_rds(meta.icicle, "circRNA/data/icicle_metadata.rds")


# EXPORT JUNCTION INFO ----------------------------------------------------
junc_info.bound <- bind_rows(junc_info.ppmi, junc_info.icicle)
write_fst(junc_info.bound, "circRNA/data/bound_juncInfo.fst")



# EXPORT JUNCTION LEVEL MATRICES ------------------------------------------
# function to create count matrix for BSJ + FSJ counts
createJuncMatrix <- function(juncInfo, junctionType) {
  juncInfo %>% pivot_wider(id_cols = coord_id, names_from = id, values_from = {{ junctionType }}, values_fill = 0)
}


# create + export the BSJ count matrices
bsj_matrix.ppmi <- createJuncMatrix(junc_info.ppmi, "bsj")
write_fst(bsj_matrix.ppmi, "circRNA/data/ppmi_BSJrawCounts.fst")
bsj_matrix.icicle <- createJuncMatrix(junc_info.icicle, "bsj")
write_fst(bsj_matrix.icicle, "circRNA/data/icicle_BSJrawCounts.fst")

# create + export the FSJ count matrices
fsj_matrix.ppmi <- createJuncMatrix(junc_info.ppmi, "fsj")
write_fst(fsj_matrix.ppmi, "circRNA/data/ppmi_FSJrawCounts.fst")
fsj_matrix.icicle <- createJuncMatrix(junc_info.icicle, "fsj")
write_fst(fsj_matrix.icicle, "circRNA/data/icicle_FSJrawCounts.fst")



# SET UP DDS OBJECTS ------------------------------------------------------

# match order of metadata and counts
matchMetaCountMatrix <- function(countMatrix, metadata) {
  # move coord id into the rowname
  countMatrix <- countMatrix %>% column_to_rownames('coord_id')
  # match orders of metadata and IDs in txi countMatrix
  countMatrix <- countMatrix[, match(metadata$id, colnames(countMatrix))]
  # check they're the same
  same <- identical(colnames(countMatrix), metadata$id)
  stopifnot(same)
  cat(paste0("Do colnames of count matrix and metadata IDs match? ", same, "\n"))
  # return count matrix
  return(countMatrix)
}

bsj_matrix.ppmi <- matchMetaCountMatrix(countMatrix = bsj_matrix.ppmi, metadata = meta.ppmi)
bsj_matrix.icicle <- matchMetaCountMatrix(countMatrix = bsj_matrix.icicle, metadata = meta.icicle)

fsj_matrix.ppmi <- matchMetaCountMatrix(countMatrix = fsj_matrix.ppmi, metadata = meta.ppmi)
fsj_matrix.icicle <- matchMetaCountMatrix(countMatrix = fsj_matrix.icicle, metadata = meta.icicle)


# CREATE DESEQDATASETS
# specify design formula (just condition before QC)
design.ppmi <- ~ condition
design.icicle <- ~ condition

# create DDS
bsj_dds.ppmi <- DESeqDataSetFromMatrix(
  countData = bsj_matrix.ppmi,
  colData = meta.ppmi,
  design = design.ppmi
)

bsj_dds.icicle <- DESeqDataSetFromMatrix(
  countData = bsj_matrix.icicle,
  colData = meta.icicle,
  design = design.icicle
)

fsj_dds.ppmi <- DESeqDataSetFromMatrix(
  countData = fsj_matrix.ppmi,
  colData = meta.ppmi,
  design = design.ppmi
)

fsj_dds.icicle <- DESeqDataSetFromMatrix(
  countData = fsj_matrix.icicle,
  colData = meta.icicle,
  design = design.icicle
)


# NORMALISATION BASED ON LINEAR SIZE FACTORS
# get size factors from linear RNA quantification
sf.ppmi <- read_rds("linear/data/ppmi_linearSizeFactors.rds")
sf.ppmi <- sf.ppmi[names(bsj_matrix.ppmi)] # ensure only samples included in BSJ detection are there

sf.icicle <- read_rds("linear/data/icicle_linearSizeFactors.rds")
sf.icicle <- sf.icicle[names(bsj_matrix.icicle)] # ensure only samples included in BSJ detection are there

# add size factors onto DDS objects
sizeFactors(bsj_dds.ppmi) <- sf.ppmi
sizeFactors(bsj_dds.icicle) <- sf.icicle

sizeFactors(fsj_dds.ppmi) <- sf.ppmi
sizeFactors(fsj_dds.icicle) <- sf.icicle


# SAVE PRE FILTERED AND DE DDS
write_rds(bsj_dds.ppmi, 'circRNA/data/ppmi_bsj_DDS_preFilter.rds')
write_rds(bsj_dds.icicle, 'circRNA/data/icicle_bsj_DDS_preFilter.rds')
write_rds(fsj_dds.ppmi, 'circRNA/data/ppmi_fsj_DDS_preFilter.rds')
write_rds(fsj_dds.icicle, 'circRNA/data/icicle_fsj_DDS_preFilter.rds')