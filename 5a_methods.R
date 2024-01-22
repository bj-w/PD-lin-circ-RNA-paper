
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


# import data -------------------------------------------------------------

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

# filter TPM counts to reduce file size
tpm_counts <- getTPM()

### circular rna
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


# participant characteristics ---------------------------------------------
meta.icicle %>% 
  group_by(condition) %>% 
  summarise(mean = mean(moca, na.rm = TRUE),
            sd = sd(moca, na.rm = TRUE))

# differences
## sex
fisher.test(table(meta.ppmi$sex, meta.ppmi$condition))
fisher.test(table(meta.icicle$sex, meta.icicle$condition))
## age at sample collection
t.test(meta.ppmi$age_at_consent ~ meta.ppmi$condition)
t.test(meta.icicle$age_at_consent ~ meta.icicle$condition)
## updrs3
t.test(meta.ppmi$updrs3 ~ meta.ppmi$condition)
## moca
t.test(meta.ppmi$moca ~ meta.ppmi$condition)
t.test(meta.icicle$moca ~ meta.icicle$condition)
# seq depth ---------------------------------------------------------------
seq_depth <- bind_rows(
  fread(here("circRNA/data/raw/rerun/ppmi/qc/untrimmed/multiqc_data/multiqc_fastqc.txt")) %>%
    janitor::clean_names() %>%
    rename(id = sample) %>%
    mutate(id = gsub(id, pattern = "_R[1-2]", replacement = "")) %>%
    cleanPPMIsampleID() %>%
    select(id, total_sequences) %>%
    distinct() %>%
    mutate(study = "PPMI"),
  fread(here("circRNA/data/raw/rerun/icicle/qc/untrimmed/multiqc_data/multiqc_fastqc.txt")) %>%
    janitor::clean_names() %>%
    rename(id = sample) %>%
    mutate(id = gsub(id, pattern = "_R[1-2]", replacement = "")) %>%
    select(id, total_sequences) %>%
    distinct() %>%
    mutate(study = "ICICLE-PD")
) %>% mutate(study = factor(study, levels = c("PPMI", "ICICLE-PD")))

# stats
seq_depth %>%
  group_by(study) %>%
  summarise(
    median_seqDepth = round(median(total_sequences) / 1e+06, 3),
    iqr_seqDepth = round(IQR(total_sequences) / 1e+06, 3)
  )

# plot
seq_depth %>%
  ggplot(aes(total_sequences, fill = study)) +
  geom_density(alpha = 0.3, colour = NA) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  scale_fill_manual(values = study_colours) +
  labs(
    x = "Paired-end reads (Million)",
    y = "Density",
    fill = "Cohort"
  ) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1))
ggsave(here("combinedLinCirc/output/figures/supp/cohortSequencingDepths.png"),
       device = agg_png,
       height = 3, width = 5, dpi = 600)
ggsave(here("combinedLinCirc/output/figures/supp/cohortSequencingDepths.pdf"),
       height = 3, width = 5, dpi = 600)


# circrna detection -------------------------------------------------------
# Overlap in the circRNA detection between CIRI2, CIRCexplorer2 and PFv2
# Get all BSJs detected from each tool
detectedBsj <- function(study) {
  # get BSJs from CIRI2
  cat("Getting CIRI2 BSJs...\n")
  ciri <- map(list.files(paste0("circRNA/data/raw/rerun/", study, "/tool_output/ciri/"), full.names = TRUE), fread) %>%
    list_rbind()
  ciri <- paste0(ciri$chr, ":", ciri$circRNA_start - 1, "-", ciri$circRNA_end, ":", ciri$strand)
  ciri <- unique(ciri)
  # get BSJs from CIRCexplorer2
  cat("Getting CIRCexplorer2 BSJs...\n")
  ce <- map(list.files(paste0("circRNA/data/raw/rerun/", study, "/tool_output/ce/"), full.names = TRUE), fread) %>%
    list_rbind()
  ce <- paste0(ce$V1, ":", ce$V2, "-", ce$V3, ":", ce$V6)
  ce <- unique(ce)
  # get BSJs from PFv2
  cat("Getting PFv2 BSJs...\n")
  pf <- map(list.files(paste0("circRNA/data/raw/rerun/", study, "/tool_output/pf/"), full.names = TRUE), fread) %>%
    list_rbind()
  pf <- paste0(pf$V1, ":", pf$V2, "-", pf$V3, ":", pf$V6)
  pf <- unique(pf)
  bsjs <- list(ciri = ciri, ce = ce, pf = pf)
}

detected_bsj.ppmi <- detectedBsj(study = "ppmi")
detected_bsj.icicle <- detectedBsj(study = "icicle")

# how many detected by each tool
map(detected_bsj.ppmi, length)
map(detected_bsj.icicle, length)

# how many detected per cohort
reduce(detected_bsj.ppmi, c) %>%
  unique() %>%
  length()
reduce(detected_bsj.icicle, c) %>%
  unique() %>%
  length()

# Venn diagrams
bsj_venn.ppmi <- ggVennDiagram(
  x = list(
    CIRI2 = detected_bsj.ppmi[["ciri"]],
    CIRCexplorer2 = detected_bsj.ppmi[["ce"]],
    PFv2 = detected_bsj.ppmi[["pf"]]
  ),
  set_color = "black"
) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  scale_colour_manual(values = c("black", "black", "black")) +
  ggtitle("PPMI")

bsj_venn.icicle <- ggVennDiagram(
  x = list(
    CIRI2 = detected_bsj.icicle[["ciri"]],
    CIRCexplorer2 = detected_bsj.icicle[["ce"]],
    PFv2 = detected_bsj.icicle[["pf"]]
  ),
  set_color = "black"
) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  scale_colour_manual(values = c("black", "black", "black")) +
  ggtitle("ICICLE-PD")

bsj_venn.ppmi | bsj_venn.icicle
ggsave(here("combinedLinCirc/output/figures/supp/bsj_venn.svg"),
       height = 8, width = 12
)