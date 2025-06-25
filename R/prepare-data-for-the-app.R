# ==============================================================================
#
# R script for preparing data for the Shiny app
#
# Author: Jeffrey Y Lee
# Date: June 2025
#
# ==============================================================================

# ==============================================================================
# 1. Full data table
# ==============================================================================

library(tidyverse)
library(qs)

supp_table_raw <- qread("~/Documents/LCLIP/analysis_data/lclip_summary.qs") %>%
  map_dfr(~.x, .id = "exp_string")

supp_table_filtered <- supp_table_raw %>%
  dplyr::select(
    "stage_RBP" = exp_string,
    gene_id,
    gene_name,
    avg_iCLIP_tpm,
    "max_log2foldchange" = l2fc_bsmax,
    n_binding_site,
    binding_feature,
    human_homologs,
    mammalian_imp_target,
    mammalian_syp_target
  ) %>%
  mutate(avg_iCLIP_tpm = round(avg_iCLIP_tpm, digits = 3)) %>%
  mutate(max_log2foldchange = round(max_log2foldchange, digits = 3))

supp_table_filtered %>% saveRDS("./app/data/full_data.RDS")

# ==============================================================================
# 2. Binding site plots
# ==============================================================================

source("./R/plot_iclip_coverage.R")
supp_table_filtered <- readRDS("./app/data/full_data.RDS")

all_target_gene_names <- unique(supp_table_filtered$gene_name)
all_target_gene_ids <- unique(supp_table_filtered$gene_id)
saveRDS(all_target_gene_names, "./app/data/target_gene_names.RDS")

plot_iclip_coverage_wrap(gene = "mamo")

all_target_gene_names %>%
  purrr::set_names() %>%
  iwalk(~{
    plot_iclip_coverage_wrap(gene = .x) %>%
      patchworkGrob() %>%
      saveRDS(paste0(
        "./app/data/plots/",
        .y,
        ".RDS"
      ))
  }, .progress = TRUE) 

# ==============================================================================
# 3. Binding site table
# ==============================================================================

library(tidyverse)
library(qs)

lclip_bs <- qread("~/Documents/LCLIP/analysis_data/lclip_bs_list.qs")
lclip_overlap <- qread("~/Documents/LCLIP/analysis_data/lclip_intersect-feature_annotated.qs")

lclip_overlap[[1]] %>% View()
  


lclip_bs %>% names()














