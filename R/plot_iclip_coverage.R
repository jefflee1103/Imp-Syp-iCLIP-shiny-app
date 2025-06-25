library(tidyverse)
library(ggtranscript)
library(patchwork)
library(valr)
library(qs)


#  - - - - - Flybase R6.28 (~= ENSEMBL v99) gtf 

## Get gtf
fb628_gtf <- valr::read_gtf("~/Documents/LCLIP/GENOMEDIR/dmel-all-r6.28.gtf")
## Get data
# lclip_rpm_list <- qread("~/Documents/LCLIP/analysis_data/lclip_rpm_list.qs")
lclip_bs_list <- qread("~/Documents/LCLIP/analysis_data/lclip_bs_list.qs")
lclip_rpm_peaks_list <- qread("~/Documents/LCLIP/analysis_data/lclip_rpm_peaks_list.qs")

# gene region function
get_gene_region <- function(gtf, gene){
  if(gene %in% gtf$gene_symbol){
    gtf %>% dplyr::filter(gene_symbol == gene)
  } else{
    print(paste0(gene, " does not exist in FlyBase R6.28"))
  }
}

# - - - - - Gene annotation

# plot gene annotation function
plot_gene_annotation <- function(
    gtf = fb628_gtf, 
    gene, 
    tx_order = NULL,
    cds_colour = "#d1ab75", 
    utr_colour = "gray80", 
    intron_colour = "gray30", 
    txid_text_colour = "gray50"
){
  
  # gene region function
  get_gene_region <- function(gtf, input_gene){
    if(input_gene %in% gtf$gene_symbol){
      gtf %>% filter(gene_symbol == input_gene)
    } else{
      print(paste0(input_gene, " does not exist in FlyBase R6.28"))
    }
  }
  
  # feature dfs
  if(is.null(tx_order)){
    gene_region <- get_gene_region(gtf, gene) %>%
      dplyr::filter(type != "start_codon") %>%
      dplyr::filter(type != "stop_codon")
  } else {
    gene_region <- get_gene_region(gtf, gene) %>%
      dplyr::filter(type != "start_codon") %>%
      dplyr::filter(type != "stop_codon") %>%
      dplyr::filter(transcript_symbol %in% tx_order) %>%
      mutate(transcript_symbol = fct_relevel(transcript_symbol, rev(tx_order)))
  }
  
  exons <- gene_region %>% dplyr::filter(type %in% c("exon", "miRNA", "pre_miRNA"))
  introns <- to_intron(exons, "transcript_symbol") 
  cds <- gene_region %>% dplyr::filter(type %in% c("CDS", "miRNA"))
  gene_start <- gene_region %>% pull(start) %>% min()
  gene_end <- gene_region %>% pull(end) %>% max()
  
  #Plot
  exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_symbol
    )) +
    geom_intron(
      data = to_intron(exons, "transcript_symbol"),
      aes(strand = strand),
      colour = intron_colour,
      arrow = grid::arrow(ends = "last", length = grid::unit(0.04, "inches")),
      arrow.min.intron.length = 100,
    ) +
    geom_range(
      fill = utr_colour, 
      height = 0.2,
      size = 0.1
    ) +
    geom_range(
      data = cds,
      fill = cds_colour,
      size = 0.1
    ) +
    coord_cartesian(xlim = c(gene_start - 100, gene_end + 100)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    theme_void(base_size = 15) +
    theme(axis.text.y = element_text(size = 14, colour = txid_text_colour, face = "italic")) 
}


# - - - - - iCLIP coverage

# # create merged df of rpm
# exp_string <- c("L1_Imp", "L2_Imp", "L2_Syp", "L3_Syp") %>% purrr::set_names()
# 
# rpm_list <- exp_string %>%
#   purrr::map(~{
#     list.files("~/Documents/LCLIP/Xlsite/merged_reps/shifted/bedgraph_sorted", pattern = .x, full.names = TRUE) %>%
#       .[!str_detect(., "SMI")] %>%
#       purrr::set_names() %>%
#       map_dfr(read_bedgraph, .id = "strand") %>%
#       mutate(strand = if_else(str_detect(strand, ".\\+.sorted.bedgraph"), "+", "-")) %>%
#       dplyr::select(chrom, start, end, strand, value)
#   })
# 
# qsave(rpm_list, "~/Documents/LCLIP/analysis_data/lclip_rpm_list.qs")
# 
# # create list of binding sites
# bs_list <- exp_string %>%
#   purrr::map(~{
#     read_bed(
#       paste0("~/Documents/LCLIP/iCount/clusters/filtered_clusters/", .x, "_peaks_clustered_filtered_enriched_over_input_curated.bed"),
#       n_fields = 6
#     )
#   })
# 
# qsave(bs_list, "~/Documents/LCLIP/analysis_data/lclip_bs_list.qs")



## more helper functions
plot_rpm <- function(df, plot_colour, y_label, rpm_max, gene_region, col_width){
  ggplot(df, aes(x = start, y = value, group = chrom)) +
    # geom_area(position = 'identity', fill = plot_colour) +
    # geom_line(colour = plot_colour, size = 0.2) +
    geom_col(width = col_width, fill = plot_colour) + 
    geom_hline(yintercept = 0, colour = "gray40", size = 0.2) + 
    coord_cartesian(xlim = c(min(gene_region$start) - 100, max(gene_region$end) + 100)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, rpm_max), breaks = c(0, rpm_max), labels = ~{paste(.x, "-")}) +
    labs(y = y_label) + 
    theme_void(base_size = 15) +
    theme(
      axis.text.y = element_text(size = 14, colour = "gray10", hjust = 1),
      axis.title.y = element_text(size = 15, angle = 90)
    )
}

plot_bs <- function(df, plot_colour, gene_region, bs_width){
  if(nrow(df) > 0){
    ggplot(df, aes(xstart = start, xend = end, y = "")) +
      geom_range(
        fill = plot_colour,
        colour = plot_colour,
        height = 2,
        size = bs_width
      ) +
      coord_cartesian(xlim = c(min(gene_region$start) - 100, max(gene_region$end) + 100)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_void(base_size = 15) 
  } else {
    ggplot() + theme_void()
  }
}

# plot iclip coverage function
plot_iclip_coverage <- function(
    gtf = fb628_gtf,
    rpm_list = lclip_rpm_peaks_list,
    bs_list = lclip_bs_list,
    gene,
    tx_order = NULL,
    custom_rpm_max = FALSE,
    custom_rpm_max_value = 20,
    colours = c("#406e3c", "#406e3c", "#83396d", "#83396d"),
    col_width = NULL,
    bs_width = 0
){
  # gene region
  if(is.null(tx_order)){
    gene_region <- get_gene_region(gtf = gtf, gene = gene) %>%
      group_by(strand) %>%
      bed_merge()
  } else {
    gene_region <- get_gene_region(gtf = gtf, gene = gene) %>%
      dplyr::filter(transcript_symbol %in% tx_order) %>%
      group_by(strand) %>%
      bed_merge()
  }
  
  # col_width
  col_width <- ceiling((max(gene_region$end) - min(gene_region$start)) / 500)
  
  # get rpm df
  rpm_df_list <- rpm_list %>%
    purrr::map(~{
      .x %>%
        group_by(strand) %>%
        bed_intersect(gene_region, suffix = c("", ".genome")) %>%
        filter(.overlap > 0)
    })
  
  # get max rpm to normalise all exp
  if(custom_rpm_max == FALSE){
    rpm_max <- rpm_df_list %>%
      purrr::map_dfr(~.x) %>% 
      pull(value) %>% max() %>% ceiling()
  } else {
    rpm_max <- custom_rpm_max_value
  }
  
  # get bs df
  bs_df_list <- bs_list %>%
    purrr::map(~{
      .x %>%
        group_by(strand) %>%
        bed_intersect(gene_region, suffix = c("", ".genome")) %>%
        filter(.overlap > 0) 
    })
  
  # colours
  colours <- colours 
  
  # L1 Imp 
  L1_Imp_rpm <- plot_rpm(
    df = rpm_df_list$L1_Imp,
    plot_colour = colours[1],
    y_label = "L1 Imp",
    rpm_max = rpm_max,
    gene_region = gene_region,
    col_width = col_width
  )
  
  L1_Imp_bs <- plot_bs(df = bs_df_list$L1_Imp, plot_colour = colours[1], gene_region = gene_region, bs_width = bs_width)
  
  # L2 Imp 
  L2_Imp_rpm <- plot_rpm(
    df = rpm_df_list$L2_Imp,
    plot_colour = colours[2],
    y_label = "L2 Imp",
    rpm_max = rpm_max,
    gene_region = gene_region,
    col_width = col_width
  )
  
  L2_Imp_bs <- plot_bs(df = bs_df_list$L2_Imp, plot_colour = colours[2], gene_region = gene_region, bs_width = bs_width)
  
  # L2 Syp 
  L2_Syp_rpm <- plot_rpm(
    df = rpm_df_list$L2_Syp,
    plot_colour = colours[3],
    y_label = "L2 Syp",
    rpm_max = rpm_max,
    gene_region = gene_region,
    col_width = col_width
  )
  
  L2_Syp_bs <- plot_bs(df = bs_df_list$L2_Syp, plot_colour = colours[3], gene_region = gene_region, bs_width = bs_width)
  
  # L3 Syp 
  L3_Syp_rpm <- plot_rpm(
    df = rpm_df_list$L3_Syp,
    plot_colour = colours[4],
    y_label = "L3 Syp",
    rpm_max = rpm_max,
    gene_region = gene_region,
    col_width = col_width
  )
  
  L3_Syp_bs <- plot_bs(df = bs_df_list$L3_Syp, plot_colour = colours[4], gene_region = gene_region, bs_width = bs_width)
  
  # transcript annotation
  gene_annotation <- plot_gene_annotation(gtf = gtf, gene = gene, tx_order = tx_order)
  
  # output
  output <- list(
    L1_Imp_rpm = L1_Imp_rpm,
    L1_Imp_bs = L1_Imp_bs,
    L2_Imp_rpm = L2_Imp_rpm,
    L2_Imp_bs = L2_Imp_bs,
    L2_Syp_rpm = L2_Syp_rpm,
    L2_Syp_bs = L2_Syp_bs,
    L3_Syp_rpm = L3_Syp_rpm,
    L3_Syp_bs = L3_Syp_bs,
    gene_annotation = gene_annotation
  )
  
  return(output)
}

# usage 
# plot_iclip_coverage(
#   gtf = fb628_gtf,
#   rpm_list = rpm_list,
#   bs_list = bs_list,
#   gene = "pros",
#   tx_order = c("pros-RK", "pros-RL"), 
#   custom_rpm_max = FALSE,
#   custom_rpm_max_value = 17
# ) %>%
#   wrap_plots(
#     ncol = 1,
#     heights = c(rep(c(12, 0.8), 4), 10)
#   )
  # & coord_cartesian(xlim = c(10486561, 10528425))


# plot and quick wrap 

plot_iclip_coverage_wrap <- function(
    gene,
    tx_order = NULL,
    start = NULL,
    end = NULL,
    rpm_height = 12,
    bs_height = 0.8,
    annotation_height = 10,
    col_width = 6,
    bs_width = 0
){
  
  if(is.null(start)){
    plot_iclip_coverage(
      gene = gene,
      tx_order = tx_order,
      col_width = col_width,
      bs_width = bs_width
    ) %>%
      patchwork::wrap_plots(
        ncol = 1,
        heights = c(rep(c(rpm_height, bs_height), 4), annotation_height)
      )
  } else {
    plot_iclip_coverage(
      gene = gene,
      tx_order = tx_order,
      col_width = col_width,
      bs_width = bs_width
    ) %>%
      patchwork::wrap_plots(
        ncol = 1,
        heights = c(rep(c(rpm_height, bs_height), 4), annotation_height)
      ) &
      coord_cartesian(xlim = c(start, end))
  }
}



# plot_iclip_coverage_wrap(gene = "pros", tx_order = c("pros-RL", "pros-RK"))


