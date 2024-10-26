# In this script, I used the maize atlas data (FPKM values downloaded from qTeller) to plot the Shoot NxTime gene expressoin.

# Output:
# zmaShootNxTAll_leafDev_km7.pdf

# Author: Ji Huang
# Date: 2023-04-17

# 0. Prep -----------------------------------------------------------------

library(here)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)

output_dir <- here("result", "DEgenes", "heatmap", "non_NxTime_data_heatmap")

## A function to read the fpkm files downloaded from qTeller
read_in_fpkm<- function(file){
    fpkm <- read_tsv(file, col_names = c("geneid","fpkm","rm1", "rm2")) %>% 
        select("geneid", "fpkm")
    colnames(fpkm)[2] <- eval(file)
    return(fpkm)
}

## a function to add color to gene-labeling data.
add_gene_color <- function(input_df) {
    
    color_pal <- ggsci::pal_nejm()(5)
    
    input_df <- input_df %>% 
        mutate(gene_color = case_when(pathway_step %in% c("Ammonia Transporters",
                                                          "Nitrate Transporters") ~ color_pal[1],
                                      pathway_step %in% c("Nitrite Reductase",
                                                          "Nitrate Reductase") ~ color_pal[2],
                                      pathway_step %in% c("Glutamate Decarboxylase",
                                                          "Glutamate Dehydrogenase",
                                                          "Glutamate Synthase",
                                                          "Glutamine Synthetase") ~ color_pal[3],
                                      pathway_step %in% c("Asparaginase",
                                                          "Asparagine Synthetase",
                                                          "Aspartate Aminotransferase") ~ color_pal[4],
                                      pathway_step %in% c("Carbamoylphosphate Synthase",
                                                          "Carbonic Anhydrase",
                                                          "Nitrilase") ~ color_pal[5]))
    
    
    input_df <- input_df %>% 
        arrange(factor(gene_color, levels = color_pal))
    
    return(input_df)
}
# 1. Load data ------------------------------------------------------------

## 1.1 Load zma-NxTime genes.
zma_shoot <- read_tsv(here("result", "DEgenes", "edge", 
                           "zma_NxTgenes_df5_edge.tsv"))

## 1.2 Load all the fpkm files in this folder.
list_of_files <- list.files(path = here("data", "external_data", 
                                        "qTeller_v4", "atlas_leaf", "time"), 
                            full.names = TRUE)


raw_fpkm <- map(list_of_files, read_in_fpkm)
raw_fpkm_df<- purrr::reduce(raw_fpkm, inner_join, by = "geneid")

## 1.3 Clean the fpkm df
raw_fpkm_df <- as.data.frame(raw_fpkm_df)
colnames(raw_fpkm_df) <- colnames(raw_fpkm_df) %>% 
    str_remove(".*/time/") %>% 
    str_remove(".fpkm_tracking.gz") %>% 
    str_remove("B73_")

rownames(raw_fpkm_df) <- raw_fpkm_df$geneid
raw_fpkm_df$geneid <- NULL

## keep only zma-shoot-NxTime genes
fpkm_shoot <- raw_fpkm_df[zma_shoot$gene,]

## remove NA and sd=0 row.
fpkm_shoot <- fpkm_shoot[complete.cases(fpkm_shoot),]
fpkm_shoot <- fpkm_shoot[apply(fpkm_shoot, 1, function(x) 
    {all(sd(x, na.rm = TRUE) != 0)}),]

## re-order column from young to adult
correct_col_order <- c("6_DAS_GH_Coleoptile",
                       "V1_4D_PE_Pooled_Leaves", "V1_4D_PE_Stem__SAM",
                       "V3_Topmost_leaf", "V3_Stem_andSAM", 
                       "V5_Bottom_of_transition_leaf", 
                       "V5_Shoot_tip", "V5_Tip_ofStage_2_Leaf",
                       "V7_Bottom_of_transition_leaf", 
                       "V7_Tip_of_transition_leaf",
                       "V9_Eighth_Leaf", "V9_Eleventh_Leaf",
                       "V9_Immature_Leaves", "V9_Thirteenth_Leaf",
                       "VT_Thirteenth_Leaf", "R2_Thirteenth_Leaf")

fpkm_shoot <- fpkm_shoot[correct_col_order]

## row scale the fpkm (z-score)
fpkm_shoot_scale <- fpkm_shoot %>% pheatmap:::scale_rows()

## 1.4 Load the Npathway genes in maize
Npath_genes <- read_tsv(here("result", "misc", 
                             "zma_Npathway_genes_with_names_20230317.tsv")) %>% 
    distinct()

Npath_genes <- Npath_genes %>% 
    filter(!pathway_step %in% c("cyanate degradation",
                                "Glycolate Oxidase"))

Npath_genes_shoot <- Npath_genes %>% 
    filter(sgene %in% zma_shoot$gene) %>% 
    select(pathway_step, sgene, zma_gene_symbol) %>% 
    distinct()

# 2. Set up the matrix and column (time) color ---------------------------------

## Matrix color palette.
rdbu <- rev(brewer.pal(9, "RdBu"))
col_fun <- circlize::colorRamp2(c(-2, 0, 2), 
                                c(rdbu[1], rdbu[5], rdbu[9]))
##  Add time anno
time_color <- rev(mako(n = 16))
names(time_color) <- colnames(fpkm_shoot_scale)

ha = HeatmapAnnotation(DevTime = colnames(fpkm_shoot_scale), 
                       col = list(DevTime = time_color),
                       border = TRUE,
                       annotation_name_side = "left",
                       annotation_legend_param = list(
                           at = colnames(fpkm_shoot_scale)
                       )
)

# 3. Plot all the NxTime genes and label Npathway genes -------------------

## The Npathway color matches the other heatmaps.

## 3.1 Prepare the gene-labeling df
all_Npath_data <- tibble(rowid = which(rownames(fpkm_shoot_scale) %in% Npath_genes_shoot$sgene),
                         geneid = rownames(fpkm_shoot_scale)[which(rownames(fpkm_shoot_scale) %in% Npath_genes_shoot$sgene)]) %>% 
    left_join(Npath_genes_shoot, by = c("geneid"="sgene"))

## add gene color
all_Npath_data <- add_gene_color(all_Npath_data)

## add gene mark
gene_mark <- anno_mark(
    at = all_Npath_data$rowid,
    labels = all_Npath_data$zma_gene_symbol,
    which = "row",
    side = "right",
    labels_gp = grid::gpar(fontface = "italic",
                           fontsize = 10,
                           col = all_Npath_data$gene_color),
    link_width = unit(6, "mm"),
    padding = unit(1, "mm")
)

## Determine the number of kmeans clusters.
factoextra::fviz_nbclust(fpkm_shoot_scale, FUNcluster = kmeans, method = "wss")

set.seed(1655)
htmap_all <- Heatmap(fpkm_shoot_scale, name = "z-scaled \nexpression",
                     km = 7,
                     row_gap = unit(1, "mm"),
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "complete",
                     cluster_rows = TRUE,
                     cluster_columns = FALSE,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     col = col_fun,
                     show_row_dend = FALSE,
                     row_title = NULL,
                     column_title = "ZmaShootNxTimeAll LeafDev",
                     top_annotation = ha,
                     right_annotation = rowAnnotation(gene = gene_mark)
)



pdf(file = here("result", "DEgenes", "heatmap",
                "zmaShootNxTAll_leafDev_km7.pdf"),
    width = 8, height = 8)
draw(htmap_all,
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
