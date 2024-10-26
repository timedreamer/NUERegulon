# This script is to calculate AUPR and plot PR curve for DFG network
# for MTS dataset (edge) gene list. 

# Date: 2022-01-25.

# 0. Prep -----------------------------------------------------------------

## Prep libraries
library(here)
library(tidyverse)
library(precrec)
source(here("src", "functions", "ggplot2_image_settings.R"))

source(here("src", "functions", "aupr_functions.R"))

output_dir = here("result", "DFG")

# 1. Load DFG networks -------------------------------------

# Load expressed genes in TARGET experiment.
express_gene <- read_tsv(here("Expressed-gene-20TARGETs.union-202111-v2.txt"),
                         col_names = "gene", skip = 1)

## Load the DFG network.DFG was calculated on HPC. This is the edge gene list.
dfg <- read_tsv(here("result", "DFG", "best_models_jh",
                     "edges_best_zma_gamma1_lambda_w0.05_tau6_edge_list.tsv.gz"), skip=1,
                col_names = c("regulator", "target", "weight")) %>%
    filter(regulator != target) %>%  # Discard TF to itself.
    filter(target %in% express_gene$gene)

# 2. Load TARGET result for validation ------------------------------------

ptarget <- read_tsv(here("result", "zma_TARGET",
                         "23TFs_TARGET_used4paper.tsv"))

all_regulator <- unique(ptarget$regulator)
filter_edge <- ptarget$edge

# 3. Confirm whether edges are validated. ---------------------------------

dfg_6 <- confirm_edge(dfg)

dfg_6 %>% separate(edge, into = c("regulator", "target"), 
                   sep = "-") %>% 
    mutate(label = as.integer(label)) %>% group_by(regulator) %>% 
    summarise(total_positive = sum(label)) %>% arrange(desc(total_positive))

rand_aucs <- c()

dfg_val_tf <- dfg %>% filter(regulator %in% all_regulator)

# Repeat the random process for 1000 times.
rand_ntwk <- map_dbl(1:1000, confirm_edge_random, network = dfg_val_tf)

## Get the max random AUPR
temp1 <- mutate(dfg_val_tf, edge = paste0(regulator, "-", target))
set.seed(which.max(rand_ntwk))
temp2 <- temp1 %>% mutate(new_weight = sample(weight)) %>%
    arrange(desc(weight)) %>%
    filter(regulator %in% all_regulator) %>%
    mutate(label = if_else(edge %in% filter_edge, "1", "0"))
random_max_weight <- temp2$new_weight
random_max_label <- temp2$label

# Get the min random AUPR.
temp1 <- mutate(dfg_val_tf, edge = paste0(regulator, "-", target))
set.seed(which.min(rand_ntwk))
temp2 <- temp1 %>% mutate(new_weight = sample(weight)) %>%
    arrange(desc(weight)) %>%
    filter(regulator %in% all_regulator) %>%
    mutate(label = if_else(edge %in% filter_edge, "1", "0"))
random_min_weight <- temp2$new_weight
random_min_label <- temp2$label

## Combine scores and labels for three methods.

lengh_random <- length(dfg_6$weight)
length(random_max_weight[1:lengh_random])

score_list <- list(dfg_6$weight, 
                   random_max_weight[1:lengh_random], 
                   random_min_weight[1:lengh_random])
label_list <- list(dfg_6$label,  
                   random_max_label[1:lengh_random], 
                   random_min_label[1:lengh_random])
# Calculate the model for three networks.
mmcurves <- evalmod(scores = score_list, labels = label_list, 
                    modnames = c("DFG","random_max", "random_min"),
                    dsids = c(1,2,3))

# Display AUPR values.
auc(mmcurves)
mean(rand_ntwk)

# 4. Plotting total PR curve. ----------------------------------------

## Get a pvalue
(1000 - sum(auc(mmcurves)[2, 4] > rand_ntwk)) / 1000

p1 <- plot_PRcurve(input_data = mmcurves)

cutoff <- 0.33
p1_cut <- p1 + geom_hline(yintercept = cutoff, linetype = "dashed" )
p1_cut

pdf(file = here(output_dir,
                "edge_2732_DFG_PRcurve_21TFs.pdf"),
    width = 4, height = 4)
p1
p1_cut
dev.off()

# 5. Choose the Precision cutoff value. -----------------------------------

mmcurves <- evalmod(scores = score_list, labels = label_list, 
                    modnames = c("DFG","random_max", "random_min"),
                    dsids = c(1,2,3), mode = "basic")
mmcurves_df <- as_tibble(mmcurves)

mmcurves_df_score <- mmcurves_df %>% 
    filter(type == "score" & modname == "DFG")

mmcurves_df_precision <- mmcurves_df %>% 
    filter(type == "precision" & modname == "DFG")

mmcurves_combine <- left_join(mmcurves_df_score, mmcurves_df_precision, 
                              by = c("x", "modname"))
View(mmcurves_combine)

# 6. Save Pruned network --------------------------------------------------

## I chose AUPR=0.33 as cutoff. The `weight` then is 1.267024. This leaves 15762 edges.
dfg_pruned <- dfg %>% filter(weight > 1.267024)

write_tsv(dfg_pruned, file = here(output_dir,
                                  "zma_edge2732_pruned_ntwk_21TF_weight1.267_Precision0.33.tsv"))


# 7. TF and Genes in Pruned Network ---------------------------------------

dfg_pruned <- read_tsv(here(output_dir,
                            "zma_edge2732_pruned_ntwk_21TF_weight1.267_Precision0.33.tsv"))

zma_TFlist <- read_tsv(here("data", "TF_lists",
                            "ptfdb-grassius_maizeTF_list_orgainzed_v4.txt"))

# TF and Genes in dfg_pruned_shoot
length(unique(dfg_pruned$regulator))
length(unique(dfg_pruned$target))
