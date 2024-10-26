# This script calculate AUPR and plot PR curves for DFG, GENIE3, dynGENIE3 and PCC networks.

# Conclusion: DFG is way better than GENIE3, dynGENIE3 and PCC.

# Author: Ji Huang
# Date: 2023-04-07

# 0. Prep -----------------------------------------------------------------
## Prep libraries
library(here)
library(tidyverse)
library(precrec)
source(here("src", "functions","ggplot2_image_settings.R"))
source(here("src", "functions", "aupr_functions.R"))

output_dir = here("result", "GENIE3_PCC", "root")

rep_times <- 1000

## Load network.
load_ntwk <- function(ntwk_file) {
    
    library("tools")
    
    # Load expressed genes in TARGET experiment. Provided by CYC.
    # CYC filtered edges with only expressed genes and get higher AUPR.
    express_gene <- read_tsv(here("data", "TARGET_results", "CYC",
                                  "TARGET.N-noN-separately", "background",
                                  "Expressed-gene-20TARGETs.union-202111-v2.txt"),
                             col_names = "gene", skip = 1)
    
    if (file_ext(ntwk_file) == "gz") {
        ntwk <- read_tsv(ntwk_file, skip=1,
                         col_names = c("regulator", "target", "weight"))
        
    } else if (file_ext(ntwk_file) == "RDS") {
        ntwk <- readRDS(file = ntwk_file)
        colnames(ntwk) <- c("regulator", "target", "weight")
        ntwk <- ntwk %>% 
            mutate(regulator = as.character(regulator),
                   target = as.character(target),
                   weight = as.numeric(weight))
    }
    
    ntwk <- ntwk %>% 
        filter(regulator != target) %>%  # Discard TF to itself.
        filter(regulator %in% express_gene$gene) %>% 
        filter(target %in% express_gene$gene)
    
    return(ntwk)
}


# 1. Load data ------------------------------------------------------------

dfg <- load_ntwk(here("result", "DFG", "root_DFG",
                      "best_zma_MTSroot_20210501_gamma1_lambda_w0.005_tau6.sif.gz"))

g3 <- load_ntwk(here("result", "GENIE3_PCC", "root",
                     "mts_root_genie3_196tf_all.RDS"))

d3 <- load_ntwk(here("result", "GENIE3_PCC", "root",
                     "mts_root_dyngenie3_196tf_All_edge.RDS"))

pcc <- load_ntwk(here("result", "GENIE3_PCC", "root",
                      "mts_root_pcc_196tf_all.RDS"))

# Load TARGET result for validation 
ptarget <- load_zma_target()

## No JH's TF are in the NxT genes.

all_regulator <- unique(dfg$regulator)[unique(dfg$regulator) %in%
                                           unique(ptarget$regulator)]

length(all_regulator) # 11 TFs in the root NxT genes have TARGET

ptarget <- ptarget %>% 
    filter(regulator %in% all_regulator)

filter_edge <- ptarget$edge

# 2. Confirm edges --------------------------------------------------------

dfg_6 <- confirm_edge(dfg)

g3_6 <- confirm_edge(g3) 

d3_6 <- confirm_edge(d3) 

pcc_6 <- confirm_edge(pcc)


# 3. Calculate AUPR -------------------------------------------------------

dfg_val_tf <- dfg %>% filter(regulator %in% all_regulator)

test1 <- map_dbl(1:rep_times, confirm_edge_random, network = dfg_val_tf)

temp1 <- mutate(dfg_val_tf, edge = paste0(regulator, "-", target)) 
set.seed(which.max(test1))
temp2 <- temp1 %>% mutate(new_weight = sample(weight)) %>% 
    arrange(desc(weight)) %>% 
    filter(regulator %in% all_regulator) %>% 
    mutate(label = if_else(edge %in% filter_edge, "1", "0"))
random_max_weight <- temp2$new_weight
random_max_label <- temp2$label


temp1 <- mutate(dfg_val_tf, edge = paste0(regulator, "-", target)) 
set.seed(which.min(test1))
temp2 <- temp1 %>% mutate(new_weight = sample(weight)) %>% 
    arrange(desc(weight)) %>% 
    filter(regulator %in% all_regulator) %>% 
    mutate(label = if_else(edge %in% filter_edge, "1", "0"))
random_min_weight <- temp2$new_weight
random_min_label <- temp2$label

score_list <- list(dfg_6$weight, g3_6$weight, 
                   d3_6$weight, pcc_6$weight, 
                   random_max_weight, random_min_weight)
label_list <- list(dfg_6$label, g3_6$label, 
                   d3_6$label, pcc_6$label, 
                   random_max_label, random_min_label)

# Calculate the model for three networks.
mmcurves <- evalmod(scores = score_list, labels = label_list, 
                    modnames = c("DFG", "GENIE3", 
                                 "dynGENIE3", "PCC", "random_max", "random_min"),
                    dsids = c(1,2,3,4,5,6))


# Display AUPR values.
auc(mmcurves) %>% 
    filter(curvetypes == "PRC")

# modnames dsids curvetypes      aucs
# 1        DFG     1        PRC 0.3507644
# 2     GENIE3     2        PRC 0.2937120
# 3  dynGENIE3     3        PRC 0.2808599
# 4        PCC     4        PRC 0.2866665
# 5 random_max     5        PRC 0.2783315
# 6 random_min     6        PRC 0.2597041

# 4. Plot and Save PR curves ----------------------------------------------

sos_df <- fortify(mmcurves)
sos_df <- sos_df %>% filter(curvetype == "PRC")

p1 <- ggplot()+
    geom_line(data= subset(fortify(sos_df), curvetype == "PRC"), 
              aes(x = x, y = y, color = modname)) +
    ylim(c(0,1)) +
    ylab("Precision") + xlab("Recall") +
    scale_color_manual(values=c('#e41a1c', "#377eb8", 
                                "#984ea3","#4daf4a", "grey56", "grey80"))+
    labs(color = "Model") + 
    theme(legend.position = c(0.7, 0.7)) +
    coord_equal()

p1 <- p1 + 
    ggtitle("PR curves from four methods")


cutoff <- 0.45
p1_cut <- p1 + geom_hline(yintercept = cutoff, linetype = "dashed" )

pdf(file = here(output_dir,
                "edge_2294_4methods_PRcurve_11TFs.pdf"),
    width = 4, height = 4)
p1
p1_cut
dev.off()


