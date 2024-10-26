# This script is to calculate the network using GENIE3, dynGENIE3 and PCC method.
# The DFG outperforms all three methods. The dynGENIE3 ran on HPC.

# Author: Ji Huang

# 0. Prep -----------------------------------------------------------------

# Load libraries.
library(tidyverse)
library(GENIE3)
library(here)

n_cores = 2

# Load expression matrix.
expr_file <- here("result", "DFG", "mts_nxt2732_edge_quant.tsv")

expr <- read.delim(expr_file, sep = "\t", stringsAsFactors = F, header = F)

expr[1:4,1:4]
rownames(expr) <- expr$V1
expr$V1 <- NULL
expr <- as.matrix(expr) # GENIE3 only accept matrix input.

tf_file <- here("result", "DFG", "mts_203_edge_tfs.tsv")
tf_list <- read_tsv(tf_file, col_names = "tf") %>% pull(tf)


# 1. Run GENIE3 -----------------------------------------------------------

# Compute GENIE3 network by providing 145 TFs only.
set.seed(123) # For reproducibility of results
# For weigh matrix, the row is regulator, the columns are targets.
weigh_mat_all <- GENIE3(expr, regulators = tf_list, nCores = n_cores)
link_list_all <- getLinkList(weigh_mat_all)

head(link_list_all)
colnames(link_list_all) <- c("regulator","target","weight")

saveRDS(link_list_all, file = here("result","GENIE3_PCC",
                                   "mts_genie3_203tf_All_edge.RDS"))

# 3. Run PCC --------------------------------------------------------------
# PCC Compute correlation networks.
pcc_matrix <- cor(t(expr), method = "p")

link_list<- GENIE3::getLinkList(pcc_matrix, threshold = -1)

# keep the link list for the 145 TFs only. And get the absolute correlation values.
link_list_tf <- link_list %>% filter(regulatoryGene %in% tf_list) %>% 
    rename(regulator = regulatoryGene, target = targetGene) %>% 
    mutate(weight =abs(weight)) %>% 
    arrange(desc(weight))

# Save the network result.
saveRDS(link_list_tf, file = here("result","GENIE3_PCC", 
                                  "mts_pcc_203tf_All_edge.RDS"))


# 4. Run dynGENIE3 --------------------------------------------------------
expr_file <- here("result", "DFG", "mts_nxt2732_edge_quant.tsv")
expr <- read.delim(expr_file, sep = "\t", stringsAsFactors = F)
expr <- expr[-2422,]
colnames(expr)[1] <- "geneName"
rownames(expr) <- expr$geneName
expr$geneName <- NULL
expr <- as.matrix(expr) # GENIE3 only accept matrix input.

# Prepare the timepoint.
timepoint <- c(0,5,10,15,20,30,45,60,90,120)
timepoint_list <- list(timepoint, timepoint, timepoint,
                       timepoint, timepoint, timepoint)


# Separate three time-series replicate.
ts1 <- expr[,seq(1,30,3)]
ts2 <- expr[,seq(2,30,3)]
ts3 <- expr[,seq(3,30,3)]

ts4 <- expr[,c(1,seq(31,57,3))]
ts5 <- expr[,c(2,seq(32,57,3))]
ts6 <- expr[,c(3,seq(33,57,3))]

ts_data <- list(ts1,ts2,ts3, ts4, ts5, ts6)

save(ts_data, timepoint_list, tf_list,
     file = here("result","GENIE3_PCC", "dynGENIE3_input.RData"))

## Then run dynGENIE3 on HPC as below.
# module load r/intel/4.0.4
# source("dynGENIE3.R")
# load("dynGENIE3_input.RData")
# result <- dynGENIE3(TS.data = ts_data, time.points = timepoint_list, regulators = tf_list, ncores = 8)
# # keep the link list for the 145 TFs only.
# link_list <- get.link.list(result$weight.matrix[tf_list,])
# colnames(link_list) <- c("regulator","target","weight")
# saveRDS(link_list, file = "mts_shoot_dyngenie3_203tf_All_edge.RDS")

