# This script is to combine all maize TARGET results.

# Author: Ji Huang

# 0. Prep -----------------------------------------------------------------

library(tidyverse)
library(here)
library(fs)
library(scales)
library(cowplot)

FDR_cutoff <- 0.05
FC_cutoff <- 1 

# 1. Load Plasmid Key Table -----------------------------------------------

key_table <- read_tsv(here("data", "TARGET_results",
                           "maize_TARGET_table.tsv")) %>% 
    janitor::clean_names() %>% 
    select(grassius_id, v4id) %>% 
    distinct()


# 2. Process all No Nitrogen results --------------------------------------

## Get all the -N result files.
noN_files <- dir_ls(here("data", "TARGET_results",
                         "TARGET.N-noN-separately"),
       regexp = "ZmTARGET.*-gene-.", recurse = TRUE)

## Read all files into a big row-bind dataframe.
noN_raw_data <- noN_files %>% 
    map_dfr(read_tsv, skip = 1, col_names = F, .id = "source")

## Clean up the name.
noN_data <- noN_raw_data %>% 
    mutate(source = str_extract(.$source, "ZmTARGET.*edgeR") %>% 
               str_remove("ZmTARGET_") %>% str_remove(".edgeR")
           ) %>% 
    rename(target = X1, logFC = X2, logCPM = X3, 
           Fstat = X4, PValue = X5, FDR = X6) %>% 
    separate(col = source, into = c("Temp1", "time"), sep = "_") %>% 
    separate(col = Temp1, into = c("plasmid", "name"), sep = "-")

noN_data <- noN_data %>% 
    left_join(key_table, by = c("plasmid" = "grassius_id")) %>% 
    select(plasmid, name, regulator = v4id, everything()) %>% 
    mutate(regulator = case_when(plasmid == "UT5628G2mutant" ~ "Zm00001d039260",
                                 plasmid == "KN1" ~ "Zm00001d033859",
                                 TRUE ~ regulator)) %>% 
    mutate(name = if_else(plasmid == "KN1", "KN1", name))

## Filter FDR and logFC
noN_data <- noN_data %>% 
    filter(FDR < FDR_cutoff) %>% 
    filter(logFC > log2(FC_cutoff) | logFC < -log2(FC_cutoff))

noN_summary <- noN_data %>% count(plasmid, sort = TRUE)

# 3. Process all with Nitrogen Results ------------------------------------

## Same code as the -N before.
N_files <- dir_ls(here("data", "TARGET_results",
                         "TARGET.N-noN-separately"),
                    regexp = "ZmTARGET.*-geneN-.", recurse = TRUE)

N_raw_data <- N_files %>% 
    map_dfr(read_tsv, skip = 1, col_names = F, .id = "source")

N_data <- N_raw_data %>% 
    mutate(source = str_extract(.$source, "ZmTARGET.*edgeR") %>% 
               str_remove("ZmTARGET_") %>% str_remove(".edgeR")
    ) %>% 
    rename(target = X1, logFC = X2, logCPM = X3, 
           Fstat = X4, PValue = X5, FDR = X6) %>% 
    separate(col = source, into = c("Temp1", "time"), sep = "_") %>% 
    separate(col = Temp1, into = c("plasmid", "name"), sep = "-")

N_data <- N_data %>% 
    left_join(key_table, by = c("plasmid" = "grassius_id")) %>% 
    select(plasmid, name, regulator = v4id, everything()) %>% 
    mutate(regulator = case_when(plasmid == "UT5628G2mutant" ~ "Zm00001d039260",
                                 plasmid == "KN1" ~ "Zm00001d033859",
                                 TRUE ~ regulator)) %>% 
    mutate(name = if_else(plasmid == "KN1", "KN1", name))

N_data <- N_data %>% 
    filter(FDR < FDR_cutoff) %>% 
    filter(logFC > log2(FC_cutoff) | logFC < -log2(FC_cutoff))

N_summary <- N_data %>% count(plasmid, sort = TRUE)

# 4. Combine -N and +N ----------------------------------------------------

## Combined two tables and keep the edge that has the lowest FDR.
DEall <- bind_rows(noN_data, N_data) %>% 
    mutate(edge = paste0(plasmid, "_", target)) %>% 
    group_by(edge) %>% top_n(n = 1, wt = -FDR) %>% 
    distinct(edge, .keep_all = T) %>% 
    ungroup() %>% 
    select(-edge)

combine_summary <- DEall %>% count(plasmid, sort = TRUE)

# 5. Save the cleaned tables ----------------------------------------------

## Save the combined DE genes table (lowest FDR gene kept).
write_tsv(DEall, file = here("result", "zma_TARGET", "CYC", 
                              "JH_combined_data", 
                             paste0("combined_clean_table_FDR_", 
                                    FDR_cutoff, "_FC_", 
                                    FC_cutoff, ".tsv.gz"))
          )