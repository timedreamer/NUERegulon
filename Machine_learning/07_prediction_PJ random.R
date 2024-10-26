# This script builds Random Forest models to predict phenotypes using gene expression data from the Li et al. (2023) study. 
# It repeats the RF model 100 times with random 24 genes from the 353 conserved genes. 


# Outputs:
# 1. PJ_seedling_RF_PCC_24random353genes_20240529.tsv
# 2. PJ_seedling_RF_PCC_randomCompare_20240529.pdf

# Author: Ji Huang
# Date: 2024-05-29

# 0. Preparation -----------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)
library(randomForest)
library(caret)
library(broom)
library(cowplot)
library(ggbeeswarm)

out_dir <- here("data", "external_data", "Li_2023_PJ")

# Function to run Random Forest model and return pheno_corr
run_rf_model <- function(input_data, yresponse) {
  set.seed(1329)
  trainIndex <- createDataPartition(input_data[[yresponse]], p = 0.8, list = FALSE)
  trainData <- input_data[trainIndex, ]
  testData <- input_data[-trainIndex, ]
  
  trainControl <- trainControl(method = "cv", number = 5)
  
  # Use only columns that start with 'Zm'
  predictors <- grep("Zm", names(trainData), value = TRUE)
  
  rfModel <- train(reformulate(predictors, response = yresponse), 
                   data = trainData, method = "rf", 
                   trControl = trainControl,
                   tuneLength = 5,
                   metric = "RMSE")
  
  predictions <- predict(rfModel, testData)
  
  pheno_corr <- cor.test(testData[[yresponse]], predictions, method = "pearson") |> 
    broom::tidy() |> mutate(phenotype = yresponse)
  
  return(pheno_corr)
}

# These are the 353 conserved genes. Choose 24 randomly.
both_nxt_shoot <- readr::read_tsv(here("result", "network_comparison",
                                       "orthofinder2_group", "shoot",
                                       "both_nxt_shoot.tsv"))
gene_list <- unique(both_nxt_shoot$zma_gene)

# 1. Load phenotype -----------------------------------------------------------

sample_info <- read_excel(here(out_dir, "tpj16260-sup-0005-tables1-s11.xls"),
                          sheet = "Table S2", skip = 1) |>
  dplyr::rename(genotype = "Inbreds")  |> 
  mutate(Time = case_when(Time == "CN_1960&70s" ~ "CN6070",
                          Time == "CN_1980&90s" ~ "CN8090",
                          Time == "CN_2000&10s" ~ "CN0010",
                          TRUE ~ Time))

phenotype  <- read_excel(here(out_dir, "tpj16260-sup-0005-tables1-s11.xls"), skip = 1, sheet = "Table S1", na = "NA")  |> 
  left_join(sample_info, by = c("Inbred_EN" = "genotype")) 

column_names_with_na <- colnames(phenotype)[apply(phenotype, 2, anyNA)]
print(column_names_with_na)

# Fill NA values with group mean based on the "Time" column
phenotype2 <- phenotype %>%
  group_by(Time) %>%
  mutate(across(all_of(column_names_with_na), 
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# 2. Load seedling gene expression ---------------------------------------------

seedling_cpm <- read_tsv(here(out_dir, "seedling", 
                              "seedling_log2cpm_matrix.tsv.gz"))

# Repeating the random sampling 100 times
results_list <- vector("list", 100)

for (i in 1:100) {
  set.seed(968 + i)
  rand1 <- sample(gene_list, 24)
  
  myb_cpm  <- seedling_cpm  |> 
    filter(geneid %in% rand1) |> 
    as.data.frame()
  
  myb_cpm2 <- t(myb_cpm)
  colnames(myb_cpm2) <- myb_cpm2[1,]
  myb_cpm2 <- myb_cpm2[-1,]
  rownames(myb_cpm2) <- sub("^[^_]*_", "", rownames(myb_cpm2))
  
  phenotype_seedling  <- phenotype2[phenotype2$Inbred_EN %in% rownames(myb_cpm2), ]
  myb_cpm2 <- myb_cpm2[phenotype_seedling$Inbred_EN, ]
  
  if (sum(phenotype_seedling$Inbred_EN != rownames(myb_cpm2)) == 0) {
    myb_cpm2 <- apply(myb_cpm2, 2, as.numeric)
    print("Inbred name match continue...")
  } else {
    stop("Inbred_EN does not match.")
  }   
  
  input_seedling <- cbind(phenotype_seedling, myb_cpm2)
  
  # 3. Run multiple seedling RF models on all phenotypes -----------------------------------------------------------
  
  phenotypes <- grep("BLUP", colnames(phenotype_seedling), value = TRUE)
  
  results_seedling <- map_df(phenotypes, run_rf_model, input_data = input_seedling)
  results_list[[i]] <- results_seedling
}

# Combine all results into one dataframe
final_results_random <- bind_rows(results_list, .id = "iteration") |> 
    mutate(type = "Random 24 genes") |> 
    mutate(iteration = as.numeric(iteration))
# print(final_results)

write_tsv(final_results_random, here("result", "machine_learning", "extra_dataset", "PJ_seedling_RF_PCC_24random353genes_20240529.tsv"))

# 4. Compare with the real data and plot -----------------------------------------------------------

final_results_random <- read_tsv(here("result", "machine_learning", "extra_dataset", "PJ_seedling_RF_PCC_24random353genes_20240529.tsv"))

results_seedling <- read_tsv(here("result", "machine_learning", "extra_dataset", "PJ_seedling_RF_PCC_24mybgenes_20240528.tsv")) |> 
    mutate(type = "Zma NUENet regulon genes (24)", iteration = 0)

results_seedling <- results_seedling[, colnames(final_results_random)]

final_results_all <- bind_rows(final_results_random, results_seedling) |> 
    mutate(type = factor(type, levels = c("Zma NUENet regulon genes (24)", "Random 24 genes")))

final_results_all <- final_results_all |> 
    mutate(phenotype = str_remove(phenotype, "_BLUP"))

p_seedling_random <- ggplot(final_results_all, aes(x = phenotype, y = estimate, 
    color = type, size = type)) +
    geom_beeswarm() +
    scale_color_manual(values = c("#cd8d00", "grey80")) +
    scale_size_manual(values = c(3, 1)) +
    scale_alpha_manual(values = c(1, 0.5)) +
    labs(y = "PCC") +
    cowplot::theme_minimal_vgrid() +
    theme(legend.position = "bottom")

ggsave2(here("result", "machine_learning", "extra_dataset", "PJ_seedling_RF_PCC_randomCompare_20240529.pdf"), 
    p_seedling_random, width = 5, height = 4)  
