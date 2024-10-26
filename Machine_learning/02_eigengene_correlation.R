# This script calculate the correlation between conserved NxTime genes' eigen genes and maize field phenotype from Cheng, et al 2021 NC.

# Conclusion: the zma-shoot consverved NxTime genes have significant positive correlation with NUE phenotypes.

# Author: Ji Huang
# Date: 2023-10-09

# Output:
# eigen_corr_zma_shootConserved_NxTg.pdf

# 0. Prep -----------------------------------------------------------------

library(caret)
library(tidyverse)
library(here)
library(cowplot)
library(gridExtra)

source(here("src", "functions", "machine_learning.R"))

# Function to calculate eigen genes for a list of genes.
calc_eigen_zma_leaf <- function(data, gene_list) {
    
    # norm_zma_leaf
    norm_zma_leaf_filter <- data %>% 
        filter(geneid %in% gene_list)
    
    train_data <- prepare_data(input_data = norm_zma_leaf,
                               input_gene_list = gene_list,
                               lowN_keyword = "low_n")
    
    temp1 <- train_data %>% 
        select(-class) %>% 
        rownames_to_column(var = "sample_name") %>% 
        mutate(sample_name = str_remove(sample_name, "_n_rep._.*"))
    
    temp1 <- temp1 %>% 
        pivot_longer(cols = !c("sample_name"), 
                     names_to = "geneid", values_to = "gene_expr") %>% 
        group_by(sample_name, geneid) %>% 
        summarise(mean_expr = mean(gene_expr)) %>% 
        pivot_wider(names_from = geneid, values_from = mean_expr) %>% 
        mutate(sample_name = str_replace(sample_name, "x_", "x")) %>% 
        separate(sample_name, sep = "_", into = c("geno","N")) %>% 
        ungroup()
    
    temp1_w_pheno <- mean_pheno %>% 
        left_join(temp1, by = c("geno", "N")) %>% 
        drop_na()
    
    temp2_expr <- temp1_w_pheno %>% 
        mutate(rname = paste0(geno, "_", N)) %>% 
        as.data.frame() %>% 
        column_to_rownames(var = "rname") %>% 
        select(starts_with("Zm"))
    
    temp2_pheno <- temp1_w_pheno[,1:13]
    
    pca_result <- prcomp(temp2_expr, center = F, scale. = F)
    PC1_values <- pca_result$x[, 1]
    
    temp2_pheno <- temp2_pheno %>% 
        bind_cols(PC1 = PC1_values)
    
    return(temp2_pheno)
    
}

# Function to create a scatterplot with correlation using cor.test
create_scatterplot_with_corr <- function(data, x_value, y_value) {
    
    # Calculate Pearson's correlation coefficient using cor.test
    x <- as.numeric(data %>% pull(x_value))
    y <- as.numeric(data%>% pull(y_value))
    
    cor_result <- cor.test(x, y, method = "p")
    
    pcc_value <- round(cor_result[["estimate"]],3)
    p_value <- cor_result[["p.value"]]
    p_value <- 11*p_value # Befferoni correction
    
    # Create the scatterplot
    p <- ggplot(data, aes(x = .data[[x_value]], 
                          y = .data[[y_value]],
                          color = N)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "grey60") +
        theme_cowplot() +
        labs(caption = paste("PCC =", round(pcc_value, 3),
                             "pval = ", format(p_value, scientific = TRUE, trim = TRUE, digit = 3)))
    
    return(p)
}

# Function to arrange and print the plots in multiple pages
arrange_and_print_plots <- function(plot_list, plots_per_page = 6) {
    num_plots <- length(plot_list)
    num_pages <- ceiling(num_plots / plots_per_page)
    
    for (page in 1:num_pages) {
        start_idx <- (page - 1) * plots_per_page + 1
        end_idx <- min(page * plots_per_page, num_plots)
        
        plots_to_print <- plot_list[start_idx:end_idx]
        
        # Arrange plots in a grid with 3 rows and 3 columns
        arranged_plots <- grid.arrange(grobs = plots_to_print, ncol = 3, nrow = 2)
        
        # Print the arranged plots
        print(arranged_plots)
    }
}

# 1. Load data -----------------------------------------------------------

## 1.1 load conserved gene list.
load(here("result", "DEgenes", "N_response", "nxt_gene_list.RData"))

## 1.2 Load field leaf expression data.
norm_zma_leaf <- read_tsv(here("data", "external_data", "Cheng_NC", 
                               "zma_field_norm_expr.tsv.gz")) %>% 
    janitor::clean_names() %>% 
    as.data.frame()

## 1.3 Load maize field phenotype data.
pheno <- read_tsv(here("data", "external_data", "Cheng_NC",
                       "Trait - ANOVA - Corrplot","maize",
                       "20142016-phenotype.txt"))

# 2. Clean up phenotype and expression matrix -----------------------------

## 2.1 Average and clean phenotype table.
mean_pheno <- pheno %>% 
    group_by(Genotype, NRate) %>% 
    summarise(mSB = mean(StoverBiomass),
              mSN = mean(StoverN),
              mGB = mean(GrainBiomass),
              mGN = mean(GrainN),
              mTB = mean(TOTALBIOMASS),
              mT  = mean(TOTAL),
              mNUT= mean(Nutil),
              mTBY= mean(TotalBiomassYield),
              mGY = mean(GrainYield),
              mSNUE = mean(StoverNUE),
              mTNUE = mean(TotalNUE)) %>% 
    rename(geno = Genotype, N = NRate) %>% 
    ungroup() %>% 
    mutate(geno = tolower(geno),
           N = if_else(N == "H", "high", "low"))


# 3. Calculate and plot zma-shoot conserved genes -------------------------

gene_list <- unique(both_species_nxt_shoot$zma_gene)

leaf_eigen <- calc_eigen_zma_leaf(data = norm_zma_leaf, 
                                  gene_list = gene_list) %>% 
    mutate(N = if_else(N=="high", "highN", "lowN"))

colnames(leaf_eigen)
leaf_eigen <- leaf_eigen |> 
    select(geno, N, mNUT, mSNUE, mTNUE, everything())

all_corrplots <- map(colnames(leaf_eigen[3:13]), 
           .f = create_scatterplot_with_corr, data = leaf_eigen, x_value = "PC1")

# Call the function to print the plots in multiple pages
plot_path <- here("result", "machine_learning", "shoot",
                  "eigen_corr_zma_shootConserved_NxTg.pdf")
pdf(plot_path, width = 10, height = 6)
arrange_and_print_plots(all_corrplots)
dev.off()
