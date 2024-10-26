# Functions to process and plot the MEME-SEA output.

# Author: Ji Huang
# Date: 2023-01-27

# Function1: Load SEA output
load_SEA_output <- function(SEA_output_dir) {
    
    stopifnot(dir.exists(SEA_output_dir))
    
    f_files<- list.files(SEA_output_dir, pattern = "sea.tsv", 
                         recursive = TRUE, full.names = T)
    
    motif_result <- f_files %>%
        purrr::set_names() %>% 
        purrr::map_dfr(., ~read_tsv(., col_names = TRUE, comment = "#"),
                       .id = "fname") %>% 
        dplyr::mutate(JIT_cat = 
                          str_remove(.$fname, pattern = SEA_output_dir) %>% 
                          str_remove("/sea.tsv") %>% 
                          str_remove("/SEA_")
        )
    
    motif_result <- motif_result %>% 
        dplyr::select(JIT_cat, everything(), -fname) %>% 
        mutate(JIT_cat = factor(JIT_cat, levels = c("J5", "J10", "J15", "J20", 
                                                    "J30", "J45", "J60", 
                                                    "J90", "J120"))) %>% 
        arrange(JIT_cat)
    
    return(motif_result)
    
}

# Function2: prepare_motif_data_4heatmap
prepare_motif_data_4heatmap <- function(motif_result) {
    
    pdat1 <- motif_result %>% 
        select(JIT_cat, ID, QVALUE) %>% 
        mutate(LOG10_QVALUE = -log10(QVALUE)) %>%
        select(-QVALUE) %>% 
        pivot_wider(names_from = ID, values_from = LOG10_QVALUE) %>% 
        replace(is.na(.), 0) %>% 
        as.data.frame() %>% 
        t()
    
    colnames(pdat1) <- pdat1["JIT_cat",]
    pdat1 <- pdat1[-1,]
    
    pdat2 <- matrix(as.numeric(pdat1),
                    ncol = ncol(pdat1))
    rownames(pdat2) <- rownames(pdat1)  
    colnames(pdat2) <- colnames(pdat1)
    
    return(pdat2)
}

# Function3: prepare_motif_data_4upset
prepare_motif_data_4upset <- function(motif_result) {
    
    ups_dat <- motif_result %>% 
        select(JIT_cat, ID)
    
    ups_dat1 <- motif_result %>% group_by(JIT_cat) %>% group_split()
    names(ups_dat1) <- ups_dat %>% group_by(JIT_cat) %>% 
        group_keys() %>% pull(JIT_cat)
    
    ups_dat1 <- map(ups_dat1, pull, ID)
    
    return(ups_dat1)
}

# Function4. Plot motif heatmaps.
# Date: 2023-02-07
plot_heatmap_motif <- function(outdir,
                                   motif_db = c("cluster80", "JASPAR656"),
                                   qval_cutoff = 0.05, ENR_cutoff = 2) {
    
    
    stopifnot(dir.exists(outdir))
    motif_db <- match.arg(motif_db)
    
    output_dir = outdir
    
    # Set color.
    color1 <- brewer.pal(n=9, name = "Blues")
    col_fun = colorRamp2(c(0, 3, 30), c(color1[1],color1[5], color1[9]))
    
    # Load sea reports.
    f_files<- list.files(output_dir, pattern = "sea.tsv", 
                         recursive = TRUE, full.names = T)
    motif_cluster <- load_SEA_output(SEA_output_dir = output_dir)
    
    
    
    if (motif_db == "cluster80") {
        # Settings for cluster80
        # Add family name
        motif_cluster <- motif_cluster %>% left_join(., cluster80, by = "ID") %>%
            mutate(ID = paste0(family, "_", cluster))
        show_row_names = TRUE
        plot_width = 5; plot_height = 5
        plot_name <- "SEA_clusterMotif_heatmap.pdf"
        
        motif_cluster_qval <- motif_cluster %>% 
            filter(QVALUE < qval_cutoff)  %>% filter(ENR_RATIO > ENR_cutoff)
        
        pheatmap_data <- prepare_motif_data_4heatmap(motif_result = motif_cluster_qval)
        
        p_motif <- Heatmap(pheatmap_data, col=col_fun,
                           cluster_rows = TRUE, 
                           cluster_columns = FALSE, 
                           show_row_names = show_row_names,
                           name = "-log10(qval)")
        
        p_motif_noRowClustering <- Heatmap(pheatmap_data, col=col_fun,
                                           cluster_rows = FALSE, 
                                           cluster_columns = FALSE, 
                                           show_row_names = show_row_names,
                                           name = "-log10(qval)")
        
        
    } else if (motif_db == "JASPAR656"){
        # Settings for JASPAR656
        motif_cluster <- motif_cluster
        show_row_names = FALSE
        plot_width = 4; plot_height = 5
        plot_name <- "SEA_JASPAMotif_heatmap.pdf"
        
        # cluster80 <- readxl::read_xlsx(here("data", "external_data", 
        #                                     "ClusterDesc_NCor045.xlsx")) %>% 
        #     janitor::clean_names() %>% 
        #     mutate(ID = paste0("cluster_", cluster))
        
        motif_cluster_qval <- motif_cluster %>% 
            filter(QVALUE < qval_cutoff)  %>% filter(ENR_RATIO > ENR_cutoff)
        
        pheatmap_data <- prepare_motif_data_4heatmap(motif_result = motif_cluster_qval)
        
        # Add N-associated motif markds
        Nrespons_motif <- readxl::read_xlsx(here("data", "external_data",
                                                 "JASPAR2022_plantCore.xlsx"),
                                            sheet = "Nresponsive_TF", skip = 1)
        
        mark_data <- tibble(rowid = which(rownames(pheatmap_data) %in% Nrespons_motif$Motif),
                            Motif = rownames(pheatmap_data)[which(rownames(pheatmap_data) %in% Nrespons_motif$Motif)]) %>% 
            left_join(Nrespons_motif, by = "Motif")
        
        motif_mark <- anno_mark(
            at = mark_data$rowid,
            labels = mark_data$Name,
            which = "row", side = "left",
            labels_gp = gpar(fontsize = 10),
            link_width = unit(2, "mm"),
            padding = unit(0.5, "mm")
        )
        
        p_motif <- Heatmap(pheatmap_data, col=col_fun,
                           cluster_rows = TRUE, 
                           cluster_columns = FALSE, show_row_names = show_row_names,
                           name = "-log10(qval)",
                           left_annotation = rowAnnotation(motif = motif_mark))
        
        p_motif_noRowClustering <- Heatmap(pheatmap_data, col=col_fun,
                                           cluster_rows = FALSE, 
                                           cluster_columns = FALSE, show_row_names = FALSE,
                                           name = "-log10(qval)",
                                           left_annotation = rowAnnotation(motif = motif_mark))
    }
    
    # Save plots
    pdf(file = here(output_dir, plot_name),
        width = plot_width, height = plot_height)
    draw(p_motif)
    draw(p_motif_noRowClustering)
    dev.off()
    
    draw(p_motif_noRowClustering)
    
}