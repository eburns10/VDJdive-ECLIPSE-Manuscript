library(scRepertoire)
library(Seurat)
library(tidyverse)
library(ECLIPSE)
library(VDJdive)
library(ggpubr)
source("/home/ecb86/Desktop/R_tools/VDJdive_ECLIPSE_paperFunctions.R")
source("/home/ecb86/Desktop/ECLIPSE/revisions/analysis/revisions_functions.R")

braun <- readRDS("/nfs/roberts/pi/pi_db2475/lab_data/sc_assays/braun_cancer_cell_2021/processed/braun_cancer_cell_seurat.rds")

braun@meta.data <- braun@meta.data %>% separate(Sample, into = c("donor", "tissue"), remove = FALSE, sep = "_")
braun@meta.data <- braun@meta.data %>% mutate(tissue = case_when(
  tissue %in% c("T", "T2") ~ "Tumor",
  tissue == "N" ~ "Normal",
  tissue %in% c("M", "V")  ~ "Other"))

cd8 <- subset(braun, subset = ClusterName_TCells %in% c("Treg", "Activated CD4.1", "Activated CD4.2", "CD4 TCM.1", "CD4 TCM.2", "Mito-high T cells", "-1"), invert = TRUE)

dirs <- list.files("/nfs/roberts/pi/pi_db2475/lab_data/sc_assays/braun_cancer_cell_2021/raw/TCR/cellranger_output_indexHoppingRemoved", recursive = TRUE, pattern = "all_contig_annotations\\.csv$", full.names = TRUE)
dirs <- dirs[!str_detect(dirs, "SC_VDJ_ASSEMBLER_CS")]
dirs <- str_replace(dirs, "all.*$", "")

matching <- read.csv("/nfs/roberts/pi/pi_db2475/ecb86/sc_datasets/braun_cancer_cell/braun_matching2.csv")
matching <- matching %>% filter(Sample != "")

cd8@meta.data <- cd8@meta.data %>% left_join(matching[, 2:3], by = "Sample")
cd8@meta.data <- as.data.frame(cd8@meta.data)
rownames(cd8@meta.data) <- cd8$Cell_barcode
cd8$Cell_Barcode <- str_replace(cd8$Cell_barcode, "\\-.*$", "-1")

#ElbowPlot(cd8, ndims = 50)
#cd8 <- RunUMAP(cd8, dims = 1:30)

cd8$ClusterName_TCells[cd8$ClusterName_TCells == "PMCH+ terminally Exhausted CD8.1"] <- "PMCH+ terminally \n Exhausted CD8.1"
cd8$ClusterName_TCells[cd8$ClusterName_TCells == "Terminally exhausted CD8.2"] <- "Terminally \n exhausted CD8.2"
cd8$ClusterName_TCells[cd8$ClusterName_TCells == "Terminally exhausted CD8.3"] <- "Terminally \n exhausted CD8.3"

cd8s <- ECLIPSE(cd8,
                contig_folder_paths = dirs,
                contig_type = "all",
                original_barcode = "Cell_Barcode",
                batch = "index",
                donor = "donor",
                group = "donor")

cd8s <- tcrDoubletDetect(cd8s)


cd8s <- countChains(cd8s)


### Defining the cells that are suitable for this. They must have 2 chains, with 1 TCRa and 1 TCRb. They also cannot be from a clone with 3 chains.
set.seed(16)
cd8s@meta.data <- cd8s@meta.data %>% mutate(suitable = case_when(
  scr_count == 2 & !is.na(CTaa) & !str_detect(scR_CTaa, ";") & !str_detect(CTaa, ":::") ~ TRUE,
  TRUE ~ FALSE))
#cd8s@meta.data <- cd8s@meta.data %>% filter(suitable) %>% distinct(Cell_barcode) %>% mutate(pos = sample.int(n())) %>% right_join(cd8s@meta.data, by = "Cell_barcode")

### Running the pro-processing of ECLIPSE. Pre is just all the cleaning, and then returns the output
threeA <- pre2(contig_folder_paths = dirs, 
               contig_type = "all", 
               seurat_object = cd8, 
               batch = "index", 
               donor = "donor", 
               group = "donor", 
               original_barcode = "Cell_Barcode")

### Finding suitable barcodes
suitable_barcodes <- (cd8s@meta.data %>% filter(suitable) %>% select(Cell_barcode))$Cell_barcode

d1 <- data.frame()
d2 <- data.frame()
d3 <- data.frame()
d4 <- data.frame()


for (i in 1:10) {
  
    message("Starting iteration ", i)
    
    ### Extracting the outputs from pre.
    contigs <- threeA[[1]]
    obj <- threeA[[2]]
    
    ### Giving them a unique key for joining, before removing unsuitable cells.
    contigs <- contigs %>% unite(c(barcode, cdr3), sep = "_", remove = FALSE, col = "key")
    n_distinct(contigs$key) == nrow(contigs)
    contigs <- contigs %>% distinct(key) %>% mutate(pos = sample.int(n())) %>% left_join(contigs, by = "key")
    contigs$suitable <- contigs$barcode %in% suitable_barcodes
    
    ### Considering only contigs that aren't doublets and have the lowest pos variable in their cell which is used for randomization. This ensures 2 contigs aren't removed from the same cell.
    ### Then selecting only 500 contigs to remove and retaining only the key
    contigs2 <- contigs %>% filter(barcode %in% cd8s$Cell_barcode) %>% 
      group_by(barcode) %>% filter(pos == min(pos)) %>% ungroup() %>%
      filter(suitable) %>% arrange(desc(pos))
    contigs2 <- contigs2 %>% filter(row_number() <= 500)
    barcodes <- contigs2$barcode
    contigs2 <- contigs2$key
    
    ### Removing those 500 contigs from the original contigs file, and removing extra columns
    contigs3 <- contigs %>% filter(key %in% contigs2 == FALSE)
    contigs3 <- contigs3 %>% select(-pos, -key)
    
    
    ### Taking the new contigs (reduced by 500) and running them through the second half ECLIPSE pipeline
    pst <- post2(contig_folder_paths = dirs,
                 contig_type = "all",
                 all = contigs3, 
                 seurat_object = obj, 
                 batch = "index", 
                 donor = "donor", 
                 group = "donor", 
                 original_barcode = "Cell_Barcode")
    
    ### Extracting the new clone call and whether it was called by VDJdive or not. Then joining to the regular seurat object for comparison with/without the 500 contigs.
    pst2 <- pst %>% select(Cell_barcode, "new_CTaa" = CTaa, "new_scrCTaa" = scR_CTaa, "newCloneCall" = cloneCallSource)
    join <- cd8s@meta.data %>% left_join(pst2, by = "Cell_barcode")
    
    ### Defining whether each cell was changed and the condition it was in
    altered_barcodes <- str_replace(contigs2, "_.*$", "")
    join <- join %>% mutate(altered = case_when(
      Cell_barcode %in% altered_barcodes ~ TRUE,
      !Cell_barcode %in% altered_barcodes ~ FALSE),
      steady = case_when(
        CTaa == new_CTaa ~ "Same Call",
        CTaa != new_CTaa & newCloneCall == "scRepertoire" ~ "Call Now\nAmbiguous",
        CTaa != new_CTaa & newCloneCall == "VDJdive" ~ "New Call"))
    
    ### Calculating the percentages of each for each group
    join2 <- join %>% filter(!is.na(steady) & !is.na(altered)) %>% group_by(altered) %>% mutate(altered_count = n()) %>% ungroup() %>% 
      group_by(altered, steady) %>% summarise(pct = 100 * n() / first(altered_count),
                                              count = n(),
                                              altered_count = first(altered_count))
    
    ### Making graph
    join2$steady <- factor(join2$steady, levels = c("New Call", "Call Now\nAmbiguous", "Same Call"))
    
    join2$iteration <- as.character(i)
    d1 <- rbind(d1, join2)
    
    
    ### Cell can't be from 3 chain clone originally
    gt <- join$CTaa[join$altered]
    after_removal <- join$new_scrCTaa[join$altered]
    after_eclipse <- join$new_CTaa[join$altered]
    
    gt_shannon <- shannon(gt)
    after_removal_shannon <- shannon(after_removal)
    after_eclipse_shannon <- shannon(after_eclipse)
    
    
    gt_shannon
    after_removal_shannon
    after_eclipse_shannon
    
    diversity <- data.frame("gt" = gt_shannon[1, 1], "sim" = after_removal_shannon[1, 1], "sim_e" = after_eclipse_shannon[1, 1], "iteration" = as.character(i))
    d2 <- rbind(d2, diversity)
    
    
    mod_clones <- join %>% filter(CTaa %in% scR_CTaa[join$Cell_barcode %in% barcodes])
    mod_clones <- mod_clones %>% select(scR_CTaa, CTaa, new_scrCTaa, new_CTaa, CTaa, clonalProportion, clonalFrequency, cloneSize, newCloneCall, altered, steady)
    mod_clones <- mod_clones %>% mutate(cloneSize2 = factor(case_when(
      clonalProportion >= 0.05 ~ "Large",
      clonalProportion < 0.05 & clonalProportion >= 0.005 ~ "Medium",
      clonalProportion < 0.005 ~ "Small")))
    
    
    
    mod_clones2 <- mod_clones %>% filter(altered == TRUE) %>%
      group_by(cloneSize2) %>%
      mutate(count = n()) %>%
      group_by(cloneSize2, steady) %>%
      summarise(pct = 100 * n() / first(count),
                ct = first(count))
    
    mod_clones2$iteration <- as.character(i)
    d3 <- rbind(d3, mod_clones2)
    
    mod_clones3 <- mod_clones %>% filter(altered == TRUE)
    ctaa_counts <-  mod_clones3 %>% group_by(CTaa) %>% summarise(ctaa_count = n())
    new_scr_counts <- mod_clones3 %>% group_by(new_scrCTaa) %>% summarise(new_scr_count = n())
    new_ctaa_counts <- mod_clones3 %>% group_by(new_CTaa) %>% summarise(new_ctaa_count = n())
    
    all_counts <- ctaa_counts %>% left_join(new_scr_counts, by = c("CTaa" = "new_scrCTaa"), keep = FALSE)
    all_counts <- all_counts %>% left_join(new_ctaa_counts, by = c("CTaa" = "new_CTaa"), keep = FALSE)
    all_counts[is.na(all_counts)] <- 0
    
    all_counts2 <- all_counts %>% pivot_longer(cols = c("ctaa_count", "new_scr_count", "new_ctaa_count"), names_to = "class", values_to = "count")
    
    all_counts2$iteration <- as.character(i)
    d4 <- rbind(d4, all_counts2)
}

write_csv(d1, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/simulations/remove/d1.csv")
write_csv(d2, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/simulations/remove/d2.csv")
write_csv(d3, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/simulations/remove/d3.csv")
write_csv(d4, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/simulations/remove/d4.csv")
