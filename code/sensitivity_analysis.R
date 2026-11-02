library(scRepertoire)
library(Seurat)
library(tidyverse)
library(ECLIPSE)
library(VDJdive)
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
                group = "donor",
                verbose = FALSE)

cd8s <- tcrDoubletDetect(cd8s)


#One Time
# Standard ECLIPSE to get suitable barcodes
# First half of removal test.
#
#Loop:
# 1X: % of cells assigned by eclipse in normal run
# 5X times run post() to get accuracy of calling
# In the end this should somewhat get coverage vs. accuracy

cd8s <- countChains(cd8s)


### Defining the cells that are suitable for this. They must have 2 chains, with 1 TCRa and 1 TCRb. They also cannot be from a clone with 3 chains.
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
               original_barcode = "Cell_Barcode",
               verbose = FALSE)

### Finding suitable barcodes
suitable_barcodes <- (cd8s@meta.data %>% filter(suitable) %>% select(Cell_barcode))$Cell_barcode

thresholds_list <- list(c(1, 0.995, 0.99),
                        c(0.99, 0.95, 0.9),
                        c(0.99, 0.9, 0.8),
                        c(0.95, 0.9, 0.8),
                        c(0.95, 0.8, 0.6),
                        c(0.9, 0.8, 0.6),
                        c(0.9, 0.7, 0.5),
                        c(0.8, 0.7, 0.5),
                        c(0.8, 0.6, 0.4),
                        c(0.7, 0.6, 0.4),
                        c(0.7, 0.5, 0.334),
                        c(0.6, 0.5, 0.334),
                        c(0.6, 0.4, 0.25),
                        c(0.5, 0.4, 0.25),
                        c(0.5, 0.334, 0.1))

d1 <- data.frame()
d2 <- data.frame()
d3 <- data.frame()
d4 <- data.frame()

meta <- tibble()
coverage <- tibble()
accuracy <- tibble()

### Outer loop. Each iteration is 1 threshold where we get the coverage, then the inner loop is run 5x times to get a good idea of accuracy
for (i in 1:length(thresholds_list)) {
  
    message("Starting condition ", i)
    set.seed(25)
    
    ### Running ECLIPSE 1 time to get final calls and % of cells assigned by VDJdive
    out <- post2(threeA[[2]],
                 contig_folder_paths = dirs,
                 contig_type = "all",
                 original_barcode = "Cell_Barcode",
                 batch = "index",
                 donor = "donor",
                 group = "donor",
                 all = threeA[[1]],
                 thresholds = thresholds_list[[i]],
                 verbose = FALSE)
    
    out2 <- out %>% group_by(donor) %>%
                    mutate(count = n()) %>%
                    group_by(donor, cloneCallSource) %>%
                    summarise(donorCount = first(count),
                              pct = 100 * n() / first(donorCount),
                              subsetCount = n())
    
    out2$condition <- as.character(i)
    out$condition <- as.character(i)
    
    coverage <- rbind(coverage, out2)
    meta <- rbind(meta, out)
  
    ### Inner loop. This is 10 iterations, so we have a good idea of the accuracy at each threshold
    for (j in 1:10) {
        
        message("Starting iteration ", j, " from condition ", i)
      
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
                     original_barcode = "Cell_Barcode",
                     thresholds = thresholds_list[[i]],
                     verbose = FALSE)
        
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
        
        join2$iteration <- as.character(j)
        join2$condition <- as.character(i)
        
        accuracy <- rbind(accuracy, join2)
    
    }
    
}




write_csv(meta, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/sensitivity_analysis/meta_rcc.csv")
write_csv(coverage, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/sensitivity_analysis/coverage_rcc.csv")
write_csv(accuracy, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/sensitivity_analysis/accuracy_rcc.csv")




