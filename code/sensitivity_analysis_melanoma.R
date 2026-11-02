library(scRepertoire)
library(Seurat)
library(tidyverse)
library(ECLIPSE)
library(VDJdive)
source("/home/ecb86/Desktop/R_tools/VDJdive_ECLIPSE_paperFunctions.R")
source("/home/ecb86/Desktop/ECLIPSE/revisions/analysis/revisions_functions.R")

cd8so <- readRDS("/nfs/roberts/pi/pi_db2475/lab_data/sc_assays/oliveira_melanoma/processed/tils.CD8.R0.6.harmonized.20200422.rds")
cd8so <- UpdateSeuratObject(cd8so)
DimPlot(cd8so, label = TRUE)

cd8so$seurat_barcode <- rownames(cd8so@meta.data)
# cd8so$patient holds the patient info


dirs <- list.dirs("/nfs/roberts/pi/pi_db2475/lab_data/sc_assays/oliveira_melanoma/raw/TCR")
dirs <- dirs[-1]

matching <- data.frame("index" = 1:16, "dirs" = str_extract(dirs, "Pt.*$"))

### This shows matching --> table(cd8s$sample, cd8s$patient)
matching$sample <- c("TIL_pre_CD45pos_CD3pos-2", # Pt11 (PtC)_pre_CD45pos_CD3pos_2
                     "TIL_pre_CD45pos_CD3pos-3", # Pt11 (PtC)_pre_CD45pos_CD3pos_3
                     "TIL_pre_CD45pos_CD3pos-1", # Pt11 (PtC)_til_pre_CD45pos_CD3pos_1
                     "TIL_rel-1", # Pt11 (PtC)_til_rel_1
                     "TIL_rel-2", # Pt11 (PtC)_til_rel_2
                     "TIL_CD45pos_CD3pos-1.rds", # Pt15 (PtD)TIL_CD45pos_CD3pos-1
                     "TIL_CD45pos_CD3pos-2.rds", # Pt15 (PtD)TIL_CD45pos_CD3pos-2
                     "TIL_CD45pos_CD3pos-3.rds", # Pt15 (PtD)TIL_CD45pos_CD3pos-3
                     "TIL_CD45pos_CD3pos-4.rds", # Pt15 (PtD)TIL_CD45pos_CD3pos-4
                     "TIL_CD45pos_CD3pos-1", # Pt2 (PtA)TIL_CD45pos_CD3pos-1
                     "TIL_CD45pos_CD3pos-2", # Pt2 (PtA)TIL_CD45pos_CD3pos-2
                     "TIL_CD45pos_CD3pos-3", # Pt2 (PtA)TIL_CD45pos_CD3pos-3
                     "TIL_sel_1", # Pt6 (PtB)TIL_sel_1
                     "TIL_sel_2", # Pt6 (PtB)TIL_sel_2
                     "TIL_sel_3", # Pt6 (PtB)TIL_sel_3
                     "TIL_sel_4") # Pt6 (PtB)TIL_sel_4

cd8so@meta.data <- cd8so@meta.data %>% left_join(matching, by = "sample")

cd8so@meta.data <- as.data.frame(cd8so@meta.data)
rownames(cd8so@meta.data) <- cd8so$seurat_barcode

cd8so$raw_barcode <- str_replace(cd8so$seurat_barcode, "^.*_", "")
cd8so$raw_barcode <- str_replace(cd8so$raw_barcode, "-.*$", "-1")

cd8so <- subset(cd8so, patient != "p6")
cd8s <- ECLIPSE(contig_folder_paths = dirs, contig_type = "all", seurat_object = cd8so, batch = "index", donor = "patient", group = "patient", original_barcode = "raw_barcode", verbose = FALSE)


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
               seurat_object = cd8so, 
               batch = "index", 
               donor = "patient", 
               group = "patient", 
               original_barcode = "raw_barcode",
               verbose = FALSE)

### Finding suitable barcodes
suitable_barcodes <- (cd8s@meta.data %>% filter(suitable) %>% select(seurat_barcode))$seurat_barcode

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
               original_barcode = "raw_barcode",
               batch = "index",
               donor = "patient",
               group = "patient",
               all = threeA[[1]],
               thresholds = thresholds_list[[i]],
               verbose = FALSE)
  
  out2 <- out %>% group_by(patient) %>%
    mutate(count = n()) %>%
    group_by(patient, cloneCallSource) %>%
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
    contigs2 <- contigs %>% filter(barcode %in% cd8s$seurat_barcode) %>% 
                            group_by(barcode) %>% 
                            filter(pos == min(pos)) %>% 
                            ungroup() %>%
                            filter(suitable) %>% 
                            arrange(desc(pos))
    contigs2 <- contigs2 %>% filter(row_number() <= 250)
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
                 donor = "patient", 
                 group = "patient", 
                 original_barcode = "raw_barcode",
                 thresholds = thresholds_list[[i]],
                 verbose = FALSE)
    
    ### Extracting the new clone call and whether it was called by VDJdive or not. Then joining to the regular seurat object for comparison with/without the 500 contigs.
    pst2 <- pst %>% select(seurat_barcode, "new_CTaa" = CTaa, "new_scrCTaa" = scR_CTaa, "newCloneCall" = cloneCallSource)
    join <- cd8s@meta.data %>% left_join(pst2, by = "seurat_barcode")
    
    ### Defining whether each cell was changed and the condition it was in
    altered_barcodes <- str_extract(contigs2, "^.*_.*_")
    altered_barcodes <- str_replace(altered_barcodes, ".$", "")
    
    join <- join %>% mutate(altered = case_when(
                                                seurat_barcode %in% altered_barcodes ~ TRUE,
                                                !seurat_barcode %in% altered_barcodes ~ FALSE),
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




write_csv(meta, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/sensitivity_analysis/meta_mel.csv")
write_csv(coverage, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/sensitivity_analysis/coverage_mel.csv")
write_csv(accuracy, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/sensitivity_analysis/accuracy_mel.csv")




