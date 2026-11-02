library(scRepertoire)
library(Seurat)
library(tidyverse)
library(ECLIPSE)
library(VDJdive)

### Setting up df
times <- data.frame("user" = NA, "system" = NA, "elapsed" = NA, "Function" = NA, "dataset" = NA, "batches" = NA, "donors" = NA, "groups" = NA, "cells" = NA)



### Cancer Cell
braun <- readRDS("/nfs/roberts/pi/pi_db2475/lab_data/sc_assays/braun_cancer_cell_2021/processed/braun_cancer_cell_seurat.rds")

braun@meta.data <- braun@meta.data %>% separate(Sample, into = c("donor", "tissue"), remove = FALSE, sep = "_")
braun@meta.data <- braun@meta.data %>% mutate(tissue = case_when(
  tissue %in% c("T", "T2") ~ "Tumor",
  tissue == "N" ~ "Normal",
  tissue %in% c("M", "V")  ~ "Other"))

cd8 <- subset(braun, subset = ClusterName_TCells %in% c("Treg", "Activated CD4.1", "Activated CD4.2", "CD4 TCM.1", "CD4 TCM.2", "Mito-high T cells", "-1"), invert = TRUE)
rm(braun)

dirs <- list.files("/nfs/roberts/pi/pi_db2475/lab_data/sc_assays/braun_cancer_cell_2021/raw/TCR/cellranger_output", recursive = TRUE, pattern = "all_contig_annotations\\.csv$", full.names = TRUE)
dirs <- dirs[!str_detect(dirs, "SC_VDJ_ASSEMBLER_CS")]
dirs <- str_replace(dirs, "all.*$", "")

matching <- read.csv("/nfs/roberts/pi/pi_db2475/ecb86/sc_datasets/braun_cancer_cell/braun_matching2.csv")
matching <- matching %>% filter(Sample != "")

cd8@meta.data <- cd8@meta.data %>% left_join(matching[, 2:3], by = "Sample")
cd8@meta.data <- as.data.frame(cd8@meta.data)
rownames(cd8@meta.data) <- cd8$Cell_barcode
cd8$Cell_Barcode <- str_replace(cd8$Cell_barcode, "\\-.*$", "-1")

cc_eclipse <- system.time(ECLIPSE(cd8,
                                  contig_folder_paths = dirs,
                                  contig_type = "all",
                                  original_barcode = "Cell_Barcode",
                                  batch = "index",
                                  donor = "donor",
                                  group = "donor"))

cd8s <- ECLIPSE(cd8,
                contig_folder_paths = dirs,
                contig_type = "all",
                original_barcode = "Cell_Barcode",
                batch = "index",
                donor = "donor",
                group = "donor")

cc_tcrdd <- system.time(tcrDoubletDetect(cd8s))

### Describing the size of the datatset
batch <- n_distinct(cd8$index)
donor <- n_distinct(cd8$donor)
group <- n_distinct(cd8$donor)
cells <- sum(!is.na(cd8s$CTaa))

### Saving
times[1,] <- c(cc_eclipse[1:3], "ECLIPSE", "cancer_cell", batch, donor, group, cells)
times[2,] <- c(cc_tcrdd[1:3], "tcrDoubletDetect", "cancer_cell", batch, donor, group, cells)

rm(cd8, cd8s)
rm(batch, donor, group, cells)
gc()



### Melanoma CD4
cd4s <- readRDS("/nfs/roberts/pi/pi_db2475/lab_data/sc_assays/oliveira_melanoma/processed/tils.CD4.harmonized.20200514.rds")
cd4s <- UpdateSeuratObject(cd4s)
DimPlot(cd4s, label = TRUE)

cd4s$seurat_barcode <- rownames(cd4s@meta.data)
# cd8s$patient holds the patient info


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

cd4s@meta.data <- cd4s@meta.data %>% left_join(matching, by = "sample")

cd4s@meta.data <- as.data.frame(cd4s@meta.data)
rownames(cd4s@meta.data) <- cd4s$seurat_barcode

cd4s$raw_barcode <- str_replace(cd4s$seurat_barcode, "^.*_", "")
cd4s$raw_barcode <- str_replace(cd4s$raw_barcode, "-.*$", "-1")

melcd4_eclipse <- system.time(ECLIPSE(contig_folder_paths = dirs, contig_type = "all", seurat_object = cd4s, batch = "index", donor = "patient", group = "patient", original_barcode = "raw_barcode"))

ocd4s <- ECLIPSE(contig_folder_paths = dirs, contig_type = "all", seurat_object = cd4s, batch = "index", donor = "patient", group = "patient", original_barcode = "raw_barcode")


melcd4_tcrdd <- system.time(tcrDoubletDetect(ocd4s))

### Describing the size of the datatset
batch <- n_distinct(cd4s$index)
donor <- n_distinct(cd4s$patient)
group <- n_distinct(cd4s$patient)
cells <- sum(!is.na(ocd4s$CTaa))

### Saving
times[3,] <- c(melcd4_eclipse[1:3], "ECLIPSE", "Melanoma CD4", batch, donor, group, cells)
times[4,] <- c(melcd4_tcrdd[1:3], "tcrDoubletDetect", "Melanoma CD4", batch, donor, group, cells)

rm(cd4s, ocd4s)
rm(batch, donor, group, cells)
gc()




### Melanoma CD8
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


melcd8_eclipse <- system.time(ECLIPSE(contig_folder_paths = dirs, contig_type = "all", seurat_object = cd8so, batch = "index", donor = "patient", group = "patient", original_barcode = "raw_barcode"))

ocd8s <- ECLIPSE(contig_folder_paths = dirs, contig_type = "all", seurat_object = cd8so, batch = "index", donor = "patient", group = "patient", original_barcode = "raw_barcode")


melcd8_tcrdd <- system.time(tcrDoubletDetect(ocd8s))

### Describing the size of the datatset
batch <- n_distinct(cd8so$index)
donor <- n_distinct(cd8so$patient)
group <- n_distinct(cd8so$patient)
cells <- sum(!is.na(ocd8s$CTaa))

### Saving
times[5,] <- c(melcd8_eclipse[1:3], "ECLIPSE", "Melanoma CD8", batch, donor, group, cells)
times[6,] <- c(melcd8_tcrdd[1:3], "tcrDoubletDetect", "Melanoma CD8", batch, donor, group, cells)

rm(cd8so, ocd8s)
rm(batch, donor, group, cells)
gc()




### COVID
library(Matrix)
dirs <- list.files("/nfs/roberts/pi/pi_db2475/ecb86/sc_datasets/th17_covid", full.names = TRUE)
samples <- str_extract(dirs, "Bac.*|Sar.*")

ps <- list()
tcr_dirs <- c()

for (i in 1:length(dirs)) {
  
  print(i)
  
  matrix_dir <- list.files(dirs[i], pattern = "matrix", full.names = TRUE)
  matrix <- readMM(matrix_dir[1])
  
  barcode_dir <- list.files(dirs[i], pattern = "barcodes", full.names = TRUE)
  barcodes <- read_tsv(barcode_dir[1], col_names = FALSE)
  barcodes <- barcodes$X1
  
  features_dir <- list.files(dirs[i], pattern = "features", full.names = TRUE)
  features <- read_tsv(features_dir[1], col_names = FALSE)
  features <- features$X1 ### Note: This has both ENSG and names. I'm taking ENSG since I can't have duplicates
  
  tcr_dir <- list.files(dirs[i], pattern = "contig", full.names = TRUE)
  tcr_dirs[i] <- tcr_dir[1]
  
  ps[[i]] <- as.data.frame(matrix)
  colnames(ps[[i]]) <- paste(samples[i], barcodes, sep = "_")
  rownames(ps[[i]]) <- features
  
  if (length(matrix_dir) > 1 | length(barcode_dir) > 1 | length(features_dir) > 1 | length(tcr_dir) > 1) {
    message(paste("Warning, multiple files found"))
  }
  
}

p <- bind_cols(ps)
rm(ps)
gc()

covid <- CreateSeuratObject(p)
rm(p, matrix)
gc()

covid@meta.data <- covid@meta.data %>% separate(orig.ident, into = c("condition", "sample"), sep = 3, remove = FALSE)
covid$tissue <- str_extract(covid$sample, "B|L")
covid$donor <- as.character(str_replace(covid$sample, "B|L", ""))
covid$donor <- paste("D", covid$donor, sep = "")
covid$Barcode <- rownames(covid@meta.data)

matching <- data.frame("index" = c(1:22), "orig.ident" = samples)
covid@meta.data <- covid@meta.data %>% left_join(matching, by = "orig.ident")

covid@meta.data <- as.data.frame(covid@meta.data)
rownames(covid@meta.data) <- covid$Barcode

covid$original_barcode <- str_replace(covid$Barcode, ".*_", "")

covid <- JoinLayers(covid)

covid_eclipse <- system.time(ECLIPSE(contig_type = "manual", contig_file_paths = tcr_dirs, seurat_object = covid, batch = "index", donor = "donor", group = "donor", original_barcode = "original_barcode", format = "none"))

covid2 <- ECLIPSE(contig_type = "manual", contig_file_paths = tcr_dirs, seurat_object = covid, batch = "index", donor = "donor", group = "donor", original_barcode = "original_barcode", format = "none")

covid_tcrdd <- system.time(tcrDoubletDetect(covid2))

### Describing the size of the datatset
batch <- n_distinct(covid$index)
donor <- n_distinct(covid$donor)
group <- n_distinct(covid$donor)
cells <- sum(!is.na(covid2$CTaa))

### Saving
times[7,] <- c(covid_eclipse[1:3], "ECLIPSE", "COVID", batch, donor, group, cells)
times[8,] <- c(covid_tcrdd[1:3], "tcrDoubletDetect", "COVID", batch, donor, group, cells)

rm(covid, covid2)
rm(batch, donor, group, cells)
gc()




### Saving
write_csv(times, "/home/ecb86/Desktop/ECLIPSE/revisions/analysis/run_times.csv")


