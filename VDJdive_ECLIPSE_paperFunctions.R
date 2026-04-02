

compareClones <- function(object, clusters, cloneThresh, cloneClassThresh){
              
              #### This function sees whether 3 TCR chains or 2 TCR chains are more commonly detected in these clones
              #### It also finds the major and minor bitypic chain of each 3 TCR clone. 
              #### Lastly, it finds the percentage of cells of each chain type in each cluster
              
  
              #### This first part is visualization whether 2 or 3 TCR chains are most commonly detected
              df <- object@meta.data %>% filter(str_detect(CTaa, ":::") & cloneCallSource == "VDJdive" & str_count(CTaa, ":::") == 1)
              
              df2 <- df %>% group_by(CTaa) %>% mutate(group_count = n()) %>% ungroup() %>%
                            group_by(CTaa, scR_CTaa) %>% summarise(count = n(), group_Count = first(group_count))
              
              df2 <- df2 %>% group_by(CTaa) %>% mutate(new_count = case_when(
                                                                             !str_detect(scR_CTaa, "^NA_|_NA$") ~ count,
                                                                             TRUE ~ 0)) %>% ungroup()
              df2 <- df2 %>% separate(CTaa, into = c("final_alpha", "final_beta"), sep = "_", remove = FALSE) %>%
                             separate(scR_CTaa, into = c("scr_alpha", "scr_beta"), sep = "_", remove = FALSE)
              df2 <- df2 %>% mutate(type = case_when(
                                                     str_detect(final_alpha, ":::") ~ "a",
                                                     str_detect(final_beta, ":::") ~ "b",
                                                     TRUE ~ NA))
              a_major <- df2 %>% filter(type == "a" & scr_alpha != "NA") %>% group_by(CTaa, scr_alpha) %>% summarise(count = sum(count))
              a_major2 <- a_major %>% group_by(CTaa) %>% summarise(mx = scr_alpha[which.max(count)])
              a_major2 <- a_major2 %>% mutate(chain = "alpha", type = case_when(
                                                                                str_detect(mx, ";") ~ "Double Chain",
                                                                                !str_detect(mx, ";") ~ "Single Chain"))
              b_major <- df2 %>% filter(type == "b" & scr_beta != "NA") %>% group_by(CTaa, scr_beta) %>% summarise(count = sum(count))
              b_major2 <- b_major %>% group_by(CTaa) %>% summarise(mx = scr_beta[which.max(count)])
              b_major2 <- b_major2 %>% mutate(chain = "beta", type = case_when(
                                                                               str_detect(mx, ";") ~ "Double Chain",
                                                                               !str_detect(mx, ";") ~ "Single Chain"))
              major2 <- rbind(a_major2, b_major2)
              plot1 <- major2 %>% ggplot() + geom_bar(aes(x = chain, fill = type), position = "fill") + theme_classic() + 
                                             labs(fill = "Dominant Type in Clone") + ylab("Proportion") + 
                                             annotate("text", x = 2, y = 1.02, label = paste("n = ", sum(major2$chain == "beta"), sep = "")) + 
                                             annotate("text", x = 1, y = 1.02, label = paste("n = ", sum(major2$chain == "alpha"), sep = ""))
              print(plot1)
              
              
              #### This is about finding the major and minor bitypic chains
              a_major <- a_major %>% mutate(count2 = case_when(
                                                               str_detect(scr_alpha, ";") | !str_detect(CTaa, paste("^", scr_alpha, ":::|:::", scr_alpha, "_", sep = "")) ~ 0,
                                                               !str_detect(scr_alpha, ";") & str_detect(CTaa, paste("^", scr_alpha, ":::|:::", scr_alpha, "_", sep = "")) ~ count))
              a_major <- a_major %>% group_by(CTaa) %>% mutate(mx = scr_alpha[which.max(count2)]) %>% ungroup()
              a_major <- a_major %>% mutate(count3 = case_when(
                                                               mx == scr_alpha ~ 0,
                                                               mx != scr_alpha ~ count2))
              a_major <- a_major %>% group_by(CTaa, mx) %>% summarise(mx2 = scr_alpha[which.max(count3)]) %>% ungroup()
              a_major <- a_major %>% mutate(mx2 = case_when(
                                                            mx2 == mx | str_detect(mx2, ";") ~ "",
                                                            mx2 != mx & !str_detect(mx2, ";") ~ mx2))
                                                                                              
              b_major <- b_major %>% mutate(count2 = case_when(
                                                               str_detect(scr_beta, ";") | !str_detect(CTaa, paste("_", scr_beta, ":::|:::", scr_beta, "$", sep = "")) ~ 0,
                                                               !str_detect(scr_beta, ";") & str_detect(CTaa, paste("_", scr_beta, ":::|:::", scr_beta, "$", sep = "")) ~ count))
              b_major <- b_major %>% group_by(CTaa) %>% mutate(mx = scr_beta[which.max(count2)]) %>% ungroup()
              b_major <- b_major %>% mutate(count3 = case_when(
                                                               mx == scr_beta ~ 0,
                                                               mx != scr_beta ~ count2))
              b_major <- b_major %>% group_by(CTaa, mx) %>% summarise(mx2 = scr_beta[which.max(count3)]) %>% ungroup()
              b_major <- b_major %>% mutate(mx2 = case_when(
                                                            mx2 == mx | str_detect(mx2, ";") ~ "",
                                                            mx2 != mx & !str_detect(mx2, ";") ~ mx2))
              
              combo_major <- rbind(a_major, b_major)
              combo_major <- combo_major %>% filter(!str_detect(mx, ";"))
              
              df3 <- df2 %>% left_join(combo_major, by = "CTaa")
              
              df3 <- df3 %>% mutate(class = case_when(
                                                      type == "a" & final_alpha == str_replace_all(scr_alpha, ";", ":::") ~ "1°\n&\n2°",
                                                      type == "a" & scr_alpha == mx ~ "1°",
                                                      type == "a" & scr_alpha == mx2 ~ "2°",
                                                      type == "a" & scr_alpha == "NA" ~ "NA",
                                                      type == "b" & final_beta == str_replace_all(scr_beta, ";", ":::") ~ "1°\n&\n2°",
                                                      type == "b" & scr_beta == mx ~ "1°",
                                                      type == "b" & scr_beta == mx2 ~ "2°",
                                                      type == "b" & scr_beta == "NA" ~ "NA",
                                                      TRUE ~ "Other"))
              
              
              
              df4 <- df3 %>% select(CTaa, scR_CTaa, type, mx, mx2, class)
              df4 <- df4 %>% rowwise() %>% mutate(code = paste(CTaa, scR_CTaa, sep = "-"))
              df4 <- df4 %>% select(-CTaa, -scR_CTaa)
              
              df <- df %>% rowwise() %>% mutate(code = paste(CTaa, scR_CTaa, sep = "-"))
              df <- df %>% left_join(df4, by = "code", keep = FALSE)
              
              
              #### Lastly, finding the cluster percentages for each chain configuration of each clone
              df <- df %>% mutate(clusters = .data[[clusters]])
              phen <- df %>% group_by(CTaa) %>% mutate(clone_count = n()) %>% ungroup() %>% group_by(CTaa, class, clusters) %>% summarise(count = n(),
                                                                                                                                          clone_count = first(clone_count))
              phen <- phen %>% group_by(CTaa, class) %>% mutate(clone_class_count = sum(count)) %>% ungroup()
              phen <- phen %>% mutate(pct = 100 * count / clone_class_count)
              phen2 <- phen %>% filter(clone_count >= cloneThresh & class != "Other" & clone_class_count >= cloneClassThresh) 

              
              phen2 <- phen2 %>% complete(CTaa, class, clusters)
              phen2 <- phen2 %>% group_by(CTaa, class) %>% mutate(combos = sum(!is.na(count))) %>% ungroup()
              phen2 <- phen2 %>% filter(combos != 0)
              phen2 <- phen2 %>% mutate(Pct = case_when(
                                                        !is.na(pct) ~ pct,
                                                        is.na(pct) ~ 0))
              
              return(list(df, phen2))

              
}



compareRNA <- function(tcr){

              plot1 <- tcr %>% ggplot(aes(x = class, y = nCount_RNA, fill = class)) + geom_violin() + geom_boxplot() + 
                               stat_compare_means(comparison = list(c("1\n&\n2", "2"), c("1", "2"), c("1\n&\n2", "NA"), c("1", "NA"), c("1\n&\n2", "1"))) + 
                               theme_classic() + xlab("Original Chain")
              print(plot1)
              
              plot2 <- tcr %>% filter(class != "Other") %>% ggplot(aes(x = class, y = nCount_RNA, fill = class)) + 
                                                            geom_violin() + geom_boxplot() + theme_classic() + 
                                                            xlab("Original Chain") + stat_compare_means(method = 'kruskal.test')
              print(plot2)
              
              tcr_filtered <- tcr %>% filter(class != "Other")
              model <- aov(tcr_filtered$nCount_RNA ~ tcr_filtered$class)
              plot3 <- plot(TukeyHSD(model))
              print(plot3)
}




clusterPermutationTest1v2 <- function(object, tcr, minCells, clusters){
              
              ### Formatting and calculating percentages of each cluster/chain pairing
              df <- tcr %>% separate(scR_CTaa, into = c("scr_alpha", "scr_beta"), sep = "_", remove = FALSE)
              
              new <- df %>% filter(scr_alpha != "NA" & scr_beta != "NA" & !str_detect(scr_alpha, ";") & !str_detect(scr_beta, ";") & class %in% c("1°", "2°"))
              new <- new %>% group_by(scR_CTaa, scr_alpha, scr_beta, CTaa, class, donor) %>% mutate(subset_size = n()) %>% ungroup() %>% 
                             group_by(scR_CTaa, scr_alpha, scr_beta, CTaa, class, donor, clusters) %>% summarise(clone_size = first(clonalFrequency),
                                                                                                                 subset_size = first(subset_size),
                                                                                                                 cluster_subset_prop = n() / first(subset_size))
              
              #### Because summarise only summarizes groups that exist, need to run complete and add the meta data so that all groups are considered, even those that are 0%
              new <- new %>% ungroup() %>% complete(clusters, class, CTaa, donor)
              new <- new %>% group_by(CTaa, donor) %>% filter(sum(clone_size, na.rm = T) > 0) %>% ungroup()
              new <- new %>% arrange(desc(clone_size)) %>% group_by(CTaa, donor) %>% mutate(clone_size = case_when(
                                                                                                                   is.na(clone_size) ~ first(clone_size),
                                                                                                                   !is.na(clone_size) ~ clone_size)) %>% ungroup()
              
              new <- new %>% arrange(desc(subset_size)) %>% group_by(CTaa, class, donor) %>% mutate(subset_size = case_when(
                                                                                                                            is.na(subset_size) ~ first(subset_size),
                                                                                                                            !is.na(subset_size) ~ subset_size)) %>% ungroup()
              
              new <- new %>% mutate(cluster_subset_prop = case_when(
                                                                    is.na(cluster_subset_prop) ~ 0,
                                                                    !is.na(cluster_subset_prop) ~ cluster_subset_prop))
              
              #### Removing classes that don't exist and ensuring each subset (i.e. clone 1 with bitypic major chain) is big enough
              new <- new %>% group_by(class, CTaa, donor) %>% filter(sum(cluster_subset_prop) > 0) %>% ungroup()
              new <- new %>% filter(subset_size >= minCells) %>% group_by(CTaa, donor) %>% filter(n_distinct(class) > 1) %>% ungroup()
              
              
              ### Calculating difference in cluster percentage between major and minor bytpic chain cells. Then adding up all the differences from all 13 clusters
              diff <- new %>% group_by(CTaa, clusters, donor) %>% summarise(difference = abs(cluster_subset_prop[class == "2°"] - cluster_subset_prop[class == "1°"]))
              diff <- diff %>% group_by(CTaa, donor) %>% summarise(total = sum(difference) * 0.5,
                                                                   count = n(),
                                                                   diff_clusters = sum(difference > 0))
              mean_diff <- mean(diff$total)
              
              
              #### Pre-processing to assess random pairings
              object@meta.data <- object@meta.data %>% mutate(clusters = .data[[clusters]])
              object@meta.data <- object@meta.data %>% separate(scR_CTaa, into = c("scr_alpha", "scr_beta"), sep = "_", remove = FALSE)
              rando <- object@meta.data %>% filter(scr_alpha != "NA" & scr_beta != "NA" & !str_detect(scr_alpha, ";") & !str_detect(scr_beta, ";")) %>%
                                            group_by(scR_CTaa, donor) %>% mutate(clone_count = n()) %>%
                                            group_by(scR_CTaa, clusters, donor) %>% summarise(prop = n() / first(clone_count),
                                                                                              clone_size = first(clone_count))
              
              rando <- rando %>% ungroup() %>% complete(scR_CTaa, clusters, donor)
              rando <- rando %>% group_by(scR_CTaa, donor) %>% filter(sum(clone_size, na.rm = T) > 0) %>% ungroup()
              rando <- rando %>% arrange(desc(clone_size)) %>% group_by(scR_CTaa, donor) %>% mutate(clone_size = case_when(
                                                                                                                           is.na(clone_size) ~ first(clone_size),
                                                                                                                           !is.na(clone_size) ~ clone_size)) %>% ungroup()
              
              rando <- rando %>% arrange(desc(prop)) %>% mutate(prop = case_when(
                                                                                 is.na(prop) ~ 0,
                                                                                 !is.na(prop) ~ prop))
              rando <- rando %>% filter(clone_size >= minCells)
              rando <- rando %>% unite(scR_CTaa, donor, sep = " - ", remove = FALSE, col = "key")
              
              if (n_distinct(new$clusters) != n_distinct(rando$clusters)) {
                print("Permuted number of clusters != observed number of clusters")
                return()
              }
              
              perm_diffs <- c()
              i <- 1
              
              ### 10,000x permutations of paired sets of clones, with the number of paired sets = number of 3 TCR clones
              while (length(perm_diffs) < 10000) {
                    
                    if(i %% 1000 == 0){
                      message(paste(i / 1000, "0% complete", sep = ""))
                    }
                    
                    rando2 <- rando %>% distinct(key) %>% mutate(index = sample.int(n())) %>% left_join(rando, by = "key")  
    
                    rando2 <- rando2 %>% mutate(group = factor((index + 1) %/% 2),
                                                fake_class = factor((index %% 2) + 1))
         
                    
                    ### This is key. The 2 clones need to be from the same donor for a fair comparison
                    rando2 <- rando2 %>% group_by(group) %>% filter(n_distinct(donor) == 1 & n_distinct(fake_class) == 2) %>% 
                                         ungroup() %>% filter(dense_rank(index) <= (2 * nrow(diff)))
                    
                    if (n_distinct(rando2$index) != 2 * nrow(diff)) {
                      next
                    }
                    
                    perm_diff <- rando2 %>% group_by(group, clusters, donor) %>% reframe(difference = abs(prop[fake_class == 2] - prop[fake_class == 1]))
                    perm_diff <- perm_diff %>% group_by(group, donor) %>% summarise(total = sum(difference) * 0.5,
                                                                                    count = n(),
                                                                                    diff_clusters = sum(difference > 0))
                    perm_diffs[i] <- mean(perm_diff$total)
                    i <- i + 1
              }
              

              return(list(perm_diffs, mean_diff))

}




clusterPermutationTestMissing <- function(object, tcr, minCells, clusters){
                  
                  
                  df <- object@meta.data %>% filter(!is.na(scR_CTaa)) %>% mutate(class = case_when(
                                                                                                   str_detect(scR_CTaa, "^NA_|_NA$") ~ "Missing",
                                                                                                   !str_detect(scR_CTaa, "^NA_|_NA$") ~ "Non-Missing"))
                  
                  new <- df %>% mutate(clusters = .data[[clusters]])
                  
                  new <- new %>% group_by(CTaa, class, donor) %>% mutate(subset_size = n()) %>% ungroup() %>% 
                                 group_by(CTaa, class, donor, clusters) %>% summarise(clone_size = first(clonalFrequency),
                                                                                      subset_size = first(subset_size),
                                                                                      cluster_subset_prop = n() / first(subset_size))
                  
                  new <- new %>% ungroup() %>% complete(clusters, class, CTaa, donor)
                  new <- new %>% group_by(CTaa, donor) %>% filter(sum(clone_size, na.rm = T) > 0) %>% ungroup()
                  new <- new %>% arrange(desc(clone_size)) %>% group_by(CTaa, donor) %>% mutate(clone_size = case_when(
                                                                                                                       is.na(clone_size) ~ first(clone_size),
                                                                                                                       !is.na(clone_size) ~ clone_size)) %>% ungroup()
                  
                  new <- new %>% arrange(desc(subset_size)) %>% group_by(CTaa, donor, class) %>% mutate(subset_size = case_when(
                                                                                                                                is.na(subset_size) ~ first(subset_size),
                                                                                                                                !is.na(subset_size) ~ subset_size)) %>% ungroup()
                  
                  new <- new %>% mutate(cluster_subset_prop = case_when(
                                                                        is.na(cluster_subset_prop) ~ 0,
                                                                        !is.na(cluster_subset_prop) ~ cluster_subset_prop))
                  
                  new <- new %>% group_by(class, CTaa, donor) %>% filter(sum(cluster_subset_prop) > 0) %>% ungroup()
                  new <- new %>% filter(subset_size >= minCells) %>% group_by(CTaa, donor) %>% filter(n_distinct(class) > 1) %>% ungroup()
                  
                  diff <- new %>% group_by(CTaa, clusters, donor) %>% reframe(difference = abs(cluster_subset_prop[class == "Missing"] - cluster_subset_prop[class == "Non-Missing"]))
                  diff <- diff %>% group_by(CTaa, donor) %>% summarise(total = sum(difference) * 0.5,
                                                                count = n())
                  mean_diff <- mean(diff$total)
                  
                  
                  object@meta.data <- object@meta.data %>% mutate(clusters = .data[[clusters]])
                  object@meta.data <- object@meta.data %>% separate(scR_CTaa, into = c("scr_alpha", "scr_beta"), sep = "_", remove = FALSE)
                  rando <- object@meta.data %>% group_by(scR_CTaa, donor) %>% mutate(clone_count = n()) %>%
                                                group_by(scR_CTaa, donor, clusters) %>% summarise(prop = n() / first(clone_count),
                                                                                                  clone_size = first(clone_count)) %>% filter(!is.na(scR_CTaa))
                  
                  rando <- rando %>% ungroup() %>% complete(scR_CTaa, clusters, donor)
                  rando <- rando %>% group_by(scR_CTaa, donor) %>% filter(sum(clone_size, na.rm = T) > 0) %>% ungroup()
                  rando <- rando %>% arrange(desc(clone_size)) %>% group_by(scR_CTaa, donor) %>% mutate(clone_size = case_when(
                                                                                                                               is.na(clone_size) ~ first(clone_size),
                                                                                                                               !is.na(clone_size) ~ clone_size)) %>% ungroup()
              
                  rando <- rando %>% arrange(desc(prop)) %>% mutate(prop = case_when(
                                                                                    is.na(prop) ~ 0,
                                                                                    !is.na(prop) ~ prop))
                  rando <- rando %>% filter(clone_size >= minCells)
                  
                  rando <- rando %>% mutate(class = case_when(
                                                              str_detect(scR_CTaa, "^NA_|_NA$") ~ "Missing",
                                                              !str_detect(scR_CTaa, "^NA_|_NA$") ~ "Non-Missing"))
                  
                  ### Finding the number of clones that are missing 1 chain for each donor
                  rando <- rando %>% group_by(donor) %>% mutate(n_missing = n_distinct(scR_CTaa[which(class == "Missing")])) %>% ungroup()
                  
                  ### Makes a key that is the clone, donor, and class      
                  rando <- rando %>% unite(scR_CTaa, donor, class, sep = " - ", remove = FALSE, col = "key")
                  
                  
                  if (n_distinct(new$clusters) != n_distinct(rando$clusters)) {
                    print("Permuted number of clusters != observed number of clusters")
                    return()
                  }
                  
                  perm_diffs <- c()
                  
              for (i in 1:10000) {
                
                    if (i %% 1000 == 0) {
                      message(paste(i / 1000, "0% complete", sep = ""))
                    }
                
                    ### For each class/donor combo, each clone gets a random index that is then added to the rando tibble
                    rando2 <- rando %>% distinct(key, donor, class) %>% 
                                        group_by(donor, class) %>% mutate(index = sample.int(n())) %>% ungroup() %>% select(-donor, -class) %>%
                                        left_join(rando, by = "key")
                    
                    ### Filters so that there are the same amount of missing and non-missing clones for each donor. 
                    ### Missing is less common, so we filter down to the number of missing so that the number is the same
                    rando2 <- rando2 %>% filter(index <= n_missing) %>% mutate(group = factor(index),
                                                                               class = factor(class))
    
                    ### For each group/cluster/donor the difference in the proportion between the random missing and ranodom non-missing clone are calculated
                    ### Then the TVD is found grouping adding up the differences and multiplying by 0.5 for each group/donor
                    perm_diff <- rando2 %>% group_by(group, clusters, donor) %>% reframe(difference = abs(prop[class == "Missing"] - prop[class == "Non-Missing"]))
                    perm_diff <- perm_diff %>% group_by(group, donor) %>% summarise(total = sum(difference) * 0.5,
                                                                                    count = n())
                    
                    ### Randomly filtering the tibble down to the number of clones in diff which is the observed
                    ### This ensures it is a far comparison in terms of n of clones
                    perm_diff <- perm_diff[sample(nrow(perm_diff), size = nrow(diff)), ]
                    
                    
                    if (nrow(perm_diff) != nrow(diff)) {
                      next
                    }
                    perm_diffs[i] <- mean(perm_diff$total)
                    
              }
              

              return(list(perm_diffs, mean_diff))

}




scRep <- function(seurat_object, file_paths = NULL, group, file = "filtered", folders = NULL, batch, original_barcode, mode = "standard", keepConfidentContigs = FALSE){
  
              ### Reading the data in
  if (file != "manual") {
                  contigs_list <- vector("list", length(folders))
    
        for (i in 1:length(folders)) {
            
             if (file == "filtered") {        
                  files <- list.files(folders[i], pattern = "filtered_contig_annotations.csv", full.names = TRUE)
             } else if (file == "all") {
                  files <- list.files(folders[i], pattern = "all_contig_annotations.csv", full.names = TRUE)  
             }    
                
          
                if (length(files) == 0) {
                  print("No contig files found, returning nothing")
                  return()
                }
                else if (length(files) > 1) {
                  print("Multiple contig files found, returning nothing")
                  return()
                }
                else if (length(files) == 1) {
                  contigs_list[[i]] <- read.csv(files[1])
                  contigs_list[[i]]$barcode <- paste(i, "_", contigs_list[[i]]$barcode, sep = "")
                }
  
        }
                
  } else if (file == "manual") {
                contigs_list <- vector("list", length(file_paths))
                
          for (i in 1:length(file_paths)) {
                
                contigs_list[[i]] <- read.csv(file_paths[i])

                if (nrow(contigs_list[[i]]) == 0) {
                  
                  message(paste("File path ", i, " did not lead to any contigs being loaded in. Check the path name and that it is a .csv file", sep = ""))
                  return()
                  
                }
                
                contigs_list[[i]]$barcode <- paste(i, "_", contigs_list[[i]]$barcode, sep = "")
          
          }
   }
            ### Joins together all batches and adds data to them about which sample/donor/batch each contig was from
                contigs_full <- bind_rows(contigs_list)
                
                if (keepConfidentContigs == TRUE) {
                  
                  contigs_full <- contigs_full %>% filter(!cdr3 %in% c("", "None") & !cdr3_nt %in% c("", "None"))
                    contigs_full <- contigs_full %>% filter(chain %in% c("TRA", "TRB", "TRG", "TRD"))
                    
                    contigs_full <- contigs_full %>% mutate(productive = "true",
                                                            is_cell = "true",
                                                            full_length = "true")
                }
                
                contigs_full <- contigs_full %>% filter(high_confidence %in% c("true", "True", "TRUE") &
                                                        full_length %in% c("true", "True", "TRUE") &
                                                        productive %in% c("true", "True", "TRUE") &
                                                        is_cell %in% c("true", "True", "TRUE"))
                                                          
                
                seurat_object$tcrEclipse_barcode <- paste(seurat_object@meta.data[[batch]], "_", seurat_object@meta.data[[original_barcode]], sep = "")
                seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)
               
                cell_data <- seurat_object@meta.data %>% select(tcrEclipse_barcode, unique(c(batch, group)), seurat_object_barcode)
                all <- contigs_full %>% left_join(cell_data, by = c("barcode" = "tcrEclipse_barcode"), keep = FALSE)
                all <- all %>% mutate(barcode = seurat_object_barcode) %>% select(-seurat_object_barcode)
                if("barcode" %in% colnames(seurat_object@meta.data)) {
                  seurat_object@meta.data <- seurat_object@meta.data %>% select(-barcode)
                }
                
                seurat_object@meta.data <- seurat_object@meta.data %>% dplyr::rename("barcode" = seurat_object_barcode)
                message("QC Filtering:")
                message(paste(nrow(all), " intial contigs found before filtering", sep = ""))
                
                ### Renaming barcodes, filtering to only cells in seurat object, removing cells without a CDR3, and renaming QC metrics from cellranger so these contigs aren't thrown out
                all <- all %>% filter(barcode %in% seurat_object$barcode)
  
              seurat_object$tcrEclipseGroup <- seurat_object@meta.data[[group]]
              contig_list <- createHTOContigList(all, 
                                                   seurat_object,
                                                   group.by = "tcrEclipseGroup")
              
              if (mode == "standard") {
                no_vdjd <- combineTCR(contig_list)
              }
              else if (mode == "noNaFilterMulti") {
                no_vdjd <- combineTCR(contig_list, removeNA = TRUE, filterMulti = TRUE)
              }
              else if (mode == "noNaNoMulti") {
                no_vdjd <- combineTCR(contig_list, removeNA = TRUE, removeMulti = TRUE)
              }
              else if (!mode %in% c("standard", "lecoz")){
                message("Unknown mode. Returning nothing")
                return()
              }
              
              obj <- combineExpression(no_vdjd, seurat_object, cloneCall = "aa", proportion = TRUE)
              return(obj)
}




vdjdive <- function(seurat_object, group, folders, donor, batch, original_barcode, write_folder) {
              
  
              contigs_list <- vector("list", length(folders))
    
        for (i in 1:length(folders)) {
          
                  files <- list.files(folders[i], pattern = "filtered_contig_annotations.csv", full.names = TRUE)
                
                
          
                if (length(files) == 0) {
                  print("No contig files found, returning nothing")
                  return()
                }
                else if (length(files) > 1) {
                  print("Multiple contig files found, returning nothing")
                  return()
                }
                else if (length(files) == 1) {
                  contigs_list[[i]] <- read.csv(files[1])
                  contigs_list[[i]]$barcode <- paste(i, "_", contigs_list[[i]]$barcode, sep = "")
                }
  
        }
                
                ### Joins together all batches and adds data to them about which sample/donor/batch each contig was from
                contigs_full <- bind_rows(contigs_list)
                
                seurat_object$tcrEclipse_barcode <- paste(seurat_object@meta.data[[batch]], "_", seurat_object@meta.data[[original_barcode]], sep = "")
                seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)
               
                cell_data <- seurat_object@meta.data %>% select(tcrEclipse_barcode, unique(c(batch, donor, group)), seurat_object_barcode)
                all <- contigs_full %>% left_join(cell_data, by = c("barcode" = "tcrEclipse_barcode"), keep = FALSE)
                all <- all %>% mutate(barcode = seurat_object_barcode) %>% select(-seurat_object_barcode)
                if("barcode" %in% colnames(seurat_object@meta.data)) {
                  seurat_object@meta.data <- seurat_object@meta.data %>% select(-barcode)
                }
                
                seurat_object@meta.data <- seurat_object@meta.data %>% dplyr::rename("barcode" = seurat_object_barcode)
                message("QC Filtering:")
                message(paste(nrow(all), " intial contigs found before filtering", sep = ""))
                
                ### Renaming barcodes, filtering to only cells in seurat object, removing cells without a CDR3, and renaming QC metrics from cellranger so these contigs aren't thrown out
                all <- all %>% filter(barcode %in% seurat_object$barcode)
  
                write.csv(all, paste(write_folder, "/filtered_contig_annotations.csv", sep = ""), row.names = FALSE)
                
                contigs <- readVDJcontigs(write_folder)
                vdj <- suppressMessages(clonoStats(contigs, method = "EM", type = "TCR", assignment = TRUE, group = donor))
                
                message("Extracting results from VDJdive...")
                df <- vdj@assignment
                df <- as(df, "RsparseMatrix")
                colnames(df) <- clonoNames(vdj)
                
                rm(vdj)
                
                #### Making a loop that that finds the most likely clone for each cell as well as the second most likely and the ratio of likelihoods between the two
                ### Putting all that data into a df called em
                barcode <- rownames(df)
                clone <- vector(length = nrow(df))
                maxes <- vector(length = nrow(df))
                second_max <- vector(length = nrow(df))
                
                for (i in 1:nrow(df)) {
                  
                  maxes[i] <- max(df[i, ])
                  clone[i] <- names(which.max(df[i, ]))
                  second_max[i] <- max(df[i, colnames(df) != clone[i]])
                  
                }
                
                rm(df)
                
                ### Summarizing the VDJdive EM data and then filtering cells that don't meet QC thresholds
                em <- data.frame(barcode, maxes, clone, second_max)
                em$ratio <- em$maxes / em$second_max
                em2 <- em %>% arrange(maxes) %>% mutate(nm = row_number(),
                                                        assigned_by_VDJdive = maxes >= 0.8 | (maxes >= 0.7 & ratio > 5) | (maxes > 0.5 & ratio > 10))
                plot2 <- em2 %>% ggplot() + geom_point(aes(x = nm, y = maxes, color = assigned_by_VDJdive)) +
                                            geom_hline(yintercept = 0.8) +
                                            xlab("Ranked Position") + ylab("Proportion Assigned to Most Likely Clone") + ggtitle("VDJdive EM Assignment QC")
                print(plot2)
                final <- em %>% filter(maxes >= 0.8 | (maxes >= 0.7 & ratio > 5) | (maxes > 0.5 & ratio > 10))
                
                final <- final %>% mutate(clone2 = str_replace(clone, " ", "_"))
                final <- final %>% select(barcode,
                                          "CTaa" = clone2,
                                          "EM_max_prop" = maxes,
                                          "EM_2ndmax_prop" = second_max,
                                          "VDJ_assign_ratio" = ratio)
                
                obj <- seurat_object
                obj@meta.data <- obj@meta.data %>% left_join(final, by = "barcode")
                
                ### Making columns so that scRepertoire visualizations can be used. Grouping by the group specified in the arguments, not the donor
                obj@meta.data <- obj@meta.data %>% group_by(.data[[group]]) %>% mutate(eclipseTcrGroupSize = sum(!is.na(CTaa))) %>% ungroup()
                obj@meta.data <- obj@meta.data %>% group_by(CTaa, .data[[group]]) %>% mutate(clonalProportion = case_when(
                                                                                                                          !is.na(CTaa) ~ n() / first(eclipseTcrGroupSize),
                                                                                                                           is.na(CTaa) ~ NA),
                                                                                             clonalFrequency = case_when(
                                                                                                                         !is.na(CTaa) ~ n(),
                                                                                                                          is.na(CTaa) ~ NA)) %>% ungroup() %>% select(-eclipseTcrGroupSize)
                obj@meta.data <- obj@meta.data %>% mutate(cloneSize = factor(case_when(
                                                                                       clonalProportion <= 1 & clonalProportion > 0.1 ~ "Hyperexpanded (0.1 < X <= 1)",
                                                                                       clonalProportion <= 0.1 & clonalProportion > 0.01 ~ "Large (0.01 < X <= 0.1)",
                                                                                       clonalProportion <= 0.01 & clonalProportion > 0.001 ~ "Medium (0.001 < X <= 0.01)",
                                                                                       clonalProportion <= 0.001 & clonalProportion > 1e-04 ~ "Small (1e-04 < X <= 0.001)",
                                                                                       clonalProportion <= 1e-04 & clonalProportion > 0 ~ "Rare (0 < X <= 1e-04)",
                                                                                       clonalProportion <= 0 ~ "None ( < X <= 0)",
                                                                                       is.na(clonalProportion) ~ NA)))
                
                return(obj)

}


pre <- function(folders = NULL, file_type = "all", file_paths = NULL, seurat_object, batch, donor, group, original_barcode, write_folder, mode = "standard", allowFourChains = TRUE, format = "Blank") {
  
  #### THIS FIRST PART IS READING THE DATA IN AND DOING QC
                ### Checking the input data is usable
                if (write_folder %in% folders & file_type != "manual") {
                  message("Folders cannot be the same as the write folder. This would result in the writing of a new .csv file over the raw data, removing the raw data from the computer.")
                  return()
                }
  
                if (is.null(folders) & file_type != "manual") {
                  message("Please provide folders that have the contigs annotations in them.")
                  return()
                }
                
                if (is.null(file_paths) & file_type == "manual") {
                  message("No file paths provided. Please provide a vector of file path names to file_paths")
                  return()
                }
  
                if (!batch %in% colnames(seurat_object@meta.data)) {
                  message("Batch column name provided is not in the Seurat object. Please try again.")
                  return()
                }
  
                if (!donor %in% colnames(seurat_object@meta.data)) {
                  message("Donor column name provided is not in the Seurat object. Please try again.")
                  return()
                }
  
                if (!group %in% colnames(seurat_object@meta.data)) {
                  message("Group column name provided is not in the Seurat object. Please try again.")
                  return()
                }
                
                if (!original_barcode %in% colnames(seurat_object@meta.data)) {
                  message("Original barcode column name provided is not in the Seurat object. Please try again.")
                  return()
                }
  
                if (!allowFourChains %in% c(TRUE, FALSE)) {
                  message("allowFourChains must be TRUE or FALSE. Default is TRUE.")
                  return()
                }
                
                if (!file_type %in% c("ALL", "All", "all", "Filtered", "FILTERED", "filtered", "all_contig_annotations.csv", "filtered_contig_annotations.csv", "manual")) {
                  print("Incorrect file type. Should be all or filtered. Returning nothing")
                  return()
                }
                
                ### Reading the data in
                
    if (file_type != "manual") {
                
                contigs_list <- vector("list", length(folders))
      
        for (i in 1:length(folders)) {
          
                if (file_type %in% c("ALL", "All", "all", "all_contig_annotations.csv")) {
                  files <- list.files(folders[i], pattern = "all_contig_annotations.csv", full.names = TRUE)
                }
                else if (file_type %in% c("Filtered", "filtered", "FILTERED", "filtered_contig_annotations.csv")) {
                  files <- list.files(folders[i], pattern = "filtered_contig_annotations.csv", full.names = TRUE)
                }
                
          
                if (length(files) == 0) {
                  print("No contig files found, returning nothing")
                  return()
                }
                else if (length(files) > 1) {
                  print("Multiple contig files found, returning nothing")
                  return()
                }
                else if (length(files) == 1) {
                  contigs_list[[i]] <- read.csv(files[1])
                  contigs_list[[i]]$barcode <- paste(i, "_", contigs_list[[i]]$barcode, sep = "")
                }
  
        }
      
    } else if (file_type == "manual") {
      
               contigs_list <- vector("list", length(file_paths))
        
        for (i in 1:length(file_paths)) {
                
                contigs_list[[i]] <- read.csv(file_paths[i])

                if (nrow(contigs_list[[i]]) == 0) {
                  
                  message(paste("File path ", i, " did not lead to any contigs being loaded in. Check the path name and that it is a .csv file", sep = ""))
                  return()
                  
                }
                
                contigs_list[[i]]$barcode <- paste(i, "_", contigs_list[[i]]$barcode, sep = "")
          
          }
          
    }
                
                ### Joins together all batches and adds data to them about which sample/donor/batch each contig was from
                contigs_full <- bind_rows(contigs_list)
                
                if("sample" %in% colnames(contigs_full)) {
                  
                  contigs_full <- contigs_full %>% select(-sample)
                  
                }
                
                if (format %in% c("None", "none")) {
                  contigs_full[contigs_full == "None"] <- ""
                }
                
                seurat_object$tcrEclipse_barcode <- paste(seurat_object@meta.data[[batch]], "_", seurat_object@meta.data[[original_barcode]], sep = "")
                seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)
               
                cell_data <- seurat_object@meta.data %>% select(tcrEclipse_barcode, unique(c(batch, donor, group)), seurat_object_barcode)
                all <- contigs_full %>% left_join(cell_data, by = c("barcode" = "tcrEclipse_barcode"), keep = FALSE)
                all <- all %>% mutate(barcode = seurat_object_barcode) %>% select(-seurat_object_barcode)
                if("barcode" %in% colnames(seurat_object@meta.data)) {
                  seurat_object@meta.data <- seurat_object@meta.data %>% select(-barcode)
                }
                
                seurat_object@meta.data <- seurat_object@meta.data %>% dplyr::rename("barcode" = seurat_object_barcode)
                message("QC Filtering:")
                message(paste(nrow(all), " intial contigs found before filtering", sep = ""))
                
                ### Renaming barcodes, filtering to only cells in seurat object, removing cells without a CDR3, and renaming QC metrics from cellranger so these contigs aren't thrown out
                all <- all %>% filter(barcode %in% seurat_object$barcode)
                message(paste(nrow(all), " contigs were from true cells", sep = ""))
                if(nrow(all) == 0){
                  message("Check that your barcodes are the original (i.e. likely end with -1 and match the contig files)")
                  return()
                }
                all <- all %>% filter(!cdr3 %in% c("", "None") & !cdr3_nt %in% c("", "None"))
                all <- all %>% filter(chain %in% c("TRA", "TRB", "TRG", "TRD"))
                message(paste(sum(all$high_confidence %in% c("true", "True", "TRUE") == FALSE), " contigs were removed for not being high confidence", sep = ""))
                all <- all %>% filter(high_confidence %in% c("true", "True", "TRUE"))
                all <- all %>% mutate(productive = "true",
                                      is_cell = "true",
                                      full_length = "true")
          
                ### Removing cells with abnormally long CDR3s
                all <- all %>% mutate(cdr3_length = str_length(cdr3))
                plot <- all %>% ggplot() + geom_histogram(aes(x = cdr3_length, fill = chain), bins = 30) + geom_vline(aes(xintercept = 30)) + xlab("CDR3 Length (aa)") + ylab("Count")
                print(plot)
                message(paste(sum(all$cdr3_length > 30), " contigs had a CDR3 longer than 30 amino acids and were removed", sep = ""))
                all <- all %>% filter(cdr3_length <= 30)
                
                ### There are weird TCRs with TRDV and then TRAJ and TRAC. Cellranger calls them as TRD, but they are actually TRA since they pair with a beta chain and have no J chain.
                all <- all %>% mutate(ad_hybrid_tcr = case_when(
                                                                chain %in% c("TRD", "Multi") & str_detect(v_gene, "TRDV") & (str_detect(j_gene, "TRAJ") | str_detect(c_gene, "TRAC")) ~ TRUE,
                                                                TRUE ~ FALSE))
                message(paste(sum(all$ad_hybrid_tcr), " contigs were found that were hybrid TRD/TRA chains. These were converted to TRA chains", sep = ""))
                all <- all %>% mutate(chain = case_when(
                                                        chain %in% c("TRD", "Multi") & str_detect(v_gene, "TRDV") & (str_detect(j_gene, "TRAJ") | str_detect(c_gene, "TRAC")) ~ "TRA",
                                                        TRUE ~ chain))
                message(paste(sum(all$chain %in% c("TRG", "TRD", "Multi")), " true TRG/TRD contigs were found and removed", sep = ""))
                all <- all %>% filter(chain %in% c("TRA", "TRB"))
                
                ### Removing cells with a * at the end of their cdr3. These indicate stop codons, so removing these contigs
                if (sum(str_detect(all$cdr3, "\\*")) > 0) {
                  message(paste("Warning, ",  sum(str_detect(all$cdr3, "\\*")), " contigs have a stop codon at the end of their CDR3 sequence, so removing them", sep = ""))
                  all <- all %>% filter(!str_detect(cdr3, "\\*"))
                }
                
                ### Removing duplicated chains. Some datasets (particularly contig files made with older version of Cellranger) have lots of duplicated rows that have identical CDR3 and gene segments for the same barcode. This messes up future analysis
                duplicates <- all %>% arrange(desc(umis), desc(reads)) %>% group_by(barcode, chain, cdr3) %>% filter(row_number() > 1) %>% ungroup()
                message(paste(nrow(duplicates), " duplicated contigs were removed", sep = ""))
                all <- all %>% arrange(desc(umis), desc(reads)) %>% group_by(barcode, chain, cdr3) %>% filter(row_number() == 1) %>% ungroup()
                
                return(list(all, seurat_object))
}


post <- function(folders, all, seurat_object, batch, donor, group, original_barcode, write_folder, mode = "standard", allowFourChains = TRUE, format = "Blank") {
  
  
                ### Noting cells with only 1 type of chain (alpha or beta) and that chain isn't seen in any other cell. 
                all2 <- all %>% group_by(barcode) %>% mutate(chainTypes = n_distinct(chain)) %>% ungroup()
                all2 <- all2 %>% group_by(cdr3, chain, .data[[donor]]) %>% mutate(maxChains = max(chainTypes)) %>% ungroup()
                message(paste(sum(all2$maxChains == 1), " orphan contigs were found with no observed pairing.", sep = ""))
                orphanContigs <- all2 %>% filter(maxChains == 1) %>% select(barcode, chain, cdr3)
                
                message(paste("After all filtering, ", nrow(all2), " contigs from ", n_distinct(all2$barcode), " true cells remained", sep = ""))
                   
                ### Gives better spacing
                message()
                message()
                message()
                
                message("Finding clones with extra chain...")
                message()
                
                donors <- unique(as.vector(all[, donor])[[1]])
                new_contigs <- vector("list", length(donors))
         for (i in 1:length(donors)) {
                  
                ### THIS PART IS FINDING THE COMBO CLONES AND REWRITING THE CONTIG FILES
                message(paste("Running donor ", i, "/", length(donors), " (", donors[i], ")", sep = ""))
                slim <- all2 %>% filter(.data[[donor]] == donors[i])
                
      
                ### Calling clones with multiple alpha chains and then rewriting the contig files to address this
                fused_a <- findSpecialClones2(contigs = slim, chain = "a")
    
                ### Calling clones with multiple beta chains and then rewriting the contig files to address this
                new_contigs[[i]] <- findSpecialClones2(contigs = fused_a, chain = "b")
                
         }
                
                contigs_post <- bind_rows(new_contigs)
                
   ### THIS PART IS FEEDING THE EDITED CONTIG FILES INTO VDJdive WHICH CALLS AMBIGUOUS CELLS
                
                message("Running VDJdive clonoStats() for clone EM assignment...")
                ### Writing a .csv file as filtered_contig_annotations.csv
                ### Reading that into VDJdive
                ### Using clonoStats for EM assignment
                ### Extracting output and changing format
                
                if (!write_folder %in% folders) {
                  write.csv(contigs_post, paste(write_folder, "/filtered_contig_annotations.csv", sep = ""), row.names = FALSE)
                }
                contigs <- readVDJcontigs(write_folder)
                vdj <- suppressMessages(clonoStats(contigs, method = "EM", type = "TCR", assignment = TRUE, group = donor))
                
                message("Extracting results from VDJdive...")
                df <- vdj@assignment
                df <- as(df, "RsparseMatrix")
                colnames(df) <- clonoNames(vdj)
                
                rm(vdj)
                
                #### Making a loop that that finds the most likely clone for each cell as well as the second most likely and the ratio of likelihoods between the two
                ### Putting all that data into a df called em
                barcode <- rownames(df)
                clone <- vector(length = nrow(df))
                maxes <- vector(length = nrow(df))
                second_max <- vector(length = nrow(df))
                
                for (i in 1:nrow(df)) {
                  
                  maxes[i] <- max(df[i, ])
                  clone[i] <- names(which.max(df[i, ]))
                  second_max[i] <- max(df[i, colnames(df) != clone[i]])
                  
                }
                
                rm(df)
                
                ### Summarizing the VDJdive EM data and then filtering cells that don't meet QC thresholds
                em <- data.frame(barcode, maxes, clone, second_max)
                em$ratio <- em$maxes / em$second_max
                em2 <- em %>% arrange(maxes) %>% mutate(nm = row_number(),
                                                        assigned_by_VDJdive = maxes >= 0.8 | (maxes >= 0.7 & ratio > 5) | (maxes > 0.5 & ratio > 10))
                plot2 <- em2 %>% ggplot() + geom_point(aes(x = nm, y = maxes, color = assigned_by_VDJdive)) +
                                            geom_hline(yintercept = 0.8) +
                                            xlab("Ranked Position") + ylab("Proportion Assigned to Most Likely Clone") + ggtitle("VDJdive EM Assignment QC")
                print(plot2)
                final <- em %>% filter(maxes >= 0.8 | (maxes >= 0.7 & ratio > 5) | (maxes > 0.5 & ratio > 10))  
                
                ### Removing cells where an orphan TCR was called. There's no way for the algorithm to accurately predict these so we don't want them.
                orphanContigs <- orphanContigs %>% left_join(final, by = "barcode")
                orphanContigs <- orphanContigs %>% mutate(remove = case_when(
                                                                             is.na(clone) ~ FALSE,
                                                                             chain == "TRA" & str_detect(clone, paste("^", cdr3, " ", sep = "")) ~ TRUE,
                                                                             chain == "TRB" & str_detect(clone, paste(" ", cdr3, "$", sep = "")) ~ TRUE,
                                                                             TRUE ~ FALSE))
                orphanContigs <- orphanContigs %>% filter(remove) %>% select(barcode, remove)
                orphanContigs <- unique(orphanContigs)
                
                final <- final %>% left_join(orphanContigs, by = "barcode")
                message(paste(sum(final$remove, na.rm = T), " cells were removed for being assigned an orphan chain", sep = ""))
                final <- final %>% filter(is.na(remove))
                
                ### Adding the VDJdive and scRepertoire calls to the Seurat object
                seurat_object$tcrEclipseGroup <- seurat_object@meta.data[[group]]
                contig_list <- createHTOContigList(all, 
                                                   seurat_object,
                                                   group.by = "tcrEclipseGroup")
                no_vdjd <- combineTCR(contig_list)
                obj <- combineExpression(no_vdjd, seurat_object, cloneCall = "aa", proportion = TRUE)
                rm(seurat_object)
                
                obj@meta.data <- obj@meta.data %>% rename(scR_CTgene = CTgene,
                                                          scR_CTnt = CTnt,
                                                          scR_CTaa = CTaa,
                                                          scR_CTstrict = CTstrict,
                                                          scR_clonalProportion = clonalProportion,
                                                          scR_clonalFrequency = clonalFrequency,
                                                          scR_cloneSize = cloneSize)
                final <- final %>% mutate(clone2 = str_replace(clone, " ", "_"))
                final <- final %>% mutate(clone2 = str_replace_all(clone2, ";", ":::"))
                final <- final %>% select(barcode,
                                          "VDJdive_clone" = clone2,
                                          "EM_max_prop" = maxes,
                                          "EM_2ndmax_prop" = second_max,
                                          "VDJ_assign_ratio" = ratio)
                
                obj@meta.data <- obj@meta.data %>% left_join(final, by = "barcode")
                
                ### Making CTaa our final clone call which uses VDJdive if available but if not the scRepertoire call
                obj@meta.data <- obj@meta.data %>% mutate(cloneCallSource = case_when(
                                                                                      !is.na(VDJdive_clone) ~ "VDJdive",
                                                                                      is.na(VDJdive_clone) & !is.na(scR_CTaa) ~ "scRepertoire",
                                                                                      TRUE ~ NA))
                if (allowFourChains == TRUE) {
                    obj@meta.data <- obj@meta.data %>% mutate(CTaa = case_when(
                                                                               !is.na(VDJdive_clone) ~ VDJdive_clone,
                                                                               is.na(VDJdive_clone) & !is.na(scR_CTaa) ~ scR_CTaa,
                                                                               TRUE ~ NA))
                } else if (alllowFourChains == FALSE) {
                    obj@meta.data <- obj@meta.data %>% mutate(CTaa = case_when(
                                                                               !is.na(VDJdive_clone) & str_count(VDJdive_clone, ":::") > 1 ~ scR_CTaa,
                                                                               !is.na(VDJdive_clone) & str_count(VDJdive_clone, ":::") <= 1 ~ VDJdive_clone,
                                                                               is.na(VDJdive_clone) & !is.na(scR_CTaa) ~ scR_CTaa,
                                                                               TRUE ~ NA))
                }
                
                obj@meta.data <- obj@meta.data %>% mutate(CTaa = str_replace_all(CTaa, ";", ":::"))
                
                ### Combining clones. Cells that were ambiguous may be assigned to a a combo clone by VDJdive but only be listed with one chain
                ### Here those cells have their CDR3 renamed as the combo CDR3
                obj@meta.data <- obj@meta.data %>% ungroup() %>% mutate(CTaa_donor = case_when(
                                                                                               !is.na(CTaa) & !is.na(.data[[donor]]) ~ paste(CTaa, .data[[donor]], sep = " "),
                                                                                                TRUE ~ NA))
 
                custom_a <- obj@meta.data %>% filter(!is.na(CTaa) & !is.na(.data[[donor]])) %>% group_by(CTaa, .data[[donor]], cloneCallSource) %>% summarise(count = n()) %>% arrange(desc(count))
                custom_a <- custom_a %>% separate(CTaa, into = c("alpha", "beta"), remove = FALSE, sep = "_")
                custom_a <- custom_a %>% mutate(special_a = case_when(
                                                                      str_detect(alpha, ":::") & cloneCallSource == "VDJdive" ~ alpha,
                                                                      !str_detect(alpha, ":::") | cloneCallSource != "VDJdive" ~ NA))
                
                custom_a <- custom_a %>% group_by(beta, .data[[donor]]) %>% arrange(desc(count)) %>% mutate(special_aChain = first(special_a)) %>% ungroup()
                custom_a <- custom_a %>% mutate(new_clone = case_when(
                                                                      str_detect(special_aChain, paste("^", alpha, ":::", "|", ":::", alpha, "$", sep = "")) & alpha != special_aChain ~ paste(special_aChain, beta, sep = "_"),
                                                                      TRUE ~ NA))
    
                custom_a <- custom_a %>% unite(col = "CTaa_donor", c(CTaa, .data[[donor]]), sep = " ", remove = FALSE)
                custom_a <- custom_a %>% filter(!is.na(new_clone)) %>% select(CTaa_donor, new_clone)
                obj@meta.data <- obj@meta.data %>% left_join(custom_a, by = "CTaa_donor")
                obj@meta.data <- obj@meta.data %>% mutate(CTaa = case_when(
                                                                           !is.na(new_clone) ~ new_clone,
                                                                            is.na(new_clone) ~ CTaa)) %>% select(-new_clone)
                
                custom_b <- obj@meta.data %>% filter(!is.na(CTaa) & !is.na(.data[[donor]])) %>% group_by(CTaa, .data[[donor]], cloneCallSource) %>% summarise(count = n()) %>% arrange(desc(count))
                custom_b <- custom_b %>% separate(CTaa, into = c("alpha", "beta"), remove = FALSE, sep = "_")
                custom_b <- custom_b %>% mutate(special_b = case_when(
                                                                      str_detect(beta, ":::") & cloneCallSource == "VDJdive" ~ beta,
                                                                      !str_detect(beta, ":::") | cloneCallSource != "VDJdive" ~ NA))
                custom_b <- custom_b %>% group_by(alpha, .data[[donor]]) %>% arrange(desc(count)) %>% mutate(special_bChain = first(special_b)) %>% ungroup()
                custom_b <- custom_b %>% mutate(new_clone = case_when(
                                                                      str_detect(special_bChain, paste("^", beta, ":::", "|", ":::", beta, "$", sep = "")) & beta != special_bChain ~ paste(alpha, special_bChain, sep = "_"),
                                                                      TRUE ~ NA))

                custom_b <- custom_b %>% unite(col = "CTaa_donor", c(CTaa, .data[[donor]]), sep = " ", remove = FALSE)
                custom_b <- custom_b %>% filter(!is.na(new_clone)) %>% select(CTaa_donor, new_clone)
                obj@meta.data <- obj@meta.data %>% left_join(custom_b, by = "CTaa_donor")
                obj@meta.data <- obj@meta.data %>% mutate(CTaa = case_when(
                                                                           !is.na(new_clone) ~ new_clone,
                                                                            is.na(new_clone) ~ CTaa)) %>% select(-new_clone, -CTaa_donor)
                
                ### Making columns so that scRepertoire visualizations can be used. Grouping by the group specified in the arguments, not the donor
                obj@meta.data <- obj@meta.data %>% group_by(.data[[group]]) %>% mutate(eclipseTcrGroupSize = sum(!is.na(CTaa))) %>% ungroup()
                obj@meta.data <- obj@meta.data %>% group_by(CTaa, .data[[group]]) %>% mutate(clonalProportion = case_when(
                                                                                                                          !is.na(CTaa) ~ n() / first(eclipseTcrGroupSize),
                                                                                                                           is.na(CTaa) ~ NA),
                                                                                             clonalFrequency = case_when(
                                                                                                                         !is.na(CTaa) ~ n(),
                                                                                                                          is.na(CTaa) ~ NA)) %>% ungroup() %>% select(-eclipseTcrGroupSize)
                obj@meta.data <- obj@meta.data %>% mutate(cloneSize = factor(case_when(
                                                                                       clonalProportion <= 1 & clonalProportion > 0.1 ~ "Hyperexpanded (0.1 < X <= 1)",
                                                                                       clonalProportion <= 0.1 & clonalProportion > 0.01 ~ "Large (0.01 < X <= 0.1)",
                                                                                       clonalProportion <= 0.01 & clonalProportion > 0.001 ~ "Medium (0.001 < X <= 0.01)",
                                                                                       clonalProportion <= 0.001 & clonalProportion > 1e-04 ~ "Small (1e-04 < X <= 0.001)",
                                                                                       clonalProportion <= 1e-04 & clonalProportion > 0 ~ "Rare (0 < X <= 1e-04)",
                                                                                       clonalProportion <= 0 ~ "None ( < X <= 0)",
                                                                                       is.na(clonalProportion) ~ NA)))
                
                
                ### Making a column that notes whether the cells in the clone appear to be MAIT or iNKT based on TCR segment usage
                obj@meta.data <- obj@meta.data %>% mutate(scR_mait_conv = case_when(
                                                                                    str_detect(scR_CTgene, "TRAV1-2\\.TRAJ33") ~ "MAIT Conventional (TRAV1-2 TRAJ33)",
                                                                                    TRUE ~ "Normal"),
                                                          scR_mait_unconv = case_when(
                                                                                      str_detect(scR_CTgene, "TRAV1-2\\.TRAJ12|TRAV1-2\\.TRAJ20") ~ "MAIT Unconventional (TRAV1-2 TRAJ12 or TRAJ20)",
                                                                                      TRUE ~ "Normal"),
                                                          scR_inkt = case_when(
                                                                               str_detect(scR_CTgene, "TRAV10\\.TRAJ18") & str_detect(scR_CTgene, "TRBV25-1") ~ "iNKT",
                                                                               TRUE ~ "Normal"))
                obj@meta.data <- obj@meta.data %>% mutate(scR_unconventional_count = 3 - (scR_mait_conv == "Normal") - (scR_mait_unconv == "Normal") - (scR_inkt == "Normal"))
                obj@meta.data <- obj@meta.data %>% mutate(cell_unconventional_subset = case_when(
                                                                                                 scR_unconventional_count == 0 ~ "Conventional",
                                                                                                 scR_unconventional_count == 1 & scR_mait_conv != "Normal" ~ scR_mait_conv,
                                                                                                 scR_unconventional_count == 1 & scR_mait_unconv != "Normal" ~ scR_mait_unconv,
                                                                                                 scR_unconventional_count == 1 & scR_inkt != "Normal" ~ scR_inkt,
                                                                                                 scR_unconventional_count > 1 ~ "Multiple Unconventional Types"))
                                                                          
                df <- obj@meta.data %>% group_by(CTaa, cell_unconventional_subset, .data[[group]]) %>% summarise(count = n()) %>% ungroup() %>% filter(cell_unconventional_subset != "Conventional") %>% arrange(desc(count))
                
                if (nrow(df) > 0) {
                     obj@meta.data <- obj@meta.data %>% ungroup() %>% mutate(CTaa_group = case_when(
                                                                                                    !is.na(CTaa) & !is.na(.data[[group]]) ~ paste(CTaa, .data[[group]], sep = " "),
                                                                                                     TRUE ~ NA))
                     df <- df %>% group_by(CTaa, .data[[group]]) %>% summarise(clone_unconventional_subset = first(cell_unconventional_subset))
                     df <- df %>% ungroup() %>% mutate(CTaa_group = case_when(
                                                                              !is.na(CTaa) & !is.na(.data[[group]]) ~ paste(CTaa, .data[[group]], sep = " "),
                                                                               TRUE ~ NA))
                     df <- df %>% ungroup() %>% select(CTaa_group, clone_unconventional_subset)
                     obj@meta.data <- obj@meta.data %>% left_join(df, by = "CTaa_group")
                     obj@meta.data <- obj@meta.data %>% mutate(clone_unconventional_subset = case_when(
                                                                                                        is.na(clone_unconventional_subset) ~ "Conventional",
                                                                                                       !is.na(clone_unconventional_subset) ~ clone_unconventional_subset)) %>% select(-CTaa_group)
                     message(paste(sum(obj$clone_unconventional_subset != "Conventional"), " cells were detected in clones that appear to be MAIT or iNKT cells", sep = ""))
                }
                
                obj@meta.data <- obj@meta.data %>% select(-scR_mait_conv, -scR_mait_unconv, -scR_inkt, -scR_unconventional_count, -tcrEclipseGroup, -tcrEclipse_barcode)
                
                obj@meta.data <- as.data.frame(obj@meta.data)
                rownames(obj@meta.data) <- obj$barcode
                
                message(paste(round(sum(obj$cloneCallSource == "VDJdive", na.rm = T) / nrow(obj@meta.data) * 100, 1), "% of cells in the Seurat object were assigned with high confidence by VDJdive clonoStats()", sep = ""))
                message(paste(round(sum(!is.na(obj$CTaa)) / nrow(obj@meta.data) * 100, 1), "% of cells in the Seurat object were assigned in the end (VDJdive or scRepertoire assignment)", sep = ""))
                message("These final clonal annotations are stored in CTaa")
                message("Summary of the top 5 largest clones:")
                big_clones <- obj@meta.data %>% group_by(.data[[donor]]) %>% mutate(donor_count = sum(!is.na(CTaa))) %>% ungroup() %>% 
                                                group_by(CTaa, .data[[donor]]) %>% summarise(Clone_Size = n(), Percent_of_Repertoire = 100 * n() / first(donor_count)) %>%
                                                select(CTaa, Clone_Size, Percent_of_Repertoire) %>% arrange(desc(Clone_Size)) %>% filter(!is.na(CTaa)) %>% head(n = 5)
                print(big_clones)
                
                
                
                return(obj@meta.data)
          
}




singlePairedTest <- function(data, group1, group2, makePaired = TRUE) {
  
  data <- as.data.frame(data)
  if (makePaired == TRUE) {
    data <- data[!is.na(data[,group1]) & !is.na(data[, group2]), ]
  }
  
  test <- wilcox.test(data[, group1], data[, group2], paired = makePaired)
  return(test)
}



runPairedTests <- function(data, groups, adjustment_method, runPaired = TRUE) {
  
  results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(results) <- c("group1", "group2", "p")
  
  for (i in 1:(length(groups) - 1)) {
    
    for (j in (i + 1):length(groups)) {
      
      pval <- singlePairedTest(data = data, group1 = groups[i], group2 = groups[j], makePaired = runPaired)$p.value
      results[nrow(results) + 1, ]<- c(groups[i], groups[j], pval)
      
    }
  }
  
  results$p.adj <- p.adjust(results$p, method = adjustment_method)
  results$p.adj2 <- signif(results$p.adj, 2)
  
  return(results)
}
