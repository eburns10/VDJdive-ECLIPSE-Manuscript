

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
              print(paste0("Observed n = ", nrow(diff)))
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
                  print(paste0("Observed n = ", nrow(diff)))
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




scRep <- function(seurat_object, file_paths = NULL, group, file = "filtered", folders = NULL, batch, original_barcode, mode = "standard", keepConfidentContigs = FALSE, cc = "aa"){
  
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
              
              obj <- combineExpression(no_vdjd, seurat_object, cloneCall = cc, proportion = TRUE)
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










getAssignmentMatrix <- function(seurat_object,
                    contig_type = "all",
                    contig_folder_paths = NULL,
                    contig_file_paths = NULL,
                    batch,
                    donor,
                    group,
                    original_barcode,
                    format = "Blank",
                    detailedEclipseStats = FALSE,
                    combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
                    allowDoubleChains = "both",
                    returnFilteredContigs = FALSE,
                    verbose = TRUE) {
    #### THIS FIRST PART IS READING THE DATA IN AND DOING QC
    ### Checking the input data is usable
    if (!contig_type %in% c(
        "ALL", "All", "all",
        "Filtered", "FILTERED", "filtered",
        "all_contig_annotations.csv",
        "filtered_contig_annotations.csv",
        "manual"
    )) {
        stop(
            "Incorrect contig type. Should be all or filtered. ",
            "Returning nothing"
        )
    }

    if (is.null(contig_folder_paths) & contig_type != "manual") {
        stop(
            "Please provide folders to contig_folder_paths that have the ",
            "contig annotations in them."
        )
    }

    if (is.null(contig_file_paths) & contig_type == "manual") {
        stop(
            "No contig file paths provided. ",
            "Please provide a vector of file path names to contig_file_paths"
        )
    }

    if (contig_type != "manual" &
        (!is.vector(contig_folder_paths) | is.list(contig_folder_paths) |
        length(contig_folder_paths) == 0)) {
        stop("contig_folder_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (contig_type == "manual" &
        (!is.vector(contig_file_paths) | is.list(contig_file_paths) |
        length(contig_file_paths) == 0)) {
        stop("contig_file_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (!is.null(contig_file_paths) & !is.null(contig_folder_paths)) {
        stop("Please do not use contig_file_paths and contig_folder_paths. ",
            "Please try again.")
    }

    if (!batch %in% colnames(seurat_object@meta.data)) {
        stop(
            "Batch column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!donor %in% colnames(seurat_object@meta.data)) {
        stop(
            "Donor column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!group %in% colnames(seurat_object@meta.data)) {
        stop(
            "Group column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!original_barcode %in% colnames(seurat_object@meta.data)) {
        stop(
            "Original barcode column name provided is not ",
            "in the Seurat object. ",
            "Please try again."
        )
    }

    if (!is.logical(detailedEclipseStats)) {
        stop(
            "TRUE or FALSE not provided to detailedEclipseStats. ",
            "Please try again."
        )
    }

    if (format %in% c("None", "none", "Blank") == FALSE) {
        stop(
            "None, none, or Blank not provided to format. ",
            "Please try again."
        )
    }

    if (!is.list(combineTcrArgs) | length(combineTcrArgs) != 4 |
        !is.logical(combineTcrArgs[[1]]) | !is.logical(combineTcrArgs[[2]]) |
        !is.logical(combineTcrArgs[[3]]) | !is.logical(combineTcrArgs[[4]])) {
        stop(
            "Incorrect input to combineTcrArgs. ",
            "Please provide a list of 4 logical (TRUE or FALSE) elements."
        )
    }

    if (allowDoubleChains %in% c(
        "Alpha", "alpha",
        "beta", "Beta",
        "Both", "both",
        "Neither", "neither"
    ) == FALSE) {
        stop(
            "Alpha, Beta, Both, or Neither not provided to ",
            "allowDoubleChains. ",
            "Please try again."
        )
    }

    if (!is.logical(verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )
    }

    if (!is.logical(returnFilteredContigs)) {
        stop(
            "TRUE or FALSE not provided to returnFilteredContigs. ",
            "Please try again."
        )
    }


    ### Reading the data in
    if (verbose == TRUE) message("\nReading TCR contig data in...\n")

    if (contig_type != "manual") {
        contigs_list <- vector("list", length(contig_folder_paths))

        for (i in seq_along(contig_folder_paths)) {
            if (contig_type %in% c(
                "ALL", "All", "all",
                "all_contig_annotations.csv"
            )) {
                files <- list.files(contig_folder_paths[i],
                    pattern = "all_contig_annotations.csv",
                    full.names = TRUE
                )
            } else if (contig_type %in%
                c(
                    "Filtered", "filtered", "FILTERED",
                    "filtered_contig_annotations.csv"
                )) {
                files <- list.files(contig_folder_paths[i],
                    pattern = "filtered_contig_annotations.csv",
                    full.names = TRUE
                )
            }


            if (length(files) == 0) {
                stop("No contig files found, returning nothing")
            } else if (length(files) > 1) {
                stop("Multiple contig files found, returning nothing")
            } else if (length(files) == 1) {
                contigs_list[[i]] <- read.csv(files[1])
                contigs_list[[i]]$barcode <- paste(i, "_",
                    contigs_list[[i]]$barcode,
                    sep = ""
                )
            }
        }
    } else if (contig_type == "manual") {
        contigs_list <- vector("list", length(contig_file_paths))

        for (i in seq_along(contig_file_paths)) {
            contigs_list[[i]] <- read.csv(contig_file_paths[i])

            if (nrow(contigs_list[[i]]) == 0) {
                stop(
                    "File path ", i, " did not lead to any ",
                    "contigs being loaded in. ",
                    "Check the path name and that it is a .csv file"
                )
            }

            contigs_list[[i]]$barcode <-
                paste(i, "_", contigs_list[[i]]$barcode, sep = "")
        }
    }

    ### Joins together all batches and adds data to them
    ### about which sample/donor/batch each contig was from
    contigs_full <- bind_rows(contigs_list)

    if ("sample" %in% colnames(contigs_full)) {
        contigs_full <- contigs_full %>% select(-all_of("sample"))
    }

    if (format %in% c("None", "none")) {
        contigs_full[contigs_full == "None"] <- ""
    }

    seurat_object$tcrEclipse_barcode <-
        paste(seurat_object@meta.data[[batch]],
            "_",
            seurat_object@meta.data[[original_barcode]],
            sep = ""
        )
    seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)

    cell_data <- seurat_object@meta.data %>%
        select(all_of(c(
            "tcrEclipse_barcode",
            unique(c(batch, donor, group)),
            "seurat_object_barcode"))
        )
    cell_data <- cell_data %>% filter(!is.na(.data[[batch]]))

    if (nrow(cell_data) != n_distinct(cell_data$tcrEclipse_barcode)) {
        stop(
            "There are likely multiple sequencing runs ",
            "labeled as one batch. ",
            "Please relabel the batch column such that each batch ",
            "is a unique sequncing run. ",
            "Returning nothing"
        )
    }

    all <- contigs_full %>%
        left_join(cell_data,
            by = c("barcode" = "tcrEclipse_barcode"),
            keep = FALSE
        )
    all <- all %>%
        mutate(barcode = .data$seurat_object_barcode) %>%
        select(-all_of("seurat_object_barcode"))
    if ("barcode" %in% colnames(seurat_object@meta.data)) {
        seurat_object@meta.data <- seurat_object@meta.data %>%
            mutate(originalBarcodeBeforeECLIPSE = .data$barcode) %>%
            select(-all_of("barcode"))
    }

    seurat_object@meta.data <- seurat_object@meta.data %>%
        dplyr::rename("barcode" = all_of("seurat_object_barcode"))

    if (verbose == TRUE) {
        message("QC Filtering:")
        message(" - ", nrow(all), " initial contigs found before filtering")
    }

    ### Renaming barcodes, filtering to only cells in seurat object,
    ### removing cells without a CDR3,
    ### and renaming QC metrics from cellranger so
    ### these contigs aren't thrown out
    all <- all %>% filter(.data$barcode %in% seurat_object$barcode)
    if (verbose == TRUE) {
        message(" - ", nrow(all), " contigs were from true cells ",
            "(coming from ", n_distinct(all$barcode), " cells)")
    }
    if (nrow(all) == 0) {
        stop(
            "Check that your barcodes are the original ",
            "(i.e. likely end with -1 and match the contig files)"
        )
    }

    if (verbose == TRUE) {
        message(" - ",
            sum(all$high_confidence %in% c("true", "True", "TRUE") == FALSE),
            " contigs were removed for not being high confidence as ",
            "defined by Cell Ranger"
        )
    }
    all <- all %>%
        filter(.data$high_confidence %in% c("true", "True", "TRUE"))
    all <- all %>% mutate(
        productive = "true",
        is_cell = "true",
        full_length = "true"
    )

    count_before_custom_filtering <- nrow(all)

    all <- all %>% filter(!.data$cdr3 %in% c("", "None") &
                            !.data$cdr3_nt %in% c("", "None"))
    all <- all %>% filter(.data$chain %in% c("TRA", "TRB", "TRG", "TRD"))

    ### Removing cells with abnormally long CDR3s
    all <- all %>% mutate(cdr3_length = str_length(.data$cdr3))

    all <- all %>% filter(.data$cdr3_length <= 30)

    ### There are TCRs with TRDV and then TRAJ and TRAC.
    ### Cellranger calls them as TRD,
    ### but they are actually TRA since they pair with a beta chain and
    ### have no J segment.
    all <- all %>%
        mutate(ad_hybrid_tcr = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC"))
            ~ TRUE,
            TRUE ~ FALSE
        ))

    all <- all %>%
        mutate(chain = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC")) ~ "TRA",
            TRUE ~ .data$chain
        ))

    all <- all %>% filter(.data$chain %in% c("TRA", "TRB"))

    ### Removing cells with a * at the end of their cdr3.
    ### These indicate stop codons, so removing these contigs
    all <- all %>% filter(!str_detect(.data$cdr3, "\\*"))

    ### Removing duplicated chains. S
    ### Some datasets (particularly contig files made with
    ### older version of Cellranger)
    ### have lots of duplicated rows that have identical CDR3
    ### and gene segments for the same barcode.
    ### This messes up future analysis
    all <- all %>%
        arrange(desc(.data$umis), desc(.data$reads)) %>%
        group_by(.data$barcode, .data$chain, .data$cdr3) %>%
        filter(row_number() == 1) %>%
        ungroup()

    ### Noting cells with only 1 type of chain (alpha or beta) and
    ### that chain isn't seen in any other cell.
    ### We don't want these cells to be predicted.
    all2 <- all %>%
        group_by(.data$barcode) %>%
        mutate(chainTypes = n_distinct(.data$chain)) %>%
        ungroup()
    all2 <- all2 %>%
        group_by(.data$cdr3, .data$chain, .data[[donor]]) %>%
        mutate(maxChains = max(.data$chainTypes)) %>%
        ungroup()
    if (verbose == TRUE) {
        message(" - ",
            sum(all2$maxChains == 1),
            " orphan contigs were found with no observed pairing. ",
            "Cells containing these contigs are unlikely to have ",
            "a clonotype predicted"
        )
    }
    orphanContigs <- all2 %>%
        filter(.data$maxChains == 1) %>%
        select(all_of(c("barcode", "chain", "cdr3")))

    if (verbose == TRUE) {
        message(" - ", count_before_custom_filtering - nrow(all),
            " contigs were removed for other reasons ",
            "(ex. irregular CDR3, TCR\u03B3/\u03B4 contig, duplicated, etc...)"
        )
        message(" - ",
            "After all filtering, ", nrow(all),
            " contigs from ", n_distinct(all$barcode),
            " true cells remained"
        )
    }

    if (returnFilteredContigs == TRUE) {

        all2$orphanContig <- FALSE
        all2$orphanContig[all2$maxChains == 1] <- TRUE

        all2 <- all2 %>%
            select(-all_of(c(
                "cdr3_length",
                "ad_hybrid_tcr",
                "chainTypes",
                "maxChains"))
                )
        return(all2)
    }

    rm(contigs_full, contigs_list)

    if (verbose == TRUE) message("\nFinding clones with extra chains...\n")

    donors <- unique(as.vector(all[, donor])[[1]])
    new_contigs <- vector("list", length(donors))

    for (i in seq_along(donors)) {
        ### THIS PART IS FINDING THE COMBO CLONES AND
        ### REWRITING THE CONTIG FILES
        if (allowDoubleChains %in% c("Neither", "neither") == FALSE &
            verbose == TRUE) {
            message(
                "Running donor ", i, "/",
                length(donors), " (", donors[i], "):"
            )
        }

        slim <- all2 %>% filter(.data[[donor]] == donors[i])

        if (allowDoubleChains %in% c("Alpha", "alpha", "Both", "both")) {
            ### Calling clones with multiple alpha chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "alpha",
                printStats = verbose
            )
        }

        if (allowDoubleChains %in% c("Beta", "beta", "Both", "both")) {
            ### Calling clones with multiple beta chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "beta",
                printStats = verbose
            )
        }

        new_contigs[[i]] <- slim

        if (verbose == TRUE) message()
    }

    contigs_post <- bind_rows(new_contigs)

    ### THIS PART IS FEEDING THE EDITED CONTIG FILES INTO VDJdive
    ### WHICH CALLS AMBIGUOUS CELLS

    if (verbose == TRUE) {
        message("Running VDJdive::clonoStats() for clone EM assignment...")
    }
    ### Making a temporary directory
    ### Writing a .csv file as filtered_contig_annotations.csv
    ### Reading that into VDJdive
    ### Using clonoStats for EM assignment
    ### Extracting output and changing format

    temp_dir <- tempfile("ECLIPSE_")
    dir.create(temp_dir)

    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

    write.csv(contigs_post,
        paste(temp_dir, "/filtered_contig_annotations.csv", sep = ""),
        row.names = FALSE
    )

    contigs <- readVDJcontigs(temp_dir)
    
    ### Please note suppressMessages is used here because
    ### there is a message stating internal coercion from
    ### Tsparse* to CsparseMatrix
    ### that I do not want to be printed.
    vdj <- suppressMessages(clonoStats(contigs,
        method = "EM",
        type = "TCR",
        assignment = TRUE,
        group = donor
    ))

    if (verbose == TRUE) message("Extracting results from VDJdive...")
    df <- vdj@assignment
    df <- as(df, "RsparseMatrix")
    colnames(df) <- clonoNames(vdj)
    return(df)
    rm(vdj)

    ### Finding the most likely clonotype, its probability,
    ### and the probability of the 2nd most likely clonotype for each cell.
    p <- df@p
    j <- df@j
    x <- df@x

    all_clones <- colnames(df)
    n <- nrow(df)

    barcode <- rownames(df)
    clone <- character(n)
    maxes <- numeric(n)
    second_max <- numeric(n)


    for (i in seq_len(n)) {
        idx <- (p[i] + 1):p[i + 1]

        vals <- x[idx]
        cols <- j[idx] + 1

        if (length(vals) == 1) {
            maxes[i] <- vals
            second_max[i] <- 0
            clone[i] <- all_clones[cols]
        } else {
            best <- vals[1]
            second <- 0
            best_col <- cols[1]

            for (k in 2:length(vals)) {
                if (vals[k] > best) {
                    second <- best
                    best <- vals[k]
                    best_col <- cols[k]
                } else if (vals[k] > second) {
                    second <- vals[k]
                }
            }

            maxes[i] <- best
            second_max[i] <- second
            clone[i] <- all_clones[best_col]
        }
    }


    ### Summarizing the VDJdive EM data per cell
    ### and then filtering cells that don't meet QC thresholds
    em <- data.frame(barcode = rownames(df), maxes, clone, second_max)
    rm(df, p, j, x)

    if (verbose == TRUE) {
        message("Finalizing clonal assignments...\n")
    }
    em$ratio <- em$maxes / em$second_max

    if (verbose == TRUE) {
        em2 <- em %>%
            arrange(.data$maxes) %>%
            mutate(
                nm = row_number(),
                assigned_by_VDJdive = .data$maxes >= 0.8 |
                    (.data$maxes >= 0.7 & .data$ratio > 5) |
                    (.data$maxes > 0.5 & .data$ratio > 10)
            )
        plot2 <- em2 %>%
            ggplot() +
            geom_point(aes(
                x = .data$nm, y = .data$maxes,
                color = .data$assigned_by_VDJdive
            )) +
            xlab("Ranked position of cell") +
            ylab(paste("Probability of cell being\n",
                "in most likely clone", sep = "")) +
            ggtitle("Clonal Prediction QC") +
            theme_bw() +
            labs(color = "Source\nof Final\nClonotype") +
            scale_color_manual(values = c("cyan4", "darkorange2"),
                labels = c("Traditional\nAnalysis\n",
                    "High\nConfidence\nPrediction")) +
            theme(
                plot.title = element_text(
                    face = "bold", hjust = 0.5,
                    size = 20
                ),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 16, hjust = 0.5),
                legend.text = element_text(size = 12),
                legend.position.inside = c(0.7, 0.4)
            )
        plot(plot2)
    }

    final <- em %>% filter(.data$maxes >= 0.8 |
        (.data$maxes >= 0.7 & .data$ratio > 5) |
        (.data$maxes > 0.5 & .data$ratio > 10))

    ### Removing cells where an orphan TCR was assigned.
    ### There's no way for the algorithm to accurately predict these so
    ### we don't want them.
    orphanContigs <- orphanContigs %>% left_join(final, by = "barcode")
    orphanContigs <- orphanContigs %>%
        mutate(remove = case_when(
            is.na(.data$clone) ~ FALSE,
            .data$chain == "TRA" &
                str_detect(
                    .data$clone,
                    paste("^", .data$cdr3, " ", sep = "")
                ) ~ TRUE,
            .data$chain == "TRB" &
                str_detect(
                    .data$clone,
                    paste(" ", .data$cdr3, "$", sep = "")
                ) ~ TRUE,
            TRUE ~ FALSE
        ))
    orphanContigs <- orphanContigs %>%
        filter(.data$remove) %>%
        select(all_of(c("barcode", "remove")))
    orphanContigs <- unique(orphanContigs)

    final <- final %>% left_join(orphanContigs, by = "barcode")
    final <- final %>% filter(is.na(.data$remove))
    
    ### Adding the VDJdive and scRepertoire calls and
    ### expressed chains to the Seurat object
    seurat_object$tcrEclipseGroup <- seurat_object@meta.data[[group]]
    
    contig_list <- createHTOContigList(all,
        seurat_object,
        group.by = "tcrEclipseGroup"
    )
    
    no_vdjd <- combineTCR(contig_list)
    
    seurat_object <- combineExpression(no_vdjd,
        seurat_object,
        cloneCall = "aa",
        proportion = TRUE
    )
    
    seurat_object$highConfChains <- seurat_object$CTaa
    seurat_object@meta.data <- seurat_object@meta.data %>%
        select(-all_of(c(
            "CTgene",
            "CTnt",
            "CTaa",
            "CTstrict",
            "clonalProportion",
            "clonalFrequency",
            "cloneSize"))
        )
    
    no_vdjd <- combineTCR(contig_list,
        removeNA = combineTcrArgs[[1]],
        removeMulti = combineTcrArgs[[2]],
        filterMulti = combineTcrArgs[[3]],
        filterNonproductive = combineTcrArgs[[4]]
    )
    obj <- combineExpression(no_vdjd,
        seurat_object,
        cloneCall = "aa",
        proportion = TRUE
    )
    rm(seurat_object)
    
    columns_to_remove <- which(colnames(obj@meta.data) %in% 
        c("scR_CTgene", "scR_CTnt", "scR_CTaa", "scR_CTstrict", 
            "scR_clonalFrequency", "scR_clonalProportion", "scR_cloneSize"))
    
    if (length(columns_to_remove) > 0) {
      
        obj@meta.data <- obj@meta.data[, -columns_to_remove]       
    
    }
    
    obj@meta.data <- obj@meta.data %>%
        rename(
            scR_CTgene = all_of("CTgene"),
            scR_CTnt = all_of("CTnt"),
            scR_CTaa = all_of("CTaa"),
            scR_CTstrict = all_of("CTstrict"),
            scR_clonalProportion = all_of("clonalProportion"),
            scR_clonalFrequency = all_of("clonalFrequency"),
            scR_cloneSize = all_of("cloneSize")
        )
    final <- final %>% mutate(clone2 = str_replace(.data$clone, " ", "_"))
    final <- final %>% mutate(clone2 = str_replace_all(
        .data$clone2,
        ";", ":::"
    ))
    final <- final %>% select(all_of("barcode"),
        "VDJdive_clone" = all_of("clone2"),
        "EM_max_prop" = all_of("maxes"),
        "EM_2ndmax_prop" = all_of("second_max"),
        "EM_max_2ndmax_ratio" = all_of("ratio")
    )
    
    obj@meta.data <- obj@meta.data %>% left_join(final, by = "barcode")
    

    ### Making CTaa our final clone call which uses VDJdive
    ### if available but if not the scRepertoire call
    obj@meta.data <- obj@meta.data %>%
        mutate(cloneCallSource = case_when(
            !is.na(.data$VDJdive_clone) ~ "VDJdive",
            is.na(.data$VDJdive_clone) & !is.na(.data$scR_CTaa) ~
                "scRepertoire",
            TRUE ~ NA
        ))

    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = case_when(
            !is.na(.data$VDJdive_clone) ~ .data$VDJdive_clone,
            is.na(.data$VDJdive_clone) & !is.na(.data$scR_CTaa) ~
                .data$scR_CTaa,
            TRUE ~ NA
        ))


    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = str_replace_all(.data$CTaa, ";", ":::"))
    
    ### Combining clones. Cells that were ambiguous may be assigned to a
    ### combo clone by VDJdive but only be listed with one chain
    ### Here those cells have their CDR3 renamed as the combo CDR3
    obj@meta.data <- obj@meta.data %>%
        ungroup() %>%
        mutate(CTaa_donor_source = case_when(
            !is.na(.data$CTaa) &
                !is.na(.data[[donor]]) &
                !is.na(.data$cloneCallSource) ~
                paste(.data$CTaa, .data[[donor]],
                    .data$cloneCallSource,
                    sep = " "
                ),
            TRUE ~ NA
        ))

    custom_a <- obj@meta.data %>%
        filter(!is.na(.data$CTaa) & !is.na(.data[[donor]]) &
            .data$cloneCallSource == "VDJdive") %>%
        group_by(.data$CTaa, .data[[donor]], .data$cloneCallSource) %>%
        summarise(count = n()) %>%
        arrange(desc(.data$count))
    custom_a <- custom_a %>%
        separate(.data$CTaa,
            into = c("alpha", "beta"),
            remove = FALSE, sep = "_"
        )
    custom_a <- custom_a %>%
        mutate(special_a = case_when(
            str_detect(.data$alpha, ":::") &
                .data$cloneCallSource == "VDJdive" ~ .data$alpha,
            !str_detect(.data$alpha, ":::") |
                .data$cloneCallSource != "VDJdive"
            ~ NA
        ))

    custom_a <- custom_a %>%
        group_by(.data$beta, .data[[donor]]) %>%
        arrange(desc(.data$count)) %>%
        mutate(special_aChain = first(.data$special_a)) %>%
        ungroup()
    custom_a <- custom_a %>%
        mutate(new_clone = case_when(
            str_detect(
                .data$special_aChain,
                paste("^", .data$alpha, ":::", "|", ":::",
                    .data$alpha, "$",
                    sep = ""
                )
            ) &
                .data$alpha != .data$special_aChain ~
                paste(.data$special_aChain, .data$beta, sep = "_"),
            TRUE ~ NA
        ))

    custom_a <- custom_a %>%
        unite(
            col = "CTaa_donor_source",
            all_of(c("CTaa", donor, "cloneCallSource")),
            sep = " ", remove = FALSE
        )
    custom_a <- custom_a %>%
        filter(!is.na(.data$new_clone)) %>%
        select(all_of(c("CTaa_donor_source", "new_clone")))

    obj@meta.data <- obj@meta.data %>%
        left_join(custom_a, by = "CTaa_donor_source")
    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = case_when(
            !is.na(.data$new_clone) ~ .data$new_clone,
            is.na(.data$new_clone) ~ .data$CTaa
        )) %>%
        select(-all_of("new_clone"))

    custom_b <- obj@meta.data %>%
        filter(!is.na(.data$CTaa) & !is.na(.data[[donor]]) &
            .data$cloneCallSource == "VDJdive") %>%
        group_by(.data$CTaa, .data[[donor]], .data$cloneCallSource) %>%
        summarise(count = n()) %>%
        arrange(desc(.data$count))
    custom_b <- custom_b %>%
        separate(.data$CTaa,
            into = c("alpha", "beta"),
            remove = FALSE, sep = "_"
        )
    custom_b <- custom_b %>%
        mutate(special_b = case_when(
            str_detect(.data$beta, ":::") & .data$cloneCallSource == "VDJdive"
            ~ .data$beta,
            !str_detect(.data$beta, ":::") | .data$cloneCallSource != "VDJdive"
            ~ NA
        ))
    custom_b <- custom_b %>%
        group_by(.data$alpha, .data[[donor]]) %>%
        arrange(desc(.data$count)) %>%
        mutate(special_bChain = first(.data$special_b)) %>%
        ungroup()
    custom_b <- custom_b %>%
        mutate(new_clone = case_when(
            str_detect(
                .data$special_bChain,
                paste("^", .data$beta, ":::", "|",
                    ":::", .data$beta, "$",
                    sep = ""
                )
            ) &
                .data$beta != .data$special_bChain ~
                paste(.data$alpha, .data$special_bChain, sep = "_"),
            TRUE ~ NA
        ))

    custom_b <- custom_b %>%
        unite(
            col = "CTaa_donor_source",
            all_of(c("CTaa", donor, "cloneCallSource")),
            sep = " ", remove = FALSE
        )
    custom_b <- custom_b %>%
        filter(!is.na(.data$new_clone)) %>%
        select(all_of(c("CTaa_donor_source", "new_clone")))

    obj@meta.data <- obj@meta.data %>%
        left_join(custom_b, by = "CTaa_donor_source")
    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = case_when(
            !is.na(.data$new_clone) ~ .data$new_clone,
            is.na(.data$new_clone) ~ .data$CTaa
        )) %>%
        select(-all_of(c("new_clone", "CTaa_donor_source")))

    ### Making columns so that scRepertoire visualizations can be used.
    ### Grouping by the group specified in the arguments, not the donor
    obj@meta.data <- obj@meta.data %>%
        group_by(.data[[group]]) %>%
        mutate(eclipseTcrGroupSize = sum(!is.na(.data$CTaa))) %>%
        ungroup()
    obj@meta.data <- obj@meta.data %>%
        group_by(.data$CTaa, .data[[group]]) %>%
        mutate(
            clonalProportion = case_when(
                !is.na(.data$CTaa) ~ n() / first(.data$eclipseTcrGroupSize),
                is.na(.data$CTaa) ~ NA
            ),
            clonalFrequency = case_when(
                !is.na(.data$CTaa) ~ n(),
                is.na(.data$CTaa) ~ NA
            )
        ) %>%
        ungroup() %>%
        select(-all_of("eclipseTcrGroupSize"))
    obj@meta.data <- obj@meta.data %>%
        mutate(cloneSize = factor(case_when(
            .data$clonalProportion <= 1 & .data$clonalProportion > 0.1
            ~ "Hyperexpanded (0.1 < X <= 1)",
            .data$clonalProportion <= 0.1 & .data$clonalProportion > 0.01
            ~ "Large (0.01 < X <= 0.1)",
            .data$clonalProportion <= 0.01 & .data$clonalProportion > 0.001
            ~ "Medium (0.001 < X <= 0.01)",
            .data$clonalProportion <= 0.001 & .data$clonalProportion > 1e-04
            ~ "Small (1e-04 < X <= 0.001)",
            .data$clonalProportion <= 1e-04 & .data$clonalProportion > 0
            ~ "Rare (0 < X <= 1e-04)",
            .data$clonalProportion <= 0
            ~ "None ( < X <= 0)",
            is.na(.data$clonalProportion)
            ~ NA
        )))

    ### Making a column that notes whether the cells in the clone appear
    ### to be MAIT or iNKT based on TCR segment usage
    obj@meta.data <- obj@meta.data %>%
        mutate(
            scR_mait_conv = case_when(
                str_detect(.data$scR_CTgene, "TRAV1-2\\.TRAJ33") ~
                    "MAIT Conventional (TRAV1-2 TRAJ33)",
                TRUE ~
                    "Normal"
            ),
            scR_mait_unconv = case_when(
                str_detect(
                    .data$scR_CTgene,
                    "TRAV1-2\\.TRAJ12|TRAV1-2\\.TRAJ20"
                )
                ~ "MAIT Unconventional (TRAV1-2 TRAJ12 or TRAJ20)",
                TRUE ~ "Normal"
            ),
            scR_inkt = case_when(
                str_detect(.data$scR_CTgene, "TRAV10\\.TRAJ18") &
                    str_detect(.data$scR_CTgene, "TRBV25-1")
                ~ "iNKT",
                TRUE ~ "Normal"
            )
        )
    obj@meta.data <- obj@meta.data %>%
        mutate(scR_unconventional_count = 3 -
            (.data$scR_mait_conv == "Normal") -
            (.data$scR_mait_unconv == "Normal") -
            (.data$scR_inkt == "Normal"))
    obj@meta.data <- obj@meta.data %>%
        mutate(cell_unconventional_subset = case_when(
            .data$scR_unconventional_count == 0 ~ "Conventional",
            .data$scR_unconventional_count == 1 &
                .data$scR_mait_conv != "Normal"
            ~ .data$scR_mait_conv,
            .data$scR_unconventional_count == 1 &
                .data$scR_mait_unconv != "Normal"
            ~ .data$scR_mait_unconv,
            .data$scR_unconventional_count == 1 & .data$scR_inkt != "Normal"
            ~ .data$scR_inkt,
            .data$scR_unconventional_count > 1
            ~ "Multiple Unconventional Types"
        ))

    df <- obj@meta.data %>%
        group_by(
            .data$CTaa,
            .data$cell_unconventional_subset,
            .data[[group]]
        ) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        filter(.data$cell_unconventional_subset != "Conventional") %>%
        arrange(desc(.data$count))

    if (nrow(df) > 0) {
        obj@meta.data <- obj@meta.data %>%
            ungroup() %>%
            mutate(CTaa_group = case_when(
                !is.na(.data$CTaa) & !is.na(.data[[group]])
                ~ paste(.data$CTaa, .data[[group]], sep = " "),
                TRUE ~ NA
            ))
        df <- df %>%
            group_by(.data$CTaa, .data[[group]]) %>%
            summarise(
                clone_unconventional_subset =
                    first(.data$cell_unconventional_subset)
            )
        df <- df %>%
            ungroup() %>%
            mutate(CTaa_group = case_when(
                !is.na(.data$CTaa) & !is.na(.data[[group]])
                ~ paste(.data$CTaa, .data[[group]], sep = " "),
                TRUE ~ NA
            ))
        df <- df %>%
            ungroup() %>%
            select(all_of(c("CTaa_group", "clone_unconventional_subset")))
        
        if ("clone_unconventional_subset" %in% colnames(obj@meta.data)) {
            obj@meta.data <- obj@meta.data %>%
                select(-all_of("clone_unconventional_subset"))
        }
        
        obj@meta.data <- obj@meta.data %>% left_join(df, by = "CTaa_group")
        obj@meta.data <- obj@meta.data %>%
            mutate(clone_unconventional_subset = case_when(
                is.na(.data$clone_unconventional_subset)
                ~ "Conventional",
                !is.na(.data$clone_unconventional_subset)
                ~ .data$clone_unconventional_subset
            )) %>%
            select(-all_of("CTaa_group"))
        if (verbose == TRUE) {
            message(
                sum(obj$clone_unconventional_subset != "Conventional"),
                " cells were detected in clones that appear to ",
                "be MAIT or iNKT cells"
            )
        }
    }

    obj@meta.data <- obj@meta.data %>% select(-all_of(c(
        "scR_mait_conv",
        "scR_mait_unconv",
        "scR_inkt",
        "scR_unconventional_count",
        "tcrEclipseGroup",
        "tcrEclipse_barcode"))
    )
    
    obj@meta.data <- as.data.frame(obj@meta.data)
    rownames(obj@meta.data) <- obj$barcode

    if (verbose == TRUE) {
        message(
            round(sum(obj$cloneCallSource == "VDJdive", na.rm = TRUE) /
                nrow(obj@meta.data) * 100, 1),
            "% of cells in the Seurat object were assigned a clonotype using ",
            "high confidence prediction from VDJdive::clonoStats()"
        )
        message(
            round(sum(!is.na(obj$CTaa)) / nrow(obj@meta.data) * 100, 1),
            "% of cells in the Seurat object were assigned a clonotype ",
            "in the end (either from prediction or traditional analysis)"
        )
        message("These final clonal annotations are stored in CTaa")
        message("Summary of the top 5 largest clones:")

        big_clones <- obj@meta.data %>%
            group_by(.data[[group]]) %>%
            mutate(group_count = sum(!is.na(.data$CTaa))) %>%
            ungroup() %>%
            group_by(.data$CTaa, .data[[group]]) %>%
            summarise(
                Clone_Size = n(),
                Percent_of_Repertoire = 100 * n() /
                    first(.data$group_count)
            ) %>%
            select(all_of(c(
                "CTaa",
                group,
                "Clone_Size",
                "Percent_of_Repertoire"
                ))
            ) %>%
            arrange(desc(.data$Clone_Size)) %>%
            filter(!is.na(.data$CTaa)) %>%
            head(n = 5)
        print(big_clones)
    }
    
    ### Cleaning up object before returning to user
    if ("originalBarcodeBeforeECLIPSE" %in% colnames(obj@meta.data)) {
        obj@meta.data <- obj@meta.data %>%
            mutate(barcode = .data$originalBarcodeBeforeECLIPSE) %>%
            select(-all_of("originalBarcodeBeforeECLIPSE"))
    } else if ("originalBarcodeBeforeECLIPSE" %in% colnames(obj@meta.data)
    == FALSE) {
        obj@meta.data <- obj@meta.data %>% select(-all_of("barcode"))
    }

    obj@meta.data <- obj@meta.data %>% select(-all_of("VDJdive_clone"))

    if (detailedEclipseStats == FALSE) {
        obj@meta.data <- obj@meta.data %>% select(-all_of(c(
            "EM_max_prop",
            "EM_2ndmax_prop",
            "EM_max_2ndmax_ratio"))
        )
    }

    return(list(obj, df))
}







clusterPermutationTestSimilarChains1v2 <- function(object, tcr, minCells, clusters, assignment_matrix){
              
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
              print(paste0("Observed n = ", nrow(diff)))
              mean_diff <- mean(diff$total)
              
              
              ### 2nd half, permutation now
              
              ### Extract values from EM assignment matrix
              p <- assignment_matrix@p
              j <- assignment_matrix@j
              x <- assignment_matrix@x
              
              all_clones <- colnames(assignment_matrix)
              n <- nrow(assignment_matrix)
              
              out <- list()
              
              for (i in seq_len(n)) {
                  
                  idx <- (p[i] + 1):p[i + 1]
                  
                  vals <- x[idx]
                  cols <- j[idx] + 1
                  clones <- all_clones[cols]
                  
                  clone_df <- data.frame("cell" = rownames(assignment_matrix)[i],
                                         "prediction" = clones,
                                         "props" = vals)
                  out[[i]] <- clone_df
              }
              
              out2 <- as.data.frame(do.call(rbind, out))
              
              ### Get metadata about the cells with the same chains as is predicted
              cells <- object@meta.data %>% filter(!is.na(scR_CTaa)) %>%
                                            select(Cell_barcode, scR_CTaa, cloneCallSource, CTaa)
              
              out2 <- out2 %>% left_join(cells,
                                         by = c("cell" = "Cell_barcode"),
                                         keep = FALSE)
              
              ### Filter out only clones that actually were predicted, and the clone must have a probability of less than <0.999 so that there will be multiple options available for comparison
              out3 <- out2 %>% filter(cloneCallSource == "VDJdive") %>%
                               group_by(cell) %>%
                               filter(max(props) <= 0.999) %>%
                               ungroup()
              
              out3 <- out3 %>% mutate(prediction = str_replace(prediction, " ", "_")) %>%
                               mutate(prediction = str_replace(prediction, ";", ":::"))
              
              ### Finding cases where two options actually were the same clone (if the clone has 3 total chains). These need to be reformatted
              out3 <- out3 %>% mutate(type = case_when(
                                                       !str_detect(CTaa, ":::") ~ "1A1B",
                                                       str_detect(CTaa, "^.*:::.*_") ~ "2A1B",
                                                       str_detect(CTaa, "^.*_.*:::") ~ "1A2B"))
              
              out3 <- out3 %>% mutate(potential_match = case_when(
                                                                  type == "2A1B" & ((str_detect(CTaa, paste0("_", str_replace(prediction, "^.*_", ""), "$")) &
                                                                      str_detect(CTaa, paste0("^", str_replace(prediction, "_.*$", ""), ":::|:::", str_replace(prediction, "_.*$", ""), "_"))) | prediction == CTaa) ~ "Yes",
                                                                  type == "1A2B" & ((str_detect(CTaa, paste0("^", str_replace(prediction, "_.*$", ""), "_")) &
                                                                      str_detect(CTaa, paste0("_", str_replace(prediction, "^.*_", ""), ":::|:::", str_replace(prediction, "^.*_", ""), "$"))) | prediction == CTaa) ~ "Yes",
                                                                  TRUE ~ "No"))
              
              ### Removing effective duplicates. If 2 cells have the same scR_CTaa, then the EM algorithm will predict the same clonotypes with the same frequencies for them
              out4 <- out3 %>% group_by(scR_CTaa, prediction, props) %>%
                               filter(row_number() == 1)
              
              out4 <- out4 %>% group_by(cell) %>%
                               mutate(predicted_chosen = max(props) == props) %>%
                               ungroup()
              
              ### Removing extremely unlikely prediction (these likely won't be any real observed clonotoype, just a predicted one)
              out4 <- out4 %>% filter(props >= 0.01) %>%
                               mutate(final_match = potential_match == "Yes" | predicted_chosen == TRUE)
              
              ### Making sure there are 2 classes so there can be a TVD calculation, then filtering out cells that likely are doublets
              out4 <- out4 %>% group_by(scR_CTaa) %>% 
                               filter(n_distinct(final_match) == 2) %>%
                               ungroup()
              
              out4 <- out4 %>% filter(str_count(scR_CTaa, ";") <= 2 & !str_detect(scR_CTaa, "_.*;.*;|;.*;.*_"))
              
              ### Calculating TVD for most likely vs. other potential clones for each observed chain pairing.
              clones <- unique(out4$scR_CTaa)
              
              meta <- object@meta.data
              
              max_clones <- out4 %>% filter(final_match == TRUE) %>%
                                     arrange(which(scR_CTaa %in% clones)) %>%
                                     select(scR_CTaa, prediction)
              
              
              summaries <- data.frame("scR_CTaa" = NA, "TVD" = NA)
              
              for (i in seq_along(clones)) {
                
                ### Filtering all cells to those with the proper scR_CTaa
                cluster_props <- meta %>% filter(scR_CTaa %in% out4$prediction[out4$scR_CTaa == clones[i]])
                
                cluster_props <- cluster_props %>% mutate(assigned = scR_CTaa %in% max_clones$prediction[max_clones$scR_CTaa == clones[i]])
                
                if (sum(cluster_props$assigned) >= minCells & sum(!cluster_props$assigned) >= minCells) {  
                    
                    ### Calculating cluster usage
                    cluster_props2 <- cluster_props %>% group_by(assigned) %>%
                                                        mutate(n_assigned = n()) %>%
                                                        group_by(.data[[clusters]], assigned) %>%
                                                        summarise(count = n(),
                                                                  clusterSubsetProp = n() / first(n_assigned),
                                                                  subsetSize = first(n_assigned)) %>%
                                                        ungroup()
                    
                    ### Completing so each class has 2 options, even if 1 option happened 0% of the time
                    cluster_props2 <- cluster_props2 %>% complete(.data[[clusters]], assigned)
                    cluster_props2$clusterSubsetProp[is.na(cluster_props2$clusterSubsetProp)] <- 0
                    
                    ### Getting TVD
                    diffs <- cluster_props2 %>% group_by(.data[[clusters]]) %>%
                                                summarise(diff = abs(clusterSubsetProp[assigned] - clusterSubsetProp[!assigned]))
                    
                    tvd <- sum(diffs$diff) * 0.5
                
                    summaries[i, 1] <- clones[i]
                    summaries[i, 2] <- tvd
                } else {
                    
                    summaries[i, 1] <- clones[i]
                    summaries[i, 2] <- NA
                }
                
              }
              
              summaries <- summaries %>% filter(!is.na(TVD))
              
              ### Randomly permuting mean TVD using the number of 3 chain clones
              tvd_means <- c()
              for (i in 1:10000){
              
                  rows <- sample(1:nrow(summaries), nrow(diff))
                  summaries_slim <- summaries[rows, ]
              
                  tvd_mn <- mean(summaries_slim$TVD)
                  tvd_means[i] <- tvd_mn
              }
              
              return(list(tvd_means, mean_diff))

}




clusterPermutationTestSimilarChainsMissing <- function(object, tcr, minCells, clusters, assignment_matrix){
                  
                  
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
                  print(paste0("Observed n = ", nrow(diff)))
                  mean_diff <- mean(diff$total)
                  
                  
                  ### 2nd half, permutations now
                  
                  ### Extract values from EM assignment matrix
                  p <- assignment_matrix@p
                  j <- assignment_matrix@j
                  x <- assignment_matrix@x
                  
                  all_clones <- colnames(assignment_matrix)
                  n <- nrow(assignment_matrix)
                  
                  out <- list()
                  
                  for (i in seq_len(n)) {
                      
                      idx <- (p[i] + 1):p[i + 1]
                      
                      vals <- x[idx]
                      cols <- j[idx] + 1
                      clones <- all_clones[cols]
                      
                      clone_df <- data.frame("cell" = rownames(assignment_matrix)[i],
                                             "prediction" = clones,
                                             "props" = vals)
                      out[[i]] <- clone_df
                  }
                  
                  out2 <- as.data.frame(do.call(rbind, out))
                  
                  ### Get metadata about the cells with the same chains as is predicted
                  cells <- object@meta.data %>% filter(!is.na(scR_CTaa)) %>%
                                                select(Cell_barcode, scR_CTaa, cloneCallSource, CTaa)
                  
                  out2 <- out2 %>% left_join(cells,
                                             by = c("cell" = "Cell_barcode"),
                                             keep = FALSE)
                  
                  ### Filter out only clones that actually were predicted, and the clone must have a probability of less than <0.999 so that there will be multiple options available for comparison
                  out3 <- out2 %>% filter(cloneCallSource == "VDJdive") %>%
                                   group_by(cell) %>%
                                   filter(max(props) <= 0.999) %>%
                                   ungroup()
                  
                  out3 <- out3 %>% mutate(prediction = str_replace(prediction, " ", "_")) %>%
                                   mutate(prediction = str_replace(prediction, ";", ":::"))
                  
                  ### Finding cases where two options actually were the same clone (if the clone has 3 total chains). These need to be reformatted
                  out3 <- out3 %>% mutate(type = case_when(
                                                           !str_detect(CTaa, ":::") ~ "1A1B",
                                                           str_detect(CTaa, "^.*:::.*_") ~ "2A1B",
                                                           str_detect(CTaa, "^.*_.*:::") ~ "1A2B"))
                  
                  out3 <- out3 %>% mutate(potential_match = case_when(
                                                                      type == "2A1B" & ((str_detect(CTaa, paste0("_", str_replace(prediction, "^.*_", ""), "$")) &
                                                                          str_detect(CTaa, paste0("^", str_replace(prediction, "_.*$", ""), ":::|:::", str_replace(prediction, "_.*$", ""), "_"))) | prediction == CTaa) ~ "Yes",
                                                                      type == "1A2B" & ((str_detect(CTaa, paste0("^", str_replace(prediction, "_.*$", ""), "_")) &
                                                                          str_detect(CTaa, paste0("_", str_replace(prediction, "^.*_", ""), ":::|:::", str_replace(prediction, "^.*_", ""), "$"))) | prediction == CTaa) ~ "Yes",
                                                                      TRUE ~ "No"))
                  
                  ### Removing effective duplicates. If 2 cells have the same scR_CTaa, then the EM algorithm will predict the same clonotypes with the same frequencies for them
                  out4 <- out3 %>% group_by(scR_CTaa, prediction, props) %>%
                                   filter(row_number() == 1)
                  
                  out4 <- out4 %>% group_by(cell) %>%
                                   mutate(predicted_chosen = max(props) == props) %>%
                                   ungroup()
                  
                  ### Removing extremely unlikely prediction (these likely won't be any real observed clonotoype, just a predicted one)
                  out4 <- out4 %>% filter(props >= 0.01) %>%
                                   mutate(final_match = potential_match == "Yes" | predicted_chosen == TRUE)
                  
                  ### Making sure there are 2 classes so there can be a TVD calculation, then filtering out cells that likely are doublets
                  out4 <- out4 %>% group_by(scR_CTaa) %>% 
                                   filter(n_distinct(final_match) == 2) %>%
                                   ungroup()
                  
                  out4 <- out4 %>% filter(str_count(scR_CTaa, ";") <= 2 & !str_detect(scR_CTaa, "_.*;.*;|;.*;.*_"))
                  
                  ### Calculating TVD for most likely vs. other potential clones for each observed chain pairing.
                  clones <- unique(out4$scR_CTaa)
                  
                  meta <- object@meta.data
                  
                  max_clones <- out4 %>% filter(final_match == TRUE) %>%
                                         arrange(which(scR_CTaa %in% clones)) %>%
                                         select(scR_CTaa, prediction)
                  
                  
                  summaries <- data.frame("scR_CTaa" = NA, "TVD" = NA)
                  
                  for (i in seq_along(clones)) {
                    
                    ### Filtering all cells to those with the proper scR_CTaa
                    cluster_props <- meta %>% filter(scR_CTaa %in% out4$prediction[out4$scR_CTaa == clones[i]])
                    
                    cluster_props <- cluster_props %>% mutate(assigned = scR_CTaa %in% max_clones$prediction[max_clones$scR_CTaa == clones[i]])
                    
                    if (sum(cluster_props$assigned) >= minCells & sum(!cluster_props$assigned) >= minCells) {  
                        
                        ### Calculating cluster usage
                        cluster_props2 <- cluster_props %>% group_by(assigned) %>%
                                                            mutate(n_assigned = n()) %>%
                                                            group_by(.data[[clusters]], assigned) %>%
                                                            summarise(count = n(),
                                                                      clusterSubsetProp = n() / first(n_assigned),
                                                                      subsetSize = first(n_assigned)) %>%
                                                            ungroup()
                        
                        ### Completing so each class has 2 options, even if 1 option happened 0% of the time
                        cluster_props2 <- cluster_props2 %>% complete(.data[[clusters]], assigned)
                        cluster_props2$clusterSubsetProp[is.na(cluster_props2$clusterSubsetProp)] <- 0
                        
                        ### Getting TVD
                        diffs <- cluster_props2 %>% group_by(.data[[clusters]]) %>%
                                                    summarise(diff = abs(clusterSubsetProp[assigned] - clusterSubsetProp[!assigned]))
                        
                        tvd <- sum(diffs$diff) * 0.5
                    
                        summaries[i, 1] <- clones[i]
                        summaries[i, 2] <- tvd
                    } else {
                        
                        summaries[i, 1] <- clones[i]
                        summaries[i, 2] <- NA
                    }
                    
                  }
                  
                  summaries <- summaries %>% filter(!is.na(TVD))
                  
                  ### Randomly permuting mean TVD using the number of 3 chain clones
                  tvd_means <- c()
                  for (i in 1:10000){
                      
                      ### Note this has to be replace = TRUE because nrow(diff) is larger than nrow(summaries)
                      rows <- sample(1:nrow(summaries), nrow(diff), replace = TRUE)
                      summaries_slim <- summaries[rows, ]
                  
                      tvd_mn <- mean(summaries_slim$TVD)
                      tvd_means[i] <- tvd_mn
                  }
                  
                  return(list(tvd_means, mean_diff))

}

countChains <- function(object) {
    ### Makes column in meta data that has the number of alpha/beta
    ### chains both originally and after ECLIPSE.
    object@meta.data <- object@meta.data %>%
        separate(.data$highConfChains,
            remove = FALSE,
            sep = "_", into = c("scr_alpha", "scr_beta")
        )
    object@meta.data <- object@meta.data %>%
        separate(.data$CTaa,
            remove = FALSE, sep = "_",
            into = c("ctaa_alpha", "ctaa_beta")
        )

    object@meta.data <- object@meta.data %>%
        mutate(
            scr_alphaCount = case_when(
                .data$scr_alpha == "NA" | is.na(.data$scr_alpha) ~ 0,
                .data$scr_alpha != "NA" ~ str_count(.data$scr_alpha, ";") + 1
            ),
            scr_betaCount = case_when(
                .data$scr_beta == "NA" | is.na(.data$scr_beta) ~ 0,
                .data$scr_beta != "NA" ~ str_count(.data$scr_beta, ";") + 1
            )
        )


    ### This adjusts for if the clones are called by my pipeline as supposed
    ### to be having 2 of a chain.
    ### Still if something has 2 that doesn't mean it's 100% a doublet
    ### because not everything can be called by my pipeline
    object@meta.data <- object@meta.data %>%
        mutate(
            adj_scr_alphaCount = case_when(
                .data$scr_alpha == str_replace(.data$ctaa_alpha, ":::", ";") &
                    str_detect(.data$ctaa_alpha, ":::") &
                    .data$cloneCallSource == "VDJdive"
                ~ .data$scr_alphaCount - 1,
                TRUE ~ .data$scr_alphaCount
            ),
            adj_scr_betaCount = case_when(
                .data$scr_beta == str_replace(.data$ctaa_beta, ":::", ";") &
                    str_detect(.data$ctaa_beta, ":::") &
                    .data$cloneCallSource == "VDJdive"
                ~ .data$scr_betaCount - 1,
                TRUE ~ .data$scr_betaCount
            )
        )


    ### Finds the total amount of chains, both raw and adjusted for
    ### clones called as being 2 of a chain in my pipeline
    object@meta.data <- object@meta.data %>%
        mutate(
            scr_count = .data$scr_alphaCount + .data$scr_betaCount,
            adj_scr_count = .data$adj_scr_alphaCount + .data$adj_scr_betaCount
        )
    
    return(object)
    
}





pre2 <- function(seurat_object,
                    contig_type = "all",
                    contig_folder_paths = NULL,
                    contig_file_paths = NULL,
                    batch,
                    donor,
                    group,
                    original_barcode,
                    format = "Blank",
                    detailedEclipseStats = FALSE,
                    combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
                    allowDoubleChains = "both",
                    returnFilteredContigs = FALSE,
                    verbose = TRUE) {
  
    #### THIS FIRST PART IS READING THE DATA IN AND DOING QC
    ### Checking the input data is usable
    if (!contig_type %in% c(
        "ALL", "All", "all",
        "Filtered", "FILTERED", "filtered",
        "all_contig_annotations.csv",
        "filtered_contig_annotations.csv",
        "manual"
    )) {
        stop(
            "Incorrect contig type. Should be all or filtered. ",
            "Returning nothing"
        )
    }

    if (is.null(contig_folder_paths) & contig_type != "manual") {
        stop(
            "Please provide folders to contig_folder_paths that have the ",
            "contig annotations in them."
        )
    }

    if (is.null(contig_file_paths) & contig_type == "manual") {
        stop(
            "No contig file paths provided. ",
            "Please provide a vector of file path names to contig_file_paths"
        )
    }

    if (contig_type != "manual" &
        (!is.vector(contig_folder_paths) | is.list(contig_folder_paths) |
        length(contig_folder_paths) == 0)) {
        stop("contig_folder_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (contig_type == "manual" &
        (!is.vector(contig_file_paths) | is.list(contig_file_paths) |
        length(contig_file_paths) == 0)) {
        stop("contig_file_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (!is.null(contig_file_paths) & !is.null(contig_folder_paths)) {
        stop("Please do not use contig_file_paths and contig_folder_paths. ",
            "Please try again.")
    }

    if (!batch %in% colnames(seurat_object@meta.data)) {
        stop(
            "Batch column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!donor %in% colnames(seurat_object@meta.data)) {
        stop(
            "Donor column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!group %in% colnames(seurat_object@meta.data)) {
        stop(
            "Group column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!original_barcode %in% colnames(seurat_object@meta.data)) {
        stop(
            "Original barcode column name provided is not ",
            "in the Seurat object. ",
            "Please try again."
        )
    }

    if (!is.logical(detailedEclipseStats)) {
        stop(
            "TRUE or FALSE not provided to detailedEclipseStats. ",
            "Please try again."
        )
    }

    if (format %in% c("None", "none", "Blank") == FALSE) {
        stop(
            "None, none, or Blank not provided to format. ",
            "Please try again."
        )
    }

    if (!is.list(combineTcrArgs) | length(combineTcrArgs) != 4 |
        !is.logical(combineTcrArgs[[1]]) | !is.logical(combineTcrArgs[[2]]) |
        !is.logical(combineTcrArgs[[3]]) | !is.logical(combineTcrArgs[[4]])) {
        stop(
            "Incorrect input to combineTcrArgs. ",
            "Please provide a list of 4 logical (TRUE or FALSE) elements."
        )
    }

    if (allowDoubleChains %in% c(
        "Alpha", "alpha",
        "beta", "Beta",
        "Both", "both",
        "Neither", "neither"
    ) == FALSE) {
        stop(
            "Alpha, Beta, Both, or Neither not provided to ",
            "allowDoubleChains. ",
            "Please try again."
        )
    }

    if (!is.logical(verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )
    }

    if (!is.logical(returnFilteredContigs)) {
        stop(
            "TRUE or FALSE not provided to returnFilteredContigs. ",
            "Please try again."
        )
    }


    ### Reading the data in
    if (verbose == TRUE) message("\nReading TCR contig data in...\n")

    if (contig_type != "manual") {
        contigs_list <- vector("list", length(contig_folder_paths))

        for (i in seq_along(contig_folder_paths)) {
            if (contig_type %in% c(
                "ALL", "All", "all",
                "all_contig_annotations.csv"
            )) {
                files <- list.files(contig_folder_paths[i],
                    pattern = "all_contig_annotations.csv",
                    full.names = TRUE
                )
            } else if (contig_type %in%
                c(
                    "Filtered", "filtered", "FILTERED",
                    "filtered_contig_annotations.csv"
                )) {
                files <- list.files(contig_folder_paths[i],
                    pattern = "filtered_contig_annotations.csv",
                    full.names = TRUE
                )
            }


            if (length(files) == 0) {
                stop("No contig files found, returning nothing")
            } else if (length(files) > 1) {
                stop("Multiple contig files found, returning nothing")
            } else if (length(files) == 1) {
                contigs_list[[i]] <- read.csv(files[1])
                contigs_list[[i]]$barcode <- paste(i, "_",
                    contigs_list[[i]]$barcode,
                    sep = ""
                )
            }
        }
    } else if (contig_type == "manual") {
        contigs_list <- vector("list", length(contig_file_paths))

        for (i in seq_along(contig_file_paths)) {
            contigs_list[[i]] <- read.csv(contig_file_paths[i])

            if (nrow(contigs_list[[i]]) == 0) {
                stop(
                    "File path ", i, " did not lead to any ",
                    "contigs being loaded in. ",
                    "Check the path name and that it is a .csv file"
                )
            }

            contigs_list[[i]]$barcode <-
                paste(i, "_", contigs_list[[i]]$barcode, sep = "")
        }
    }

    ### Joins together all batches and adds data to them
    ### about which sample/donor/batch each contig was from
    contigs_full <- bind_rows(contigs_list)
    
    if ("sample" %in% colnames(contigs_full)) {
        contigs_full <- contigs_full %>% select(-all_of("sample"))
    }

    if (format %in% c("None", "none")) {
        contigs_full[contigs_full == "None"] <- ""
    }

    seurat_object$tcrEclipse_barcode <-
        paste(seurat_object@meta.data[[batch]],
            "_",
            seurat_object@meta.data[[original_barcode]],
            sep = ""
        )
    seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)

    cell_data <- seurat_object@meta.data %>%
        select(all_of(c(
            "tcrEclipse_barcode",
            unique(c(batch, donor, group)),
            "seurat_object_barcode"))
        )
    cell_data <- cell_data %>% filter(!is.na(.data[[batch]]))

    if (nrow(cell_data) != n_distinct(cell_data$tcrEclipse_barcode)) {
        stop(
            "There are likely multiple sequencing runs ",
            "labeled as one batch. ",
            "Please relabel the batch column such that each batch ",
            "is a unique sequncing run. ",
            "Returning nothing"
        )
    }

    all <- contigs_full %>%
        left_join(cell_data,
            by = c("barcode" = "tcrEclipse_barcode"),
            keep = FALSE
        )
    all <- all %>%
        mutate(barcode = .data$seurat_object_barcode) %>%
        select(-all_of("seurat_object_barcode"))
    if ("barcode" %in% colnames(seurat_object@meta.data)) {
        seurat_object@meta.data <- seurat_object@meta.data %>%
            mutate(originalBarcodeBeforeECLIPSE = .data$barcode) %>%
            select(-all_of("barcode"))
    }

    seurat_object@meta.data <- seurat_object@meta.data %>%
        dplyr::rename("barcode" = all_of("seurat_object_barcode"))

    if (verbose == TRUE) {
        message("QC Filtering:")
        message(" - ", nrow(all), " initial contigs found before filtering")
    }
    
    ### Renaming barcodes, filtering to only cells in seurat object,
    ### removing cells without a CDR3,
    ### and renaming QC metrics from cellranger so
    ### these contigs aren't thrown out
    all <- all %>% filter(.data$barcode %in% seurat_object$barcode)
    if (verbose == TRUE) {
        message(" - ", nrow(all), " contigs were from true cells ",
            "(coming from ", n_distinct(all$barcode), " cells)")
    }
    if (nrow(all) == 0) {
        stop(
            "Check that your barcodes are the original ",
            "(i.e. likely end with -1 and match the contig files)"
        )
    }

    if (verbose == TRUE) {
        message(" - ",
            sum(all$high_confidence %in% c("true", "True", "TRUE") == FALSE),
            " contigs were removed for not being high confidence as ",
            "defined by Cell Ranger"
        )
    }
    all <- all %>%
        filter(.data$high_confidence %in% c("true", "True", "TRUE"))
    all <- all %>% mutate(
        productive = "true",
        is_cell = "true",
        full_length = "true"
    )

    count_before_custom_filtering <- nrow(all)

    all <- all %>% filter(!.data$cdr3 %in% c("", "None") &
                            !.data$cdr3_nt %in% c("", "None"))
    all <- all %>% filter(.data$chain %in% c("TRA", "TRB", "TRG", "TRD"))

    ### Removing cells with abnormally long CDR3s
    all <- all %>% mutate(cdr3_length = str_length(.data$cdr3))

    all <- all %>% filter(.data$cdr3_length <= 30)

    ### There are TCRs with TRDV and then TRAJ and TRAC.
    ### Cellranger calls them as TRD,
    ### but they are actually TRA since they pair with a beta chain and
    ### have no J segment.
    all <- all %>%
        mutate(ad_hybrid_tcr = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC"))
            ~ TRUE,
            TRUE ~ FALSE
        ))

    all <- all %>%
        mutate(chain = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC")) ~ "TRA",
            TRUE ~ .data$chain
        ))

    all <- all %>% filter(.data$chain %in% c("TRA", "TRB"))

    ### Removing cells with a * at the end of their cdr3.
    ### These indicate stop codons, so removing these contigs
    all <- all %>% filter(!str_detect(.data$cdr3, "\\*"))

    ### Removing duplicated chains. 
    ### Some datasets (particularly contig files made with
    ### older version of Cellranger)
    ### have lots of duplicated rows that have identical CDR3
    ### and gene segments for the same barcode.
    ### This messes up future analysis
    all <- all %>%
        arrange(desc(.data$umis), desc(.data$reads)) %>%
        group_by(.data$barcode, .data$chain, .data$cdr3) %>%
        filter(row_number() == 1) %>%
        ungroup()
    
    return(list(all, seurat_object))
    
}








post2 <- function(seurat_object,
                 all,
                 contig_type = "all",
                 contig_folder_paths = NULL,
                 contig_file_paths = NULL,
                 batch,
                 donor,
                 group,
                 original_barcode,
                 format = "Blank",
                 detailedEclipseStats = FALSE,
                 combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
                 allowDoubleChains = "both",
                 returnFilteredContigs = FALSE,
                 verbose = TRUE,
                 thresholds = c(0.8, 0.7, 0.5)) {

    ### Noting cells with only 1 type of chain (alpha or beta) and
    ### that chain isn't seen in any other cell.
    ### We don't want these cells to be predicted.
    all2 <- all %>%
        group_by(.data$barcode) %>%
        mutate(chainTypes = n_distinct(.data$chain)) %>%
        ungroup()
    all2 <- all2 %>%
        group_by(.data$cdr3, .data$chain, .data[[donor]]) %>%
        mutate(maxChains = max(.data$chainTypes)) %>%
        ungroup()
    if (verbose == TRUE) {
        message(" - ",
            sum(all2$maxChains == 1),
            " orphan contigs were found with no observed pairing. ",
            "Cells containing these contigs are unlikely to have ",
            "a clonotype predicted"
        )
    }
    orphanContigs <- all2 %>%
        filter(.data$maxChains == 1) %>%
        select(all_of(c("barcode", "chain", "cdr3")))

    if (returnFilteredContigs == TRUE) {

        all2$orphanContig <- FALSE
        all2$orphanContig[all2$maxChains == 1] <- TRUE

        all2 <- all2 %>%
            select(-all_of(c(
                "cdr3_length",
                "ad_hybrid_tcr",
                "chainTypes",
                "maxChains"))
                )
        return(all2)
    }

    if (verbose == TRUE) message("\nFinding clones with extra chains...\n")

    donors <- unique(as.vector(all[, donor])[[1]])
    new_contigs <- vector("list", length(donors))

    for (i in seq_along(donors)) {
        ### THIS PART IS FINDING THE COMBO CLONES AND
        ### REWRITING THE CONTIG FILES
        if (allowDoubleChains %in% c("Neither", "neither") == FALSE &
            verbose == TRUE) {
            message(
                "Running donor ", i, "/",
                length(donors), " (", donors[i], "):"
            )
        }

        slim <- all2 %>% filter(.data[[donor]] == donors[i])

        if (allowDoubleChains %in% c("Alpha", "alpha", "Both", "both")) {
            ### Calling clones with multiple alpha chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "alpha",
                printStats = verbose
            )
        }

        if (allowDoubleChains %in% c("Beta", "beta", "Both", "both")) {
            ### Calling clones with multiple beta chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "beta",
                printStats = verbose
            )
        }

        new_contigs[[i]] <- slim

        if (verbose == TRUE) message()
    }

    contigs_post <- bind_rows(new_contigs)

    ### THIS PART IS FEEDING THE EDITED CONTIG FILES INTO VDJdive
    ### WHICH CALLS AMBIGUOUS CELLS

    if (verbose == TRUE) {
        message("Running VDJdive::clonoStats() for clone EM assignment...")
    }
    ### Making a temporary directory
    ### Writing a .csv file as filtered_contig_annotations.csv
    ### Reading that into VDJdive
    ### Using clonoStats for EM assignment
    ### Extracting output and changing format

    temp_dir <- tempfile("ECLIPSE_")
    dir.create(temp_dir)

    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

    write.csv(contigs_post,
        paste(temp_dir, "/filtered_contig_annotations.csv", sep = ""),
        row.names = FALSE
    )

    contigs <- readVDJcontigs(temp_dir)
    
    ### Please note suppressMessages is used here because
    ### there is a message stating internal coercion from
    ### Tsparse* to CsparseMatrix
    ### that I do not want to be printed.
    vdj <- suppressMessages(clonoStats(contigs,
        method = "EM",
        type = "TCR",
        assignment = TRUE,
        group = donor
    ))

    if (verbose == TRUE) message("Extracting results from VDJdive...")
    df <- vdj@assignment
    df <- as(df, "RsparseMatrix")
    colnames(df) <- clonoNames(vdj)

    rm(vdj)

    ### Finding the most likely clonotype, its probability,
    ### and the probability of the 2nd most likely clonotype for each cell.
    p <- df@p
    j <- df@j
    x <- df@x

    all_clones <- colnames(df)
    n <- nrow(df)

    barcode <- rownames(df)
    clone <- character(n)
    maxes <- numeric(n)
    second_max <- numeric(n)


    for (i in seq_len(n)) {
        idx <- (p[i] + 1):p[i + 1]

        vals <- x[idx]
        cols <- j[idx] + 1

        if (length(vals) == 1) {
            maxes[i] <- vals
            second_max[i] <- 0
            clone[i] <- all_clones[cols]
        } else {
            best <- vals[1]
            second <- 0
            best_col <- cols[1]

            for (k in 2:length(vals)) {
                if (vals[k] > best) {
                    second <- best
                    best <- vals[k]
                    best_col <- cols[k]
                } else if (vals[k] > second) {
                    second <- vals[k]
                }
            }

            maxes[i] <- best
            second_max[i] <- second
            clone[i] <- all_clones[best_col]
        }
    }


    ### Summarizing the VDJdive EM data per cell
    ### and then filtering cells that don't meet QC thresholds
    em <- data.frame(barcode = rownames(df), maxes, clone, second_max)
    rm(df, p, j, x)

    if (verbose == TRUE) {
        message("Finalizing clonal assignments...\n")
    }
    em$ratio <- em$maxes / em$second_max

    if (verbose == TRUE) {
        em2 <- em %>%
            arrange(.data$maxes) %>%
            mutate(
                nm = row_number(),
                assigned_by_VDJdive = .data$maxes >= thresholds[1] |
                    (.data$maxes >= thresholds[2] & .data$ratio > 5) |
                    (.data$maxes > thresholds[3] & .data$ratio > 10)
            )
        plot2 <- em2 %>%
            ggplot() +
            geom_point(aes(
                x = .data$nm, y = .data$maxes,
                color = .data$assigned_by_VDJdive
            )) +
            xlab("Ranked position of cell") +
            ylab(paste("Probability of cell being\n",
                "in most likely clone", sep = "")) +
            ggtitle("Clonal Prediction QC") +
            theme_bw() +
            labs(color = "Source\nof Final\nClonotype") +
            scale_color_manual(values = c("cyan4", "darkorange2"),
                labels = c("Traditional\nAnalysis\n",
                    "High\nConfidence\nPrediction")) +
            theme(
                plot.title = element_text(
                    face = "bold", hjust = 0.5,
                    size = 20
                ),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 16, hjust = 0.5),
                legend.text = element_text(size = 12),
                legend.position.inside = c(0.7, 0.4)
            )
        plot(plot2)
    }

    final <- em %>% filter(.data$maxes >= thresholds[1] |
        (.data$maxes >= thresholds[2] & .data$ratio > 5) |
        (.data$maxes > thresholds[3] & .data$ratio > 10))

    ### Removing cells where an orphan TCR was assigned.
    ### There's no way for the algorithm to accurately predict these so
    ### we don't want them.
    orphanContigs <- orphanContigs %>% left_join(final, by = "barcode")
    orphanContigs <- orphanContigs %>%
        mutate(remove = case_when(
            is.na(.data$clone) ~ FALSE,
            .data$chain == "TRA" &
                str_detect(
                    .data$clone,
                    paste("^", .data$cdr3, " ", sep = "")
                ) ~ TRUE,
            .data$chain == "TRB" &
                str_detect(
                    .data$clone,
                    paste(" ", .data$cdr3, "$", sep = "")
                ) ~ TRUE,
            TRUE ~ FALSE
        ))
    orphanContigs <- orphanContigs %>%
        filter(.data$remove) %>%
        select(all_of(c("barcode", "remove")))
    orphanContigs <- unique(orphanContigs)

    final <- final %>% left_join(orphanContigs, by = "barcode")
    final <- final %>% filter(is.na(.data$remove))
    
    ### Adding the VDJdive and scRepertoire calls and
    ### expressed chains to the Seurat object
    seurat_object$tcrEclipseGroup <- seurat_object@meta.data[[group]]
    
    contig_list <- createHTOContigList(all,
        seurat_object,
        group.by = "tcrEclipseGroup"
    )
    
    no_vdjd <- combineTCR(contig_list)
    
    seurat_object <- combineExpression(no_vdjd,
        seurat_object,
        cloneCall = "aa",
        proportion = TRUE
    )
    
    seurat_object$highConfChains <- seurat_object$CTaa
    seurat_object@meta.data <- seurat_object@meta.data %>%
        select(-all_of(c(
            "CTgene",
            "CTnt",
            "CTaa",
            "CTstrict",
            "clonalProportion",
            "clonalFrequency",
            "cloneSize"))
        )
    
    no_vdjd <- combineTCR(contig_list,
        removeNA = combineTcrArgs[[1]],
        removeMulti = combineTcrArgs[[2]],
        filterMulti = combineTcrArgs[[3]],
        filterNonproductive = combineTcrArgs[[4]]
    )
    obj <- combineExpression(no_vdjd,
        seurat_object,
        cloneCall = "aa",
        proportion = TRUE
    )
    rm(seurat_object)
    
    columns_to_remove <- which(colnames(obj@meta.data) %in% 
        c("scR_CTgene", "scR_CTnt", "scR_CTaa", "scR_CTstrict", 
            "scR_clonalFrequency", "scR_clonalProportion", "scR_cloneSize"))
    
    if (length(columns_to_remove) > 0) {
      
        obj@meta.data <- obj@meta.data[, -columns_to_remove]       
    
    }
    
    obj@meta.data <- obj@meta.data %>%
        rename(
            scR_CTgene = all_of("CTgene"),
            scR_CTnt = all_of("CTnt"),
            scR_CTaa = all_of("CTaa"),
            scR_CTstrict = all_of("CTstrict"),
            scR_clonalProportion = all_of("clonalProportion"),
            scR_clonalFrequency = all_of("clonalFrequency"),
            scR_cloneSize = all_of("cloneSize")
        )
    final <- final %>% mutate(clone2 = str_replace(.data$clone, " ", "_"))
    final <- final %>% mutate(clone2 = str_replace_all(
        .data$clone2,
        ";", ":::"
    ))
    final <- final %>% select(all_of("barcode"),
        "VDJdive_clone" = all_of("clone2"),
        "EM_max_prop" = all_of("maxes"),
        "EM_2ndmax_prop" = all_of("second_max"),
        "EM_max_2ndmax_ratio" = all_of("ratio")
    )
    
    obj@meta.data <- obj@meta.data %>% left_join(final, by = "barcode")
    

    ### Making CTaa our final clone call which uses VDJdive
    ### if available but if not the scRepertoire call
    obj@meta.data <- obj@meta.data %>%
        mutate(cloneCallSource = case_when(
            !is.na(.data$VDJdive_clone) ~ "VDJdive",
            is.na(.data$VDJdive_clone) & !is.na(.data$scR_CTaa) ~
                "scRepertoire",
            TRUE ~ NA
        ))

    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = case_when(
            !is.na(.data$VDJdive_clone) ~ .data$VDJdive_clone,
            is.na(.data$VDJdive_clone) & !is.na(.data$scR_CTaa) ~
                .data$scR_CTaa,
            TRUE ~ NA
        ))


    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = str_replace_all(.data$CTaa, ";", ":::"))
    
    ### Combining clones. Cells that were ambiguous may be assigned to a
    ### combo clone by VDJdive but only be listed with one chain
    ### Here those cells have their CDR3 renamed as the combo CDR3
    obj@meta.data <- obj@meta.data %>%
        ungroup() %>%
        mutate(CTaa_donor_source = case_when(
            !is.na(.data$CTaa) &
                !is.na(.data[[donor]]) &
                !is.na(.data$cloneCallSource) ~
                paste(.data$CTaa, .data[[donor]],
                    .data$cloneCallSource,
                    sep = " "
                ),
            TRUE ~ NA
        ))

    custom_a <- obj@meta.data %>%
        filter(!is.na(.data$CTaa) & !is.na(.data[[donor]]) &
            .data$cloneCallSource == "VDJdive") %>%
        group_by(.data$CTaa, .data[[donor]], .data$cloneCallSource) %>%
        summarise(count = n()) %>%
        arrange(desc(.data$count))
    custom_a <- custom_a %>%
        separate(.data$CTaa,
            into = c("alpha", "beta"),
            remove = FALSE, sep = "_"
        )
    custom_a <- custom_a %>%
        mutate(special_a = case_when(
            str_detect(.data$alpha, ":::") &
                .data$cloneCallSource == "VDJdive" ~ .data$alpha,
            !str_detect(.data$alpha, ":::") |
                .data$cloneCallSource != "VDJdive"
            ~ NA
        ))

    custom_a <- custom_a %>%
        group_by(.data$beta, .data[[donor]]) %>%
        arrange(desc(.data$count)) %>%
        mutate(special_aChain = first(.data$special_a)) %>%
        ungroup()
    custom_a <- custom_a %>%
        mutate(new_clone = case_when(
            str_detect(
                .data$special_aChain,
                paste("^", .data$alpha, ":::", "|", ":::",
                    .data$alpha, "$",
                    sep = ""
                )
            ) &
                .data$alpha != .data$special_aChain ~
                paste(.data$special_aChain, .data$beta, sep = "_"),
            TRUE ~ NA
        ))

    custom_a <- custom_a %>%
        unite(
            col = "CTaa_donor_source",
            all_of(c("CTaa", donor, "cloneCallSource")),
            sep = " ", remove = FALSE
        )
    custom_a <- custom_a %>%
        filter(!is.na(.data$new_clone)) %>%
        select(all_of(c("CTaa_donor_source", "new_clone")))

    obj@meta.data <- obj@meta.data %>%
        left_join(custom_a, by = "CTaa_donor_source")
    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = case_when(
            !is.na(.data$new_clone) ~ .data$new_clone,
            is.na(.data$new_clone) ~ .data$CTaa
        )) %>%
        select(-all_of("new_clone"))

    custom_b <- obj@meta.data %>%
        filter(!is.na(.data$CTaa) & !is.na(.data[[donor]]) &
            .data$cloneCallSource == "VDJdive") %>%
        group_by(.data$CTaa, .data[[donor]], .data$cloneCallSource) %>%
        summarise(count = n()) %>%
        arrange(desc(.data$count))
    custom_b <- custom_b %>%
        separate(.data$CTaa,
            into = c("alpha", "beta"),
            remove = FALSE, sep = "_"
        )
    custom_b <- custom_b %>%
        mutate(special_b = case_when(
            str_detect(.data$beta, ":::") & .data$cloneCallSource == "VDJdive"
            ~ .data$beta,
            !str_detect(.data$beta, ":::") | .data$cloneCallSource != "VDJdive"
            ~ NA
        ))
    custom_b <- custom_b %>%
        group_by(.data$alpha, .data[[donor]]) %>%
        arrange(desc(.data$count)) %>%
        mutate(special_bChain = first(.data$special_b)) %>%
        ungroup()
    custom_b <- custom_b %>%
        mutate(new_clone = case_when(
            str_detect(
                .data$special_bChain,
                paste("^", .data$beta, ":::", "|",
                    ":::", .data$beta, "$",
                    sep = ""
                )
            ) &
                .data$beta != .data$special_bChain ~
                paste(.data$alpha, .data$special_bChain, sep = "_"),
            TRUE ~ NA
        ))

    custom_b <- custom_b %>%
        unite(
            col = "CTaa_donor_source",
            all_of(c("CTaa", donor, "cloneCallSource")),
            sep = " ", remove = FALSE
        )
    custom_b <- custom_b %>%
        filter(!is.na(.data$new_clone)) %>%
        select(all_of(c("CTaa_donor_source", "new_clone")))

    obj@meta.data <- obj@meta.data %>%
        left_join(custom_b, by = "CTaa_donor_source")
    obj@meta.data <- obj@meta.data %>%
        mutate(CTaa = case_when(
            !is.na(.data$new_clone) ~ .data$new_clone,
            is.na(.data$new_clone) ~ .data$CTaa
        )) %>%
        select(-all_of(c("new_clone", "CTaa_donor_source")))

    ### Making columns so that scRepertoire visualizations can be used.
    ### Grouping by the group specified in the arguments, not the donor
    obj@meta.data <- obj@meta.data %>%
        group_by(.data[[group]]) %>%
        mutate(eclipseTcrGroupSize = sum(!is.na(.data$CTaa))) %>%
        ungroup()
    obj@meta.data <- obj@meta.data %>%
        group_by(.data$CTaa, .data[[group]]) %>%
        mutate(
            clonalProportion = case_when(
                !is.na(.data$CTaa) ~ n() / first(.data$eclipseTcrGroupSize),
                is.na(.data$CTaa) ~ NA
            ),
            clonalFrequency = case_when(
                !is.na(.data$CTaa) ~ n(),
                is.na(.data$CTaa) ~ NA
            )
        ) %>%
        ungroup() %>%
        select(-all_of("eclipseTcrGroupSize"))
    obj@meta.data <- obj@meta.data %>%
        mutate(cloneSize = factor(case_when(
            .data$clonalProportion <= 1 & .data$clonalProportion > 0.1
            ~ "Hyperexpanded (0.1 < X <= 1)",
            .data$clonalProportion <= 0.1 & .data$clonalProportion > 0.01
            ~ "Large (0.01 < X <= 0.1)",
            .data$clonalProportion <= 0.01 & .data$clonalProportion > 0.001
            ~ "Medium (0.001 < X <= 0.01)",
            .data$clonalProportion <= 0.001 & .data$clonalProportion > 1e-04
            ~ "Small (1e-04 < X <= 0.001)",
            .data$clonalProportion <= 1e-04 & .data$clonalProportion > 0
            ~ "Rare (0 < X <= 1e-04)",
            .data$clonalProportion <= 0
            ~ "None ( < X <= 0)",
            is.na(.data$clonalProportion)
            ~ NA
        )))

    ### Making a column that notes whether the cells in the clone appear
    ### to be MAIT or iNKT based on TCR segment usage
    obj@meta.data <- obj@meta.data %>%
        mutate(
            scR_mait_conv = case_when(
                str_detect(.data$scR_CTgene, "TRAV1-2\\.TRAJ33") ~
                    "MAIT Conventional (TRAV1-2 TRAJ33)",
                TRUE ~
                    "Normal"
            ),
            scR_mait_unconv = case_when(
                str_detect(
                    .data$scR_CTgene,
                    "TRAV1-2\\.TRAJ12|TRAV1-2\\.TRAJ20"
                )
                ~ "MAIT Unconventional (TRAV1-2 TRAJ12 or TRAJ20)",
                TRUE ~ "Normal"
            ),
            scR_inkt = case_when(
                str_detect(.data$scR_CTgene, "TRAV10\\.TRAJ18") &
                    str_detect(.data$scR_CTgene, "TRBV25-1")
                ~ "iNKT",
                TRUE ~ "Normal"
            )
        )
    obj@meta.data <- obj@meta.data %>%
        mutate(scR_unconventional_count = 3 -
            (.data$scR_mait_conv == "Normal") -
            (.data$scR_mait_unconv == "Normal") -
            (.data$scR_inkt == "Normal"))
    obj@meta.data <- obj@meta.data %>%
        mutate(cell_unconventional_subset = case_when(
            .data$scR_unconventional_count == 0 ~ "Conventional",
            .data$scR_unconventional_count == 1 &
                .data$scR_mait_conv != "Normal"
            ~ .data$scR_mait_conv,
            .data$scR_unconventional_count == 1 &
                .data$scR_mait_unconv != "Normal"
            ~ .data$scR_mait_unconv,
            .data$scR_unconventional_count == 1 & .data$scR_inkt != "Normal"
            ~ .data$scR_inkt,
            .data$scR_unconventional_count > 1
            ~ "Multiple Unconventional Types"
        ))

    df <- obj@meta.data %>%
        group_by(
            .data$CTaa,
            .data$cell_unconventional_subset,
            .data[[group]]
        ) %>%
        summarise(count = n()) %>%
        ungroup() %>%
        filter(.data$cell_unconventional_subset != "Conventional") %>%
        arrange(desc(.data$count))

    if (nrow(df) > 0) {
        obj@meta.data <- obj@meta.data %>%
            ungroup() %>%
            mutate(CTaa_group = case_when(
                !is.na(.data$CTaa) & !is.na(.data[[group]])
                ~ paste(.data$CTaa, .data[[group]], sep = " "),
                TRUE ~ NA
            ))
        df <- df %>%
            group_by(.data$CTaa, .data[[group]]) %>%
            summarise(
                clone_unconventional_subset =
                    first(.data$cell_unconventional_subset)
            )
        df <- df %>%
            ungroup() %>%
            mutate(CTaa_group = case_when(
                !is.na(.data$CTaa) & !is.na(.data[[group]])
                ~ paste(.data$CTaa, .data[[group]], sep = " "),
                TRUE ~ NA
            ))
        df <- df %>%
            ungroup() %>%
            select(all_of(c("CTaa_group", "clone_unconventional_subset")))
        
        if ("clone_unconventional_subset" %in% colnames(obj@meta.data)) {
            obj@meta.data <- obj@meta.data %>%
                select(-all_of("clone_unconventional_subset"))
        }
        
        obj@meta.data <- obj@meta.data %>% left_join(df, by = "CTaa_group")
        obj@meta.data <- obj@meta.data %>%
            mutate(clone_unconventional_subset = case_when(
                is.na(.data$clone_unconventional_subset)
                ~ "Conventional",
                !is.na(.data$clone_unconventional_subset)
                ~ .data$clone_unconventional_subset
            )) %>%
            select(-all_of("CTaa_group"))
        if (verbose == TRUE) {
            message(
                sum(obj$clone_unconventional_subset != "Conventional"),
                " cells were detected in clones that appear to ",
                "be MAIT or iNKT cells"
            )
        }
    }

    obj@meta.data <- obj@meta.data %>% select(-all_of(c(
        "scR_mait_conv",
        "scR_mait_unconv",
        "scR_inkt",
        "scR_unconventional_count",
        "tcrEclipseGroup",
        "tcrEclipse_barcode"))
    )
    
    obj@meta.data <- as.data.frame(obj@meta.data)
    rownames(obj@meta.data) <- obj$barcode

    if (verbose == TRUE) {
        message(
            round(sum(obj$cloneCallSource == "VDJdive", na.rm = TRUE) /
                nrow(obj@meta.data) * 100, 1),
            "% of cells in the Seurat object were assigned a clonotype using ",
            "high confidence prediction from VDJdive::clonoStats()"
        )
        message(
            round(sum(!is.na(obj$CTaa)) / nrow(obj@meta.data) * 100, 1),
            "% of cells in the Seurat object were assigned a clonotype ",
            "in the end (either from prediction or traditional analysis)"
        )
        message("These final clonal annotations are stored in CTaa")
        message("Summary of the top 5 largest clones:")

        big_clones <- obj@meta.data %>%
            group_by(.data[[group]]) %>%
            mutate(group_count = sum(!is.na(.data$CTaa))) %>%
            ungroup() %>%
            group_by(.data$CTaa, .data[[group]]) %>%
            summarise(
                Clone_Size = n(),
                Percent_of_Repertoire = 100 * n() /
                    first(.data$group_count)
            ) %>%
            select(all_of(c(
                "CTaa",
                group,
                "Clone_Size",
                "Percent_of_Repertoire"
                ))
            ) %>%
            arrange(desc(.data$Clone_Size)) %>%
            filter(!is.na(.data$CTaa)) %>%
            head(n = 5)
        print(big_clones)
    }
    
    ### Cleaning up object before returning to user
    if ("originalBarcodeBeforeECLIPSE" %in% colnames(obj@meta.data)) {
        obj@meta.data <- obj@meta.data %>%
            mutate(barcode = .data$originalBarcodeBeforeECLIPSE) %>%
            select(-all_of("originalBarcodeBeforeECLIPSE"))
    } else if ("originalBarcodeBeforeECLIPSE" %in% colnames(obj@meta.data)
    == FALSE) {
        obj@meta.data <- obj@meta.data %>% select(-all_of("barcode"))
    }

    obj@meta.data <- obj@meta.data %>% select(-all_of("VDJdive_clone"))

    if (detailedEclipseStats == FALSE) {
        obj@meta.data <- obj@meta.data %>% select(-all_of(c(
            "EM_max_prop",
            "EM_2ndmax_prop",
            "EM_max_2ndmax_ratio"))
        )
    }

    return(obj@meta.data)
}


shannon <- function(vec) {

    vec2 <- data.frame("clone" = vec)
    vec2 <- vec2 %>% group_by(clone) %>% summarise(p = n() / length(vec))
    
    total <- 0
    
    for (i in 1:nrow(vec2)){
    
        extra <- -1 * vec2[i, "p"] * log(vec2[i, "p"])
        total <- total + extra
        
    }
    
    return(total)
}



ECLIPSE_noEM <- function(seurat_object,
                         contig_type = "all",
                         contig_folder_paths = NULL,
                         contig_file_paths = NULL,
                         batch,
                         donor,
                         group,
                         original_barcode,
                         format = "Blank",
                         detailedEclipseStats = FALSE,
                         combineTcrArgs = list(FALSE, FALSE, FALSE, TRUE),
                         allowDoubleChains = "both",
                         returnFilteredContigs = FALSE,
                         verbose = TRUE) {
    #### THIS FIRST PART IS READING THE DATA IN AND DOING QC
    ### Checking the input data is usable
    if (!contig_type %in% c(
        "ALL", "All", "all",
        "Filtered", "FILTERED", "filtered",
        "all_contig_annotations.csv",
        "filtered_contig_annotations.csv",
        "manual"
    )) {
        stop(
            "Incorrect contig type. Should be all or filtered. ",
            "Returning nothing"
        )
    }

    if (is.null(contig_folder_paths) & contig_type != "manual") {
        stop(
            "Please provide folders to contig_folder_paths that have the ",
            "contig annotations in them."
        )
    }

    if (is.null(contig_file_paths) & contig_type == "manual") {
        stop(
            "No contig file paths provided. ",
            "Please provide a vector of file path names to contig_file_paths"
        )
    }

    if (contig_type != "manual" &
        (!is.vector(contig_folder_paths) | is.list(contig_folder_paths) |
        length(contig_folder_paths) == 0)) {
        stop("contig_folder_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (contig_type == "manual" &
        (!is.vector(contig_file_paths) | is.list(contig_file_paths) |
        length(contig_file_paths) == 0)) {
        stop("contig_file_paths must be a vector and of length >=1. ",
            "Please try again.")
    }

    if (!is.null(contig_file_paths) & !is.null(contig_folder_paths)) {
        stop("Please do not use contig_file_paths and contig_folder_paths. ",
            "Please try again.")
    }

    if (!batch %in% colnames(seurat_object@meta.data)) {
        stop(
            "Batch column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!donor %in% colnames(seurat_object@meta.data)) {
        stop(
            "Donor column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!group %in% colnames(seurat_object@meta.data)) {
        stop(
            "Group column name provided is not in the Seurat object. ",
            "Please try again."
        )
    }

    if (!original_barcode %in% colnames(seurat_object@meta.data)) {
        stop(
            "Original barcode column name provided is not ",
            "in the Seurat object. ",
            "Please try again."
        )
    }

    if (!is.logical(detailedEclipseStats)) {
        stop(
            "TRUE or FALSE not provided to detailedEclipseStats. ",
            "Please try again."
        )
    }

    if (format %in% c("None", "none", "Blank") == FALSE) {
        stop(
            "None, none, or Blank not provided to format. ",
            "Please try again."
        )
    }

    if (!is.list(combineTcrArgs) | length(combineTcrArgs) != 4 |
        !is.logical(combineTcrArgs[[1]]) | !is.logical(combineTcrArgs[[2]]) |
        !is.logical(combineTcrArgs[[3]]) | !is.logical(combineTcrArgs[[4]])) {
        stop(
            "Incorrect input to combineTcrArgs. ",
            "Please provide a list of 4 logical (TRUE or FALSE) elements."
        )
    }

    if (allowDoubleChains %in% c(
        "Alpha", "alpha",
        "beta", "Beta",
        "Both", "both",
        "Neither", "neither"
    ) == FALSE) {
        stop(
            "Alpha, Beta, Both, or Neither not provided to ",
            "allowDoubleChains. ",
            "Please try again."
        )
    }

    if (!is.logical(verbose)) {
        stop(
            "TRUE or FALSE not provided to verbose. ",
            "Please try again."
        )
    }

    if (!is.logical(returnFilteredContigs)) {
        stop(
            "TRUE or FALSE not provided to returnFilteredContigs. ",
            "Please try again."
        )
    }


    ### Reading the data in
    if (verbose == TRUE) message("\nReading TCR contig data in...\n")

    if (contig_type != "manual") {
        contigs_list <- vector("list", length(contig_folder_paths))

        for (i in seq_along(contig_folder_paths)) {
            if (contig_type %in% c(
                "ALL", "All", "all",
                "all_contig_annotations.csv"
            )) {
                files <- list.files(contig_folder_paths[i],
                    pattern = "all_contig_annotations.csv",
                    full.names = TRUE
                )
            } else if (contig_type %in%
                c(
                    "Filtered", "filtered", "FILTERED",
                    "filtered_contig_annotations.csv"
                )) {
                files <- list.files(contig_folder_paths[i],
                    pattern = "filtered_contig_annotations.csv",
                    full.names = TRUE
                )
            }


            if (length(files) == 0) {
                stop("No contig files found, returning nothing")
            } else if (length(files) > 1) {
                stop("Multiple contig files found, returning nothing")
            } else if (length(files) == 1) {
                contigs_list[[i]] <- read.csv(files[1])
                contigs_list[[i]]$barcode <- paste(i, "_",
                    contigs_list[[i]]$barcode,
                    sep = ""
                )
            }
        }
    } else if (contig_type == "manual") {
        contigs_list <- vector("list", length(contig_file_paths))

        for (i in seq_along(contig_file_paths)) {
            contigs_list[[i]] <- read.csv(contig_file_paths[i])

            if (nrow(contigs_list[[i]]) == 0) {
                stop(
                    "File path ", i, " did not lead to any ",
                    "contigs being loaded in. ",
                    "Check the path name and that it is a .csv file"
                )
            }

            contigs_list[[i]]$barcode <-
                paste(i, "_", contigs_list[[i]]$barcode, sep = "")
        }
    }

    ### Joins together all batches and adds data to them
    ### about which sample/donor/batch each contig was from
    contigs_full <- bind_rows(contigs_list)

    if ("sample" %in% colnames(contigs_full)) {
        contigs_full <- contigs_full %>% select(-all_of("sample"))
    }

    if (format %in% c("None", "none")) {
        contigs_full[contigs_full == "None"] <- ""
    }

    seurat_object$tcrEclipse_barcode <-
        paste(seurat_object@meta.data[[batch]],
            "_",
            seurat_object@meta.data[[original_barcode]],
            sep = ""
        )
    seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)

    cell_data <- seurat_object@meta.data %>%
        select(all_of(c(
            "tcrEclipse_barcode",
            unique(c(batch, donor, group)),
            "seurat_object_barcode"))
        )
    cell_data <- cell_data %>% filter(!is.na(.data[[batch]]))

    if (nrow(cell_data) != n_distinct(cell_data$tcrEclipse_barcode)) {
        stop(
            "There are likely multiple sequencing runs ",
            "labeled as one batch. ",
            "Please relabel the batch column such that each batch ",
            "is a unique sequncing run. ",
            "Returning nothing"
        )
    }

    all <- contigs_full %>%
        left_join(cell_data,
            by = c("barcode" = "tcrEclipse_barcode"),
            keep = FALSE
        )
    all <- all %>%
        mutate(barcode = .data$seurat_object_barcode) %>%
        select(-all_of("seurat_object_barcode"))
    if ("barcode" %in% colnames(seurat_object@meta.data)) {
        seurat_object@meta.data <- seurat_object@meta.data %>%
            mutate(originalBarcodeBeforeECLIPSE = .data$barcode) %>%
            select(-all_of("barcode"))
    }

    seurat_object@meta.data <- seurat_object@meta.data %>%
        dplyr::rename("barcode" = all_of("seurat_object_barcode"))

    if (verbose == TRUE) {
        message("QC Filtering:")
        message(" - ", nrow(all), " initial contigs found before filtering")
    }

    ### Renaming barcodes, filtering to only cells in seurat object,
    ### removing cells without a CDR3,
    ### and renaming QC metrics from cellranger so
    ### these contigs aren't thrown out
    all <- all %>% filter(.data$barcode %in% seurat_object$barcode)
    if (verbose == TRUE) {
        message(" - ", nrow(all), " contigs were from true cells ",
            "(coming from ", n_distinct(all$barcode), " cells)")
    }
    if (nrow(all) == 0) {
        stop(
            "Check that your barcodes are the original ",
            "(i.e. likely end with -1 and match the contig files)"
        )
    }

    if (verbose == TRUE) {
        message(" - ",
            sum(all$high_confidence %in% c("true", "True", "TRUE") == FALSE),
            " contigs were removed for not being high confidence as ",
            "defined by Cell Ranger"
        )
    }
    all <- all %>%
        filter(.data$high_confidence %in% c("true", "True", "TRUE"))
    all <- all %>% mutate(
        productive = "true",
        is_cell = "true",
        full_length = "true"
    )

    count_before_custom_filtering <- nrow(all)

    all <- all %>% filter(!.data$cdr3 %in% c("", "None") &
                            !.data$cdr3_nt %in% c("", "None"))
    all <- all %>% filter(.data$chain %in% c("TRA", "TRB", "TRG", "TRD"))

    ### Removing cells with abnormally long CDR3s
    all <- all %>% mutate(cdr3_length = str_length(.data$cdr3))

    all <- all %>% filter(.data$cdr3_length <= 30)

    ### There are TCRs with TRDV and then TRAJ and TRAC.
    ### Cellranger calls them as TRD,
    ### but they are actually TRA since they pair with a beta chain and
    ### have no J segment.
    all <- all %>%
        mutate(ad_hybrid_tcr = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC"))
            ~ TRUE,
            TRUE ~ FALSE
        ))

    all <- all %>%
        mutate(chain = case_when(
            .data$chain %in% c("TRD", "Multi") &
                str_detect(.data$v_gene, "TRDV") &
                (str_detect(.data$j_gene, "TRAJ") |
                    str_detect(.data$c_gene, "TRAC")) ~ "TRA",
            TRUE ~ .data$chain
        ))

    all <- all %>% filter(.data$chain %in% c("TRA", "TRB"))

    ### Removing cells with a * at the end of their cdr3.
    ### These indicate stop codons, so removing these contigs
    all <- all %>% filter(!str_detect(.data$cdr3, "\\*"))

    ### Removing duplicated chains. 
    ### Some datasets (particularly contig files made with
    ### older version of Cellranger)
    ### have lots of duplicated rows that have identical CDR3
    ### and gene segments for the same barcode.
    ### This messes up future analysis
    all <- all %>%
        arrange(desc(.data$umis), desc(.data$reads)) %>%
        group_by(.data$barcode, .data$chain, .data$cdr3) %>%
        filter(row_number() == 1) %>%
        ungroup()

    ### Noting cells with only 1 type of chain (alpha or beta) and
    ### that chain isn't seen in any other cell.
    ### We don't want these cells to be predicted.
    all2 <- all %>%
        group_by(.data$barcode) %>%
        mutate(chainTypes = n_distinct(.data$chain)) %>%
        ungroup()
    all2 <- all2 %>%
        group_by(.data$cdr3, .data$chain, .data[[donor]]) %>%
        mutate(maxChains = max(.data$chainTypes)) %>%
        ungroup()
    if (verbose == TRUE) {
        message(" - ",
            sum(all2$maxChains == 1),
            " orphan contigs were found with no observed pairing. ",
            "Cells containing these contigs are unlikely to have ",
            "a clonotype predicted"
        )
    }
    orphanContigs <- all2 %>%
        filter(.data$maxChains == 1) %>%
        select(all_of(c("barcode", "chain", "cdr3")))

    if (verbose == TRUE) {
        message(" - ", count_before_custom_filtering - nrow(all),
            " contigs were removed for other reasons ",
            "(ex. irregular CDR3, TCR\u03B3/\u03B4 contig, duplicated, etc...)"
        )
        message(" - ",
            "After all filtering, ", nrow(all),
            " contigs from ", n_distinct(all$barcode),
            " true cells remained"
        )
    }

    if (returnFilteredContigs == TRUE) {

        all2$orphanContig <- FALSE
        all2$orphanContig[all2$maxChains == 1] <- TRUE

        all2 <- all2 %>%
            select(-all_of(c(
                "cdr3_length",
                "ad_hybrid_tcr",
                "chainTypes",
                "maxChains"))
                )
        return(all2)
    }

    rm(contigs_full, contigs_list)

    if (verbose == TRUE) message("\nFinding clones with extra chains...\n")

    donors <- unique(as.vector(all[, donor])[[1]])
    new_contigs <- vector("list", length(donors))

    for (i in seq_along(donors)) {
        ### THIS PART IS FINDING THE COMBO CLONES AND
        ### REWRITING THE CONTIG FILES
        if (allowDoubleChains %in% c("Neither", "neither") == FALSE &
            verbose == TRUE) {
            message(
                "Running donor ", i, "/",
                length(donors), " (", donors[i], "):"
            )
        }

        slim <- all2 %>% filter(.data[[donor]] == donors[i])

        if (allowDoubleChains %in% c("Alpha", "alpha", "Both", "both")) {
            ### Calling clones with multiple alpha chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "alpha",
                printStats = verbose
            )
        }

        if (allowDoubleChains %in% c("Beta", "beta", "Both", "both")) {
            ### Calling clones with multiple beta chains
            ### and then rewriting the contig files to address this
            slim <- findThreeChainClones(
                contigs = slim,
                chain = "beta",
                printStats = verbose
            )
        }

        new_contigs[[i]] <- slim

        if (verbose == TRUE) message()
    }

    contigs_post <- bind_rows(new_contigs)

    ### Now taking these modified contigs and feeding to scRepertoire
    seurat_object$tcrEclipseGroup <- seurat_object@meta.data[[group]]
    
    contig_list <- createHTOContigList(all,
                                       seurat_object,
                                       group.by = "tcrEclipseGroup"
    )
    
    no_vdjd <- combineTCR(contig_list)
    
    seurat_object <- combineExpression(no_vdjd,
                                       seurat_object,
                                       cloneCall = "aa",
                                       proportion = TRUE
    )
    
    seurat_object$scR_CTaa <- seurat_object$CTaa
    seurat_object@meta.data <- seurat_object@meta.data %>%
      select(-all_of(c(
        "CTgene",
        "CTnt",
        "CTaa",
        "CTstrict",
        "clonalProportion",
        "clonalFrequency",
        "cloneSize"))
      )
    
    contig_list <- createHTOContigList(contigs_post,
        seurat_object,
        group.by = "tcrEclipseGroup"
    )
    
    no_vdjd <- combineTCR(contig_list)
    
    obj <- combineExpression(no_vdjd,
        seurat_object,
        cloneCall = "aa",
        proportion = TRUE
    )

    rm(seurat_object)

    obj@meta.data <- as.data.frame(obj@meta.data)
    rownames(obj@meta.data) <- obj$barcode

    return(obj)
}






compareTRUST4 <- function(output_folder,
                          seurat_object,
                          library,
                          samples = NULL,
                          max_dist) {

  
    trust4_rcc_dirs <- list.files(output_folder,
                                   full.names = TRUE,
                                   recursive = TRUE,
                                   pattern = "airr\\.tsv") 
    
    trust4_rcc_dirs <- trust4_rcc_dirs[str_detect(trust4_rcc_dirs, "barcode")]
    
    trust4_rcc_files <- lapply(trust4_rcc_dirs, read_tsv)
    
    sample_map <- data.frame("dirs" = trust4_rcc_dirs,
                             "Sample" = str_replace(trust4_rcc_dirs, "^.*rcc/|^.*mel/", ""))
    
    sample_map$Sample <- str_replace(sample_map$Sample, "/.*$", "")
    
    for (i in 1:length(trust4_rcc_files)) {
      
      trust4_rcc_files[[i]]$Sample <- sample_map$Sample[i]
      
    }
    
    contigs <- as.data.frame(do.call(rbind, trust4_rcc_files))
    
    print("Contigs per sample:")
    print(table(contigs$Sample))
    
    contigs <- contigs %>% filter(locus %in% c("TRA", "TRB") & !is.na(junction_aa)) %>%
                           mutate(Sample = str_replace(Sample, "-", "_"))
    contigs$Sample[contigs$Sample == "S10_T"] <- "S10_T2"
    
    if (!is.null(samples)) {
      
      contigs <- contigs %>% ungroup() %>%
                             dplyr::rename("directory" = Sample)
      
      contigs <- contigs %>% left_join(samples, by = "directory")
      
    }
    
    contigs2 <- contigs %>% select(cell_id, Sample, junction_aa, locus) %>%
                            unite(col = "key", c(cell_id, Sample), sep = " ", remove = FALSE)
    
    called_cells <- seurat_object@meta.data %>% select(Cell_barcode, Sample, donor, cloneCallSource, scR_CTaa, CTaa) %>%
                                                filter(!is.na(CTaa) & cloneCallSource == "VDJdive") %>%
                                                mutate(trimmed_bc = str_replace(Cell_barcode, "-.*$", "")) %>%
                                                unite("key", c(trimmed_bc, Sample), sep = " ", remove = FALSE) %>%
                                                select(key, Cell_barcode, scR_CTaa, CTaa, donor)
    
    full <- contigs2 %>% left_join(called_cells, by = "key") %>%
                         filter(!is.na(CTaa))
    
    full2 <- full %>% select(junction_aa, locus, Cell_barcode, scR_CTaa, CTaa) %>%
                      pivot_longer(cols = c(scR_CTaa, CTaa), names_to = "type", values_to = "cdr3")
    
    full2 <- full2 %>% mutate(cdr3 = str_replace_all(cdr3, ";", ":::")) %>%
                       separate(cdr3, into = c("TRA", "TRB"), sep = "_", remove = FALSE)
    
    full2 <- full2 %>% pivot_longer(cols = c(TRA, TRB), names_to = "chain", values_to = "chain_cdr3")
    full2$chain_cdr3[full2$chain_cdr3 == "NA"] <- paste(rep("MISSING", 5), collapse = " ")
    
    full2 <- full2 %>% separate(chain_cdr3, into = c("chain_cdr3_1", "chain_cdr3_2"), remove = FALSE, sep = ":::")
    
    full2 <- full2 %>% pivot_longer(c(chain_cdr3_1, chain_cdr3_2), names_to = "chain_pos", values_to = "final_cdr3")
    full2 <- full2 %>% filter(!is.na(final_cdr3) & locus == chain)
    
    full2 <- full2 %>% mutate(lv_dist = stringdist(final_cdr3, junction_aa, method = "lv", weight = c(d = 1, i = 0.04, s = 0.2, t = 1)))
    
    full3 <- full2 %>% group_by(junction_aa, chain, Cell_barcode, type, cdr3) %>%
                       summarise(count = n(),
                                 mx_dist = max(lv_dist),
                                 mn_dist = min(lv_dist),
                                 max_chain = final_cdr3[which.max(lv_dist)],
                                 min_chain = final_cdr3[which.min(lv_dist)]) %>%
                       ungroup()
    
    
    full4 <- full3 %>% group_by(junction_aa, chain, Cell_barcode) %>% 
                       summarise(Count = sum(count),
                                 min_ctaa_dist = mn_dist[type == "CTaa"],
                                 max_ctaa_dist = mx_dist[type == "CTaa"],
                                 min_scr_dist = mn_dist[type == "scR_CTaa"],
                                 max_scr_dist = mx_dist[type == "CTaa"],
                                 min_ctaa_chain = min_chain[type == "CTaa"],
                                 min_scr_chain = min_chain[type == "scR_CTaa"],
                                 CTaa = cdr3[type == "CTaa"],
                                 scR_CTaa = cdr3[type == "scR_CTaa"]) %>%
                        ungroup()
     
    full4 <- full4 %>% mutate(seq = case_when(
                                              min_scr_dist > max_dist ~ "not seq",
                                              min_scr_dist <= max_dist ~ "seq"),
                              pred = case_when(
                                               min_ctaa_dist > max_dist ~ "not pred",
                                               min_ctaa_dist <= max_dist ~ "pred"),
                              seq_pred = case_when(
                                                   min_scr_dist > max_dist & min_ctaa_dist > max_dist ~ "not seq\nnot pred",
                                                   min_scr_dist <= max_dist & min_ctaa_dist > max_dist ~ "seq\nnot pred",
                                                   min_scr_dist > max_dist & min_ctaa_dist <= max_dist ~ "not seq\npred",
                                                   min_scr_dist <= max_dist & min_ctaa_dist <= max_dist ~ "seq\npred"))
    
    
    full <- full %>% unite("key", c(junction_aa, Cell_barcode), sep = " ", remove = FALSE)
    full4 <- full4 %>% unite("key", c(junction_aa, Cell_barcode), sep = " ", remove = FALSE) %>% ungroup()
    
    final <- full4 %>% select(key, min_ctaa_dist, min_scr_dist, seq, pred, seq_pred, min_ctaa_chain, min_scr_chain) %>%
                       left_join(full, by = "key", keep = FALSE)
    
    print(table(final$Sample, final$seq_pred))
    
    p1 <- final %>% ggplot() + geom_histogram(aes(x = min_ctaa_dist), bins = 100) + geom_vline(xintercept = max_dist)
    print(p1)
    
    seq <- final %>% filter(seq == "not seq") %>%
                     group_by(donor) %>%
                     mutate(total_size = n()) %>%
                     group_by(pred, donor) %>%
                     summarise(pct = 100 * n() / dplyr::first(total_size),
                               subset_size = n(),
                               totalSize = dplyr::first(total_size)) %>%
                     ungroup()
    
    seq <- seq %>% complete(pred, donor)
    
    seq <- seq %>% arrange(desc(pct)) %>%
                   group_by(donor) %>%
                   mutate(totalSize = case_when(
                                                is.na(totalSize) ~ dplyr::first(totalSize),
                                                !is.na(totalSize) ~ totalSize)) %>%
                   ungroup()
    
    seq$pct[is.na(seq$pct)] <- 0
    seq$subset_size[is.na(seq$subset_size)] <- 0
    
    seq$totalSize[is.na(seq$totalSize)] <- 0
    
    seq$method <- "TRUST4"
    seq$library <- library
    final$method <- "TRUST4"
    final$library <- library
    
    return(list(seq, final))
}



compareMiXCR <- function(output_folder,
                         seurat_object,
                         library,
                         samples = NULL,
                         max_dist) {
  
    mixcr_dirs <- list.files(output_folder,
                                   full.names = TRUE,
                                   recursive = TRUE,
                                   pattern = "clones\\.tsv") 
    

    mixcr_files <- lapply(mixcr_dirs, read_tsv)
    
    sample_map <- data.frame("dirs" = mixcr_dirs,
                             "Sample" = str_replace(mixcr_dirs, "^.*rcc/|^.*mel/", ""))
    
    sample_map$Sample <- str_replace(sample_map$Sample, "/.*$", "")
    
    for (i in 1:length(mixcr_files)) {
      
      mixcr_files[[i]]$Sample <- sample_map$Sample[i]
      
    }
    
    contigs <- as.data.frame(do.call(rbind, mixcr_files))
    
    print("Contigs per sample:")
    print(table(contigs$Sample))
    
    contigs <- contigs %>% filter(topChains %in% c("TRA", "TRB") & !is.na(aaSeqCDR3)) %>%
                           mutate(Sample = str_replace(Sample, "-", "_"))
    contigs$Sample[contigs$Sample == "S10_T"] <- "S10_T2"
    
    if (!is.null(samples)) {
      
      contigs <- contigs %>% ungroup() %>%
                             dplyr::rename("directory" = Sample)
      
      contigs <- contigs %>% left_join(samples, by = "directory")
      
    }
    
    contigs2 <- contigs %>% select(tagValueCELL, Sample, aaSeqCDR3, nSeqCDR3, topChains) %>%
                            unite(col = "key", c(tagValueCELL, Sample), sep = " ", remove = FALSE)
    
    called_cells <- seurat_object@meta.data %>% select(Cell_barcode, Sample, donor, cloneCallSource, scR_CTaa, CTaa) %>%
                                                filter(!is.na(CTaa) & cloneCallSource == "VDJdive") %>%
                                                mutate(trimmed_bc = str_replace(Cell_barcode, "-.*$", "")) %>%
                                                unite("key", c(trimmed_bc, Sample), sep = " ", remove = FALSE) %>%
                                                select(key, Cell_barcode, scR_CTaa, CTaa, donor)
    
    full <- contigs2 %>% left_join(called_cells, by = "key") %>%
                         filter(!is.na(CTaa))
    
    full <- full %>% dplyr::rename("cell_id" = tagValueCELL,
                                   "junction_aa" = aaSeqCDR3,
                                   "junction_nt" = nSeqCDR3,
                                   "locus" = topChains)

    full2 <- full %>% select(junction_aa, junction_nt, locus, Cell_barcode, scR_CTaa, CTaa) %>%
                  pivot_longer(cols = c(scR_CTaa, CTaa), names_to = "type", values_to = "cdr3")
    
    full2 <- full2 %>% mutate(cdr3 = str_replace_all(cdr3, ";", ":::")) %>%
                       separate(cdr3, into = c("TRA", "TRB"), sep = "_", remove = FALSE)
    
    full2 <- full2 %>% pivot_longer(cols = c(TRA, TRB), names_to = "chain", values_to = "chain_cdr3")
    full2$chain_cdr3[full2$chain_cdr3 == "NA"] <- paste(rep("MISSING", 5), collapse = " ")
    
    full2 <- full2 %>% separate(chain_cdr3, into = c("chain_cdr3_1", "chain_cdr3_2"), remove = FALSE, sep = ":::")
    
    full2 <- full2 %>% pivot_longer(c(chain_cdr3_1, chain_cdr3_2), names_to = "chain_pos", values_to = "final_cdr3")
    full2 <- full2 %>% filter(!is.na(final_cdr3) & locus == chain)
    
    full2 <- full2 %>% mutate(junction_nt_1 = str_replace(junction_nt, "^.", ""),
                              junction_nt_2 = str_replace(junction_nt, "^..", ""))
    
    full2$junction_aa_0 <- suppressWarnings(as.character(translate(DNAStringSet(full2$junction_nt), no.init.codon = TRUE)))
    full2$junction_aa_1 <- suppressWarnings(as.character(translate(DNAStringSet(full2$junction_nt_1), no.init.codon = TRUE)))
    full2$junction_aa_2 <- suppressWarnings(as.character(translate(DNAStringSet(full2$junction_nt_2), no.init.codon = TRUE)))
    
    full2 <- full2 %>% mutate(lv_dist_original = stringdist(final_cdr3, junction_aa, method = "lv", weight = c(d = 1, i = 0.04, s = 0.2, t = 1)),
                              lv_dist_0 = stringdist(final_cdr3, junction_aa_0, method = "lv", weight = c(d = 1, i = 0.04, s = 0.2, t = 1)),
                              lv_dist_1 = stringdist(final_cdr3, junction_aa_1, method = "lv", weight = c(d = 1, i = 0.04, s = 0.2, t = 1)),
                              lv_dist_2 = stringdist(final_cdr3, junction_aa_2, method = "lv", weight = c(d = 1, i = 0.04, s = 0.2, t = 1)))
    
    full2 <- full2 %>% rowwise() %>% 
                       mutate(lv_dist = min(c(lv_dist_original, lv_dist_0, lv_dist_1, lv_dist_2)),
                              junction_aa_right_frame = c(junction_aa, junction_aa_0, junction_aa_1, junction_aa_2)[which.min(c(lv_dist_original, lv_dist_0, lv_dist_1, lv_dist_2))])
                                      
    full3 <- full2 %>% group_by(junction_aa, chain, Cell_barcode, type, cdr3) %>%
                       summarise(count = n(),
                                 mx_dist = max(lv_dist),
                                 mn_dist = min(lv_dist),
                                 max_chain = final_cdr3[which.max(lv_dist)],
                                 min_chain = final_cdr3[which.min(lv_dist)],
                                 junction_aa_right_frame = junction_aa_right_frame[which.min(lv_dist)]) %>%
                       ungroup()
    
    
    full4 <- full3 %>% group_by(junction_aa, chain, Cell_barcode) %>% 
                       summarise(Count = sum(count),
                                 min_ctaa_dist = mn_dist[type == "CTaa"],
                                 max_ctaa_dist = mx_dist[type == "CTaa"],
                                 min_scr_dist = mn_dist[type == "scR_CTaa"],
                                 max_scr_dist = mx_dist[type == "CTaa"],
                                 min_ctaa_chain = min_chain[type == "CTaa"],
                                 min_scr_chain = min_chain[type == "scR_CTaa"],
                                 CTaa = cdr3[type == "CTaa"],
                                 scR_CTaa = cdr3[type == "scR_CTaa"],
                                 junction_aa_right_frame_ctaa = junction_aa_right_frame[type == "CTaa"],
                                 junction_aa_right_frame_scr = junction_aa_right_frame[type == "scR_CTaa"]) %>%
                        ungroup()
     
    full4 <- full4 %>% mutate(seq = case_when(
                                              min_scr_dist > max_dist ~ "not seq",
                                              min_scr_dist <= max_dist ~ "seq"),
                              pred = case_when(
                                               min_ctaa_dist > max_dist ~ "not pred",
                                               min_ctaa_dist <= max_dist ~ "pred"),
                              seq_pred = case_when(
                                                   min_scr_dist > max_dist & min_ctaa_dist > max_dist ~ "not seq\nnot pred",
                                                   min_scr_dist <= max_dist & min_ctaa_dist > max_dist ~ "seq\nnot pred",
                                                   min_scr_dist > max_dist & min_ctaa_dist <= max_dist ~ "not seq\npred",
                                                   min_scr_dist <= max_dist & min_ctaa_dist <= max_dist ~ "seq\npred"))
    
    
    full <- full %>% unite("key", c(junction_aa, Cell_barcode), sep = " ", remove = FALSE)
    full4 <- full4 %>% unite("key", c(junction_aa, Cell_barcode), sep = " ", remove = FALSE) %>% ungroup()
    
    final <- full4 %>% select(key, min_ctaa_dist, min_scr_dist, seq, pred, seq_pred, junction_aa_right_frame_ctaa, junction_aa_right_frame_scr, min_ctaa_chain, min_scr_chain) %>%
                       left_join(full, by = "key", keep = FALSE)
    
    print(table(final$Sample, final$seq_pred))
    
    p1 <- final %>% ggplot() + geom_histogram(aes(x = min_ctaa_dist), bins = 100) + geom_vline(xintercept = max_dist)
    print(p1)
    
    seq <- final %>% filter(seq == "not seq") %>%
                     group_by(donor) %>%
                     mutate(total_size = n()) %>%
                     group_by(pred, donor) %>%
                     summarise(pct = 100 * n() / dplyr::first(total_size),
                               subset_size = n(),
                               totalSize = dplyr::first(total_size)) %>%
                     ungroup()
    
    seq <- seq %>% complete(pred, donor)
    
    seq <- seq %>% arrange(desc(pct)) %>%
                   group_by(donor) %>%
                   mutate(totalSize = case_when(
                                                is.na(totalSize) ~ dplyr::first(totalSize),
                                                !is.na(totalSize) ~ totalSize)) %>%
                   ungroup()
    
    seq$pct[is.na(seq$pct)] <- 0
    seq$subset_size[is.na(seq$subset_size)] <- 0
    
    seq$totalSize[is.na(seq$totalSize)] <- 0
    
    seq$method <- "MiXCR"
    seq$library <- library
    final$method <- "MiXCR"
    final$library <- library
    
    return(list(seq, final))
  
}

makeBoxplot <- function(data, x, y, group, x_title, y_title, pval, colors, ylim, x_text, title = NULL) {

  plot <- data %>% ggplot(aes(.data[[x]], .data[[y]])) + 
                   geom_boxplot(aes(fill = .data[[x]], color = .data[[x]]), outlier.size = -1, lwd = 0.3, width = 0.5, alpha = 0.2, fatten = 5) +
                   geom_path(aes(group = .data[[group]]), alpha = 0.1, lwd = 0.2) +
                   geom_point(aes(fill = .data[[x]]), size = 1.5, color = "black", shape = 21, stroke = 0.3) + 
                   theme_bw() +
                   labs(x = x_title, y = y_title) +
                   stat_pvalue_manual(pval, tip.length = 0, side = 4, label = "p.adj2", size = 1.7, bracket.size = 0.2) +
                   theme(axis.text = element_text(size = 6, color = "black"),
                         axis.title = element_text(size = 8, color = "black"),
                         legend.position = "none",
                         panel.grid = element_blank(),
                         panel.border = element_rect(linewidth = 0.7, color = "black"),
                         axis.ticks = element_line(linewidth = 0.3, color = "black"),
                         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, margin = margin(t = 0.5)),
                         axis.title.x = element_text(hjust = 0.5, vjust = 0, size = 8, color = "black"),
                         axis.text.y = element_text(margin = margin(r = 0.5)),
                         axis.ticks.length = unit(0.05, "cm")) +
                   scale_fill_manual(values = colors) +
                   scale_color_manual(values = colors) +
                   scale_x_discrete(labels = x_text) +
                   scale_y_continuous(guide = guide_axis(minor.ticks = TRUE), limits = c(ylim[1], ylim[2]))
  
  if (!is.null(title)) {
    plot <- plot + ggtitle(title) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5))
  }
  
  print(plot)
}


makeLollipop <- function(data, lab_x, Color, clones, colors, type) {
  

    plot <- data %>% ggplot(aes(x = new_x, y = pct, color = .data[[Color]])) +
                     geom_hline(yintercept = 104, linewidth = 0.3) + 
                     geom_hline(yintercept = 0, linewidth = 0.3) + 
                     geom_vline(xintercept = 0.51, linewidth = 0.3) + 
                     geom_vline(xintercept = 5.49, linewidth = 0.3) + 
                     geom_segment(aes(xend = new_x, yend = 0, group = .data[[Color]]), color = "grey70", linetype = "dashed", linewidth = 0.3) +
                     geom_point(aes(color = .data[[Color]], group = .data[[Color]]), size = 1) +
                     theme_classic() + 
                     theme(legend.position = "none", 
                           axis.title = element_text(size = 8, color = "black"), 
                           panel.ontop = FALSE,
                           axis.line = element_blank(), 
                           axis.text.x = element_text(vjust = 1, size = 4, color = "black"), 
                           axis.ticks.length  = unit(0.05, "cm"),
                           axis.ticks.y = element_line(linewidth = 0.3),
                           axis.ticks.x = element_blank(),
                           axis.text.y = element_text(size = 6, color = "black", margin = margin(r = 1))) + 
                     scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 104), breaks = seq(0, 100, 20), guide = guide_axis(minor.ticks = TRUE)) + 
                     labs(y = "% of cells", x = "Cell fate") + 
                     coord_cartesian(clip = "off") + 
                     scale_color_manual(breaks = clones, 
                                        values = colors,
                                        labels = clones) +
                     scale_x_discrete(labels = c("Proliferating", "Cytotoxic", "Resident\nmemory", "Exhausted", "Terminally\nexhausted"), 
                                      limits = c("1", "2", "3", "4", "5"),
                                      expand = c(0.1, 0.1))
    if (type == "missing") {
        
        a <- str_replace(clones[1], "_.*$", "")
        b <- str_replace(clones[1], "^.*_", "")
        
        plot <- plot + annotate("text", x = lab_x, y = 90, label = "CDR3 key:", hjust = 0, size = 2) + 
                       annotate("text", x = lab_x, y = 83, label = paste0("alpha[1] ~ '=' ~ '", a, "'"), parse = TRUE, hjust = 0, size = 2) + 
                       annotate("text", x = lab_x, y = 76, label = paste0("beta[1] ~ '=' ~ '", b, "'"), parse = TRUE, hjust = 0, size = 2)
      
    } else if (type == "extra") {
      
        b <- str_replace(clones[1], "^.*_", "")
        alphas <- str_replace(clones[1], "_.*$", "")
        a1 <- str_replace(alphas, ";.*$", "")
        a2 <- str_replace(alphas, "^.*;", "")
        
        plot <- plot + annotate("text", x = lab_x, y = 90, label = "CDR3 key:", hjust = 0, size = 2) + 
                       annotate("text", x = lab_x, y = 83, label = paste0("alpha[1] ~ '=' ~ '", a1, "'"), parse = TRUE, hjust = 0, size = 2) + 
                       annotate("text", x = lab_x, y = 76, label = paste0("alpha[2] ~ '=' ~ '", a2, "'"), parse = TRUE, hjust = 0, size = 2) + 
                       annotate("text", x = lab_x, y = 69, label = paste0("beta[1] ~ '=' ~ '", b, "'"), parse = TRUE, hjust = 0, size = 2)
        
    }
    
    print(plot)

  
}


showClones <- function(obj, variable, clones, colors, type, count = NULL) {

  
    plot <- FetchData(obj, c("umap_1", "umap_2", variable)) %>% ggplot() + 
                                                                         geom_point(data = . %>% filter(.data[[variable]] == "Other"), aes(x = umap_1, y = umap_2), color = "grey80", size = 0.5, shape = 16, alpha = 0.2) +
                                                                         geom_point(data = . %>% filter(.data[[variable]] != "Other"), aes(x = umap_1, y = umap_2, color = .data[[variable]]), size = 0.5, shape = 16, alpha = 1) + 
                                                                         theme_void() + 
                                                                         labs(color = "TCR\nchain\nsets") +
                                                                         theme(legend.position = "bottom", 
                                                                               legend.text = element_text(color = "black", size = 6, margin = margin(l = 1, r = 1, t = -5, b = -5), hjust = 0.5), 
                                                                               legend.title = element_text(size = 6, margin = margin(r = 1), hjust = 0.5), 
                                                                               plot.title = element_text(size = 8, color = "black", hjust = 0.5), 
                                                                               legend.key.spacing.y = unit(0.01, "cm"),
                                                                               legend.box.margin = unit(-0.2, "cm"),
                                                                               legend.margin = margin(t = -1, b = 6),
                                                                               legend.title.position = "left") +
                                                                         guides(color = guide_legend(nrow = 3, override.aes = list(size = 3), keyheight = 0.1, keywidth = 0.01)) 
  
    if (type == "missing") {
      
        if (count == 3) {
            plot <- plot + scale_color_manual(breaks = clones, 
                                              labels = list(bquote(alpha[1] ~ "-" ~ beta[1]),
                                                            bquote("?" ~ "-" ~ beta[1]),
                                                            bquote(alpha[1] ~ "-" ~ "?")), 
                                              values = colors)
        } else if (count == 2){
          
            plot <- plot + scale_color_manual(breaks = clones, 
                                              labels = list(bquote(alpha[1] ~ "-" ~ beta[1]),
                                                            bquote("?" ~ "-" ~ beta[1])),  
                                              values = colors)
        }
    } else if (type == "extra") {
        
      plot <- plot + scale_color_manual(breaks = clones, 
                                        labels = list(bquote(alpha[1] * alpha[2] ~ "-" ~ beta[1]),
                                                      bquote(alpha[2] ~ "-" ~ beta[1]),
                                                      bquote(alpha[1] ~ "-" ~ beta[1])),
                                        values = colors)
    }
  
    print(plot)  
  
}

