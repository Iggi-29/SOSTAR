#############################################
### SOSTAR annotation filtration and work ###
#############################################

#### libraries
library(tidyverse, verbose = FALSE, warn.conflicts = FALSE)
library(dplyr, verbose = FALSE, warn.conflicts = FALSE)
library(tidyr, verbose = FALSE, warn.conflicts = FALSE)
library(openxlsx, verbose = FALSE, warn.conflicts = FALSE)
library(ggplot2, verbose = FALSE, warn.conflicts = FALSE)
library(mclust, verbose = FALSE, warn.conflicts = FALSE)
library(optparse, verbose = FALSE, warn.conflicts = FALSE)

#### define inputs
option_list <- list(
  make_option(opt_str = "--SOSTAR.final.data", 
              help = "Folder with the final data of sostar", dest = "sostar_data"),
  make_option(opt_str = "--Refference.trans", 
              help = "Refference transcripts parsed to Sostar", dest = "trans")
)
opt <- parse_args(OptionParser(option_list = option_list))

#### Warnings
if (is.null(opt$sostar_data) | is.null(opt$trans)) {
  stop("Check for help")
}

#### Final data folder
final_data_folder <- opt$sostar_data
refference_ids <- opt$trans

### Create the worked_annotation and plots folder
cmd_create1 <- paste0(x = c("mkdir -p ",final_data_folder,"SOSTAR_Events/"), collapse = "")
cmd_create2 <- paste0(x = c("mkdir -p ",final_data_folder,"plots/"), collapse = "")
cmd_create3 <- paste0(x = c("mkdir -p ",final_data_folder,"plots/densities"), collapse = "")
cmd_create4 <- paste0(x = c("mkdir -p ",final_data_folder,"plots/correlations"), collapse = "")

system(cmd_create1)
system(cmd_create2)
system(cmd_create3)
system(cmd_create4)

### Define paths for plots
density_plot_place <- paste0(x =c(final_data_folder,"plots/densities/"), collapse = "")
correlation_plot_place <- paste0(x =c(final_data_folder,"plots/correlations/TPM_FPKM.png"), collapse = "")
correlation_plot_place2 <- paste0(x =c(final_data_folder,"plots/correlations/TPM_event_expression.png"), collapse = "")

### Define the place with the expression data
place_of_expression <- paste0(x =c(final_data_folder,"expression/"), collapse = "")

#### Final data names
worked_annotation <- paste0(x =c(final_data_folder,"/SOSTAR_Events/Event_annotation.xlsx"), collapse = "")
event_data <- paste0(x =c(final_data_folder,"/SOSTAR_Events/Event_data.xlsx"), collapse = "")
event_data_raw <- paste0(x =c(final_data_folder,"/SOSTAR_Events/Event_data_raw.xlsx"), collapse = "")
event_data_filt <- paste0(x =c(final_data_folder,"/SOSTAR_Events/Event_data_filt.xlsx"), collapse = "")
event_data_check <- paste0(x =c(final_data_folder,"/SOSTAR_Events/Event_data_check.xlsx"), collapse = "")

#### SOSTAR data
sostar_result_file <- paste0(x = c(final_data_folder,"/SOSTAR_annotation_table_results.xlsx"), collapse = "")
sostar_result <- openxlsx::read.xlsx(xlsxFile = sostar_result_file, sheet = 1)

sostar_result <- sostar_result %>% 
  dplyr::filter(occurence > 0)  

#### refference ids
refference_ids <- readLines(con = refference_ids)

#### Reannot the data of SOSTAR to make it workable
sostar_result_mod <- sostar_result %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    ## Change the triangle
    new_annotation = gsub(pattern = "Δ", replacement = "ES", x = annot_find),
    new_annotation = gsub(pattern = "▼", replacement = "IR", x = new_annotation),
    ## Change "-" outside () by _then_
    new_annotation = gsub(pattern = "-(?=[^()]*?(?:\\(|$))", replacement = "_then_", x = new_annotation, perl = TRUE),
    ## Change "," outside [] by _then_
    # new_annotation = gsub(pattern = ",(?=(?:[^\\[]*\\[[^\\]]*\\])*[^\\[]*$)", replacement = "_then_", x = new_annotation, perl = TRUE)
    new_annotation = {
      
      tmp <- new_annotation
      
      # Extract bracket blocks
      matches <- gregexpr("\\[[^\\]]*\\]", tmp, perl = TRUE)
      brackets <- regmatches(tmp, matches)
      
      # Replace bracket blocks with placeholders
      for (i in seq_along(brackets)) {
        if (length(brackets[[i]]) > 0) {
          placeholders <- paste0("BRACKETPLACEHOLDER", seq_along(brackets[[i]]))
          regmatches(tmp[i], gregexpr("\\[[^\\]]*\\]", tmp[i], perl = TRUE)) <- placeholders
          
          # Replace commas safely
          tmp[i] <- gsub(",", "_then_", tmp[i])
          
          # Restore bracket content
          for (j in seq_along(placeholders)) {
            tmp[i] <- gsub(placeholders[j], brackets[[i]][j], tmp[i], fixed = TRUE)
          }
        } else {
          tmp[i] <- gsub(",", "_then_", tmp[i])
        }
      }
      
      tmp
    }
    ) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(new_annotation, .after = annot_find)

#### The original annotation correction for those "normal transcripts"
sostar_result_mod <- sostar_result_mod %>% 
  dplyr::mutate(new_annotation = ifelse(annot_ref == annot_find,annot_find,new_annotation))

#### Separate the annotations by events
sostar_result_mod <- sostar_result_mod %>%
  dplyr::mutate(new_annotation2 = new_annotation) %>% 
  dplyr::relocate(new_annotation2, .after = new_annotation) %>% 
  tidyr::separate_rows(new_annotation2, sep = "_then_") %>% 
  dplyr::rename(iso_event = new_annotation2)

#### Add gene information to the events
sostar_result_mod <- sostar_result_mod %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(iso_event = paste0(x = c(gene,";",iso_event), collapse = "")) %>% 
  dplyr::ungroup()

#### Filter out those events that are characterisic of known isoforms
sostar_result_mod_list <- list() 
for (i in 1:length(unique(sostar_result_mod$gene))) {
  ### data now
  gene_now <- unique(sostar_result_mod$gene)[i]
  sostar_result_mod_now <- sostar_result_mod %>% 
    dplyr::filter(gene == gene_now)
  
  canonical_events <- sostar_result_mod_now %>% 
    dplyr::filter(grepl(x = transcript_id, pattern = "ENST")) %>% 
    dplyr::filter(!transcript_id %in% refference_ids) %>% 
    dplyr::pull(var = iso_event)
  
  sostar_result_mod_now <- sostar_result_mod_now %>% 
    dplyr::filter(!iso_event %in% canonical_events)
  
  sostar_result_mod_list[[gene_now]] <- sostar_result_mod_now
  remove(gene_now);remove(sostar_result_mod_now);remove(canonical_events)
  }

sostar_result_mod <- do.call(rbind, sostar_result_mod_list)
remove(sostar_result_mod_list)

#### Save the annotation of each of the original events
sostar_result_mod_annotated <- sostar_result_mod
sostar_result_mod_annotated_list <- list()
for (i in 1:length(unique(sostar_result_mod_annotated$gene))) {
  gene_now <- unique(sostar_result_mod_annotated$gene)[i]
  data_now <- sostar_result_mod_annotated %>%
    dplyr::filter(gene == gene_now)

  data_now <- data_now %>%
    dplyr::rowwise() %>%
    dplyr::mutate(iso_event = unlist(strsplit(x = iso_event, split = "\\;")[[1]][2])) %>%
    dplyr::ungroup()

  sostar_result_mod_annotated_list[[gene_now]] <- data_now

  remove(gene_now);remove(data_now)
  }

openxlsx::write.xlsx(x = sostar_result_mod_annotated_list, file = worked_annotation)
remove(sostar_result_mod_annotated_list);remove(sostar_result_mod_annotated)

#### Expression profile work
sostar_result_mod_expression <- sostar_result_mod
sostar_result_mod_expression <- sostar_result_mod_expression %>%
  tidyr::pivot_longer(cols = 
                        grep(x = colnames(sostar_result_mod_expression), pattern = "barcode"), 
                      names_to = "sample_name", values_to = "event_expression") %>% 
  dplyr::select(-c(occurence,event_expression)) ### OJO en el futur eliminar també event_expression pero de moment ho mantenim

sostar_result_expression <- data.frame()
for (i in 1:length(unique(sostar_result_mod_expression$sample_name))) {
  ## data now 
  sample_now <- unique(sostar_result_mod_expression$sample_name)[i]
  data_now <- sostar_result_mod_expression %>% 
    dplyr::filter(sample_name == sample_now)
  expression_now <- paste0(x = c(place_of_expression,sample_now,".expression.gtf"), collapse = "")
  expression_now <- readLines(con = expression_now)
  expression_now <- expression_now[-c(1:2)]
  expression_now <- expression_now[grep(x = expression_now, pattern = "FPKM")]
  
  expression_now <- as.data.frame(
    do.call(rbind, strsplit(expression_now, "\t")),
    stringsAsFactors = FALSE
  )
  
  colnames(expression_now) <- c(
    "chr", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes")
  expression_now <- expression_now %>% 
    dplyr::select(chr, start, end, strand, attributes) 
  expression_now <- expression_now %>%  
    dplyr::mutate(row_id = row_number()) %>% 
    dplyr::mutate(attributes = str_remove(string = attributes, pattern = "\\;$")) %>% 
    tidyr::separate_rows(attributes, sep = "; ") %>%
    tidyr::separate(attributes, into = c("key", "value"), sep = " \"") %>% 
    dplyr::mutate(value = str_remove(string = value, pattern = "\"")) %>% 
    dplyr::select(c(key, value, row_id)) %>% 
    tidyr::pivot_wider(names_from = key, values_from = value, id_cols = row_id) %>% 
    dplyr::select(-c(row_id, gene_id, ref_gene_name))
  
  data_now <- data_now %>% 
    dplyr::left_join(y = expression_now, by = "transcript_id") %>% 
    dplyr::mutate(
      cov = ifelse(is.na(cov),0,cov),
      FPKM = ifelse(is.na(FPKM),0,FPKM),
      TPM = ifelse(is.na(TPM),0,TPM)) %>% 
    # dplyr::filter(cov >= 10) %>% 
    dplyr::rename(event_expression = cov)  %>% 
    dplyr::mutate(event_expression = round(x = as.numeric(event_expression), digits = 2),
                  FPKM = round(x = as.numeric(FPKM), digits = 2),
                  TPM = round(x = as.numeric(TPM), digits = 2))
  
  sostar_result_expression <- rbind(sostar_result_expression, data_now)
  remove(data_now);remove(sample_now);remove(expression_now)
  
  } 

sostar_result_expression <- sostar_result_expression %>%
  dplyr::filter(event_expression > 0)

remove(sostar_result_mod);remove(sostar_result_mod_expression)

#### TPM FPKM correlation
corrplot <- ggplot(data = sostar_result_expression) +
  geom_point(mapping = aes(x = TPM, y = FPKM)) +
  theme_classic() +
  geom_smooth(aes(x = TPM, y = FPKM),
              method = "lm",
              se = FALSE) +
  labs(title = "TPM / FPKM correlation",
       y = "FPKM", y = "TPM")


corrplot2 <- ggplot(data = sostar_result_expression) +
  geom_point(mapping = aes(x = TPM, y = event_expression)) +
  theme_classic() +
  geom_smooth(aes(x = TPM, y = event_expression),
              method = "lm",
              se = FALSE) +
  labs(title = "TPM / event_expression correlation",
       y = "FPKM", y = "event_expression")


ggsave(plot = corrplot, filename = correlation_plot_place, width = 10, height = 10)
remove(corrplot)
ggsave(plot = corrplot2, filename = correlation_plot_place2, width = 10, height = 10)
remove(corrplot2)

 
#### Event work
n_of_samples <- length(unique(sostar_result_expression$sample_name))

#### Group the data by event
sostar_result_expression_events <- sostar_result_expression %>%
  ### Summarize the expression of the events
  dplyr::select(-c(transcript_id, chr, start, end, strand, gene, start, annot_ref, annot_find, new_annotation)) %>%
  dplyr::group_by(iso_event,sample_name) %>%
  dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(gene = strsplit(x = iso_event, split = "\\;")[[1]][1]) %>% 
  dplyr::ungroup() %>% 
  ### Relative expression of the events to the sample (normalizaton)
  dplyr::group_by(sample_name, gene) %>%
  dplyr::mutate(sample_event_expression = sum(event_expression),
                sample_FPKM = sum(FPKM),
                sample_TPM = sum(TPM)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(norm_event_expression = round(x = c(event_expression/sample_event_expression), digits = 4)*100,
                norm_FPKM = round(x = c(FPKM/sample_FPKM), digits = 4)*100,
                norm_TPM = round(x = (TPM/sample_TPM), digits = 4)*100) %>% 
  dplyr::select(-c(sample_event_expression,sample_FPKM,sample_TPM)) %>% 
  ### Relative expression of the events to the sample - percentile
  dplyr::group_by(sample_name, gene) %>%
  dplyr::mutate(sample_per_event_expression = percent_rank(event_expression),
                sample_per_FPKM = percent_rank(FPKM),
                sample_per_TPM = percent_rank(TPM)) %>% 
  dplyr::ungroup() %>% 
  ### Relative expression of the events to the event - percentile
  dplyr::group_by(iso_event, gene) %>%
  dplyr::mutate(iso_event_per_event_expression = percent_rank(event_expression),
                iso_event_per_FPKM = percent_rank(FPKM),
                iso_event_per_TPM = percent_rank(TPM)) %>% 
  dplyr::ungroup() %>%
  ### Do the FoldChange
  dplyr::group_by(iso_event) %>% 
  dplyr::mutate(mean_iso_event_event_expression = mean(event_expression),
                mean_iso_FPKM = mean(FPKM),
                mean_iso_TPM = mean(TPM)) %>% 
  dplyr::ungroup() %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(logFC_iso_event = log2(event_expression/mean_iso_event_event_expression),
                logFC_FPKM = log2(FPKM/mean_iso_FPKM),
                logFC_TPM = log2(TPM/mean_iso_TPM)) %>% 
  dplyr::ungroup() %>% 
  ### Check for the nª of occurrences
  dplyr::group_by(iso_event) %>%
  dplyr::mutate(occurrences = paste0(x = c(n(),"/",n_of_samples), collapse = "")) %>%
  dplyr::ungroup() %>% 
  ### Better annotation
  dplyr::rowwise() %>% 
  dplyr::mutate(gene = strsplit(x = iso_event, split = "\\;")[[1]][1],
                iso_event = strsplit(x = iso_event, split = "\\;")[[1]][2]) %>%
  dplyr::ungroup() %>% 
  dplyr::relocate(gene, .before = 1) %>% 
  ### Remove junk things
  dplyr::select(-c(norm_FPKM,sample_per_FPKM,iso_event_per_FPKM,logFC_FPKM))  
  # dplyr::select(-c(mean_iso_event_per_event_expression,mean_iso_FPKM,mean_iso_TPM)) %>%

sostar_result_expression_events_list_w_density <- list()
##### Check a la densitat
for (i in 1:length(unique(sostar_result_expression_events$sample_name))) {
  ### data now sample
  sample_now <- unique(sostar_result_expression_events$sample_name)[i]
  sostar_result_expression_events_sample <- sostar_result_expression_events %>% 
    dplyr::filter(sample_name == sample_now)
  
  for (g in 1:length(unique(sostar_result_expression_events$gene))) {
  ### data now sample gene
  gene_now <- unique(sostar_result_expression_events_sample$gene)[g]
  sostar_result_expression_events_sample_gene <- sostar_result_expression_events_sample %>% 
    dplyr::filter(gene == gene_now)
  
  ### do the density calculation
  density_now <- density(sostar_result_expression_events_sample_gene$event_expression)
  ## check peak density data and the coordinates
  peak_idx <- which.max(density_now$y)
  peak_x <- density_now$x[peak_idx]
  peak_y <- density_now$y[peak_idx]
  ## medium half
  half_y <- peak_y/1.5
  right_point <- seq(peak_idx, length(density_now$y))
  cross_idx <- right_point[which(density_now$y[right_point] <= half_y)[1]]
  cross_x <- density_now$x[cross_idx]
  cross_y <- density_now$y[cross_idx]
  
  ### plot the density
  png(filename = paste0(density_plot_place,sample_now,"_",gene_now,"_","coverage_density.png"), 
      height = 1000, width = 1000)
  plot(density_now, main = paste0("Density plot for event coverage\n",sample_now," gene: ",gene_now,"\n",
                                  "Covrage threshold:", as.character(round(cross_x, digits = 2))))
  
  abline(h = peak_y, col = "red")
  abline(v = peak_x, col = "red")
  
  abline(h = half_y, col = "blue", lty = 2)
  abline(v = cross_x, col = "blue", lty = 2)
  
  dev.off()
  
  ### mark the density
  sostar_result_expression_events_sample_gene <- sostar_result_expression_events_sample_gene %>% 
    dplyr::mutate(event_expression_density = ifelse(event_expression >= cross_x,"right","left"))
  
  sostar_result_expression_events_sample_gene$thress <- cross_x
  
  sostar_result_expression_events_list_w_density[[paste0(sample_now,"_",gene_now)]] <- sostar_result_expression_events_sample_gene
  }
  }

sostar_result_expression_events <- do.call(rbind, sostar_result_expression_events_list_w_density)

sostar_result_mod_events_long_list <- list()
sostar_result_mod_events_long_list_raw <- list()
for (i in 1:length(unique(sostar_result_expression_events$gene))) {
  gene_now <- unique(sostar_result_expression_events$gene)[i]
  
  
  data_now_raw <- sostar_result_expression_events %>%
    dplyr::filter(gene == gene_now)
  cross_x_now <- data_now_raw$thress 
    
  data_now <- sostar_result_expression_events %>%
    dplyr::filter(gene == gene_now) %>% 
    dplyr::filter(event_expression >= cross_x_now) %>%
    dplyr::ungroup()
  

  sostar_result_mod_events_long_list[[gene_now]] <- data_now
  sostar_result_mod_events_long_list_raw[[gene_now]] <- data_now_raw
  
  remove(gene_now);remove(data_now)
}

openxlsx::write.xlsx(x = sostar_result_mod_events_long_list,
                     file = event_data)
openxlsx::write.xlsx(x = sostar_result_mod_events_long_list_raw,
                     file = event_data_raw)


#### Event work filter
sostar_result_mod_events_long_list_check <- list()
sostar_result_mod_events_long_list_filt <- list()
for(i in 1:length(sostar_result_mod_events_long_list)) {
  ### data now
  gene_now <- names(sostar_result_mod_events_long_list)[i]
  data_now <- sostar_result_mod_events_long_list[[gene_now]]
  
  ### apply the filters
  data_now <- data_now %>% 
    ### select good cols
    dplyr::select(gene, iso_event, sample_name,
                  event_expression, FPKM, TPM,
                  sample_per_TPM, iso_event_per_TPM, logFC_TPM, occurrences, 
                  event_expression_density) %>% 
    ## do the filter (actually a check)
    dplyr::rowwise() %>% 
    dplyr::mutate(sample_per_check = ifelse(sample_per_TPM >= 0.66,"yes","no"),
                  iso_event_per_check = ifelse((iso_event_per_TPM >= 0.66 | is.na(iso_event_per_TPM)),"yes","no"),
                  logFC_check = ifelse(logFC_TPM > 0,"yes","no"),
                  occurrences_ratio = as.numeric(strsplit(x = occurrences, split = "\\/")[[1]][1])/as.numeric(strsplit(x = occurrences, split = "\\/")[[1]][2]),
                  occurrences_check = ifelse(occurrences_ratio <= 0.15,
                                       "yes","no"),
                  event_expression_density_check = ifelse(event_expression_density == "right",
                                                          "yes","no")) %>% 
    dplyr::ungroup() %>%
    dplyr::ungroup()
  
  ### filter for those isoforms that meet 3 out of 5 
  data_now$summ_of_filters <- rowSums(data_now[, grep(pattern = "check", x = colnames(data_now))] == "yes", na.rm = TRUE)
  sostar_result_mod_events_long_list_check[[gene_now]] <- data_now
  data_now <- data_now %>%
    dplyr::filter((summ_of_filters >= 3 & event_expression > 10) | 
                    (occurrences_ratio <= 0.15) | 
                    (summ_of_filters >= 3 & event_expression_density_check == "yes"))
  
  sostar_result_mod_events_long_list_filt[[gene_now]] <- data_now
  }

openxlsx::write.xlsx(x = sostar_result_mod_events_long_list_filt, 
                     file = event_data_filt)

openxlsx::write.xlsx(x = sostar_result_mod_events_long_list_check,
                     file = event_data_check)

