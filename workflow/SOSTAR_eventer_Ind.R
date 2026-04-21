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
  make_option(opt_str = "--SOSTAR.eventer.folder", 
              help = "Folder with the final data of sostar", dest = "sostar_eventer_folder"),
  make_option(opt_str = "--Refference.trans", 
              help = "Refference transcripts parsed to Sostar", dest = "trans"),
  make_option(opt_str = "--Individual.results",
              help = "Individual run results", dest = "ind")
)
opt <- parse_args(OptionParser(option_list = option_list))

#### Warnings
if (is.null(opt$sostar_eventer_folder) | is.null(opt$trans) | is.null(opt$ind)) {
  stop("Check for help")
}

#### Final data folder
sostar_eventer_folder <- opt$sostar_eventer_folder
refference_ids <- opt$trans
ind <- opt$ind

# sostar_eventer_folder <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR_output_EN/SOSTAR_Events/" 
# refference_ids <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR_reffernce/trans_list_ICO.txt"
# ind <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR_output_EN/Individual/"

#### Raw data import
refference_ids <- readLines(con = refference_ids)

#### Data files for batch runs
event_data_file <- paste0(x =c(sostar_eventer_folder,"Event_data.xlsx"), collapse = "")
event_data_raw_file <- paste0(x =c(sostar_eventer_folder,"Event_data_raw.xlsx"), collapse = "")
event_data_filt_file <- paste0(x =c(sostar_eventer_folder,"Event_data_filt.xlsx"), collapse = "")
event_data_check_file <- paste0(x =c(sostar_eventer_folder,"Event_data_check.xlsx"), collapse = "")

genes_to_work_on <- readxl::excel_sheets(path = event_data_raw_file)

event_data <- lapply(X = genes_to_work_on, FUN = function(x) {
  openxlsx::read.xlsx(xlsxFile = event_data_file, sheet = x)
  })
names(event_data) <- genes_to_work_on
event_data_raw <- lapply(X = genes_to_work_on, FUN = function(x) {
  openxlsx::read.xlsx(xlsxFile = event_data_raw_file, sheet = x)
})
names(event_data_raw) <- genes_to_work_on
event_data_check <- lapply(X = genes_to_work_on, FUN = function(x) {
  openxlsx::read.xlsx(xlsxFile = event_data_check_file, sheet = x)
})
names(event_data_check) <- genes_to_work_on
event_data_filt <- lapply(X = genes_to_work_on, FUN = function(x) {
  openxlsx::read.xlsx(xlsxFile = event_data_filt_file, sheet = x)
})
names(event_data_filt) <- genes_to_work_on

#### SOSTAR data - individual
sostar_individual_data_list <- list.files(path = ind,
                                          recursive = TRUE, full.names = TRUE,
                                          pattern = "SOSTAR_annotation_table_results.xlsx$")
event_of_indiviual_analysis <- list()
for (i in 1:length(sostar_individual_data_list)) {
  ### data now
  sostar_now <- sostar_individual_data_list[i]
  sample_name <- gsub(pattern = "SOSTAR_annotation_table_results\\.xlsx$",
                      replacement = "", x = sostar_now)
  sample_name <- basename(path = sample_name)
  sample_name <- gsub(pattern = ".*_individual_", replacement = "", x = sample_name)
  sostar_result <- openxlsx::read.xlsx(xlsxFile = sostar_now)
  
  ### Start the filtering
  sostar_result <- sostar_result %>% 
    dplyr::filter(occurence > 0)  
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
  names(sostar_result_mod_annotated_list) <- paste0(sample_name,"_",names(sostar_result_mod_annotated_list))
  event_of_indiviual_analysis <- c(event_of_indiviual_analysis,sostar_result_mod_annotated_list)
  }
remove(sample_name)
remove(sostar_result_mod)
remove(i)
remove(genes_to_work_on)
remove(refference_ids)
remove(sostar_result)
remove(sostar_result_mod_annotated)
remove(sostar_now)
remove(sostar_individual_data_list)

batch_data_annotator <- function(batch_to_annotate, event_of_indiviual_analysis = event_of_indiviual_analysis) {
  
  final_list <- list()
  for (i in 1:length(batch_to_annotate)) {
    ### data now
    gene_now <- names(batch_to_annotate)[i]
    event_gene_now <- batch_to_annotate[[i]]
    
    ### individual data now
    individual_data_now <- event_of_indiviual_analysis[grep(
      x = names(event_of_indiviual_analysis),
      pattern = paste0("_",gene_now,"$"))]
  
    data_of_the_genes <- data.frame()
    for (samp in 1:length(unique(event_gene_now$sample_name))) {
      ### data now
      sample_now <- unique(event_gene_now$sample_name)[samp]
      event_gene_now_sample_now <- event_gene_now %>%
        dplyr::filter(sample_name == sample_now)
      ### individual data filering
      individual_data_now_sample_now <- individual_data_now[[grep(pattern = paste0("^",sample_now,"_"),
                                                                  x = names(individual_data_now))]]
      individual_data_now_sample_now <- individual_data_now_sample_now$iso_event
      ## annotate data now
      event_gene_now_sample_now <- event_gene_now_sample_now %>%
        dplyr::mutate(supported_by_individual_analysis = ifelse(iso_event %in% individual_data_now_sample_now,"yes","no"))
      data_of_the_genes <- rbind(data_of_the_genes,event_gene_now_sample_now)
    }
    final_list[[gene_now]] <- data_of_the_genes
  }
  return(final_list)
}

event_data <- batch_data_annotator(batch_to_annotate = event_data, 
                                    event_of_indiviual_analysis = event_of_indiviual_analysis)
event_data_raw <- batch_data_annotator(batch_to_annotate = event_data_raw, 
                                    event_of_indiviual_analysis = event_of_indiviual_analysis)
event_data_filt <- batch_data_annotator(batch_to_annotate = event_data_filt, 
                                    event_of_indiviual_analysis = event_of_indiviual_analysis)
event_data_check <- batch_data_annotator(batch_to_annotate = event_data_check, 
                                    event_of_indiviual_analysis = event_of_indiviual_analysis)

openxlsx::write.xlsx(file = event_data_file, x = event_data)
openxlsx::write.xlsx(file = event_data_raw_file, x = event_data_raw)
openxlsx::write.xlsx(file = event_data_filt_file, x = event_data_filt)
openxlsx::write.xlsx(file = event_data_check_file,x = event_data_check)















