##########################
### Sashimi plots in R ###
##########################

##### Libraries to work on
library(biomaRt)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(GenomicAlignments)
library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)

##### Function sourcing
source("/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR/workflow/functions/intron_addition.R")

##### Constants
bam_file_place <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR_output_IGTP_all_samples/alignment/realigned"
annotation_of_genes <- "/imppc/labs/eclab/Resources/MINION_ddbb/Ref_files/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf"
trans_to_work_on <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR_reffernce/trans_list.txt"

##### Imort annotation data
### Interesting transcripts
trans_to_work_on <- readLines(con = trans_to_work_on)
### GTF annotation of the genes
annotation_of_genes <- rtracklayer::import(con = annotation_of_genes)
## GTF annotation to dataframe
annotation_of_genes_df <- as.data.frame(annotation_of_genes)
# transcript information
annotation_of_genes_df_trans <- annotation_of_genes_df %>% 
  dplyr::filter(type == "transcript") %>% 
  dplyr::filter(transcript_id %in% trans_to_work_on)
# exon and intron information
annotation_of_genes_df <- annotation_of_genes_df %>% 
  dplyr::filter(type %in% c("exon","CDS")) %>% 
  dplyr::filter(transcript_id %in% trans_to_work_on)

annotation_of_genes_df <- annotation_of_genes_df %>% 
  dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id, exon_id, exon_number) %>% 
  dplyr::mutate(exon_number = as.numeric(exon_number))

## filter non-useful CDS
annotation_of_genes_df <- annotation_of_genes_df %>% 
  dplyr::group_by(transcript_id) %>% 
  dplyr::filter(type != "CDS" | exon_number == max(exon_number) | exon_number == min(exon_number))  %>% 
  dplyr::ungroup()

annotation_of_genes_df <- intron_addition(gtf_df = annotation_of_genes_df)

#### Read the alignments
list_of_bams <- list.files(path = bam_file_place, pattern = "\\.bam$", full.names = TRUE)
list_of_bams <- list_of_bams[1]
### Define regions of interest
regions_of_interest <- GRanges(seqnames = as.character(annotation_of_genes_df_trans$seqnames),
                               ranges = IRanges(
                                 start = annotation_of_genes_df_trans$start,
                                 end = annotation_of_genes_df_trans$end))
regions_of_interest <- GenomicRanges::reduce(x = regions_of_interest)
### Parameters to scan the bam file
param <- ScanBamParam(
  what = c("qname", "flag", "rname", "pos", "cigar", "strand"),
  which = regions_of_interest)

for (i in 1:length(list_of_bams)) {
  ### Data now
  bam_file_now <- list_of_bams[i]
  sample_name_now <- gsub(pattern = "\\.realigned\\.bam$", replacement = "", x = basename(bam_file_now))
  cat(paste0(x = c("Working on sample: ",i," out of ",length(list_of_bams)," which is: ",sample_name_now,"\n"), collapse = ""))
  ### Import the data
  cat(paste0(x = c("Getting the alignment data","\n"), collapse = ""))
  gal <- readGAlignments(file = bam_file_now, param = param)
  
  ### Raw junction data
  cat(paste0(x = c("Getting the junctions data","\n"), collapse = ""))
  juncs <- summarizeJunctions(x = gal)
  juncs_df <- as.data.frame(juncs)
  juncs_df <- juncs_df %>%
    dplyr::mutate(junction = paste0(seqnames,":",start,"-",end))
  ## Junction that each read has
  juncs_by_read <- cigarRangesAlongReferenceSpace(
    cigar = cigar(gal),ops = "N",
    reduce.ranges = FALSE, drop.empty.ranges = TRUE)

  n_reads <- length(juncs_by_read)
  tictoc::tic()
  grl_juncs <- GRangesList(
    lapply(X = seq_along(juncs_by_read), FUN = function(i) {
      ir <- juncs_by_read[[i]]
      cat("Processing read",i,"/",n_reads,"\r")

      if (length(ir) == 0) return(GRanges())

      GRanges(
        seqnames = seqnames(gal)[i],
        ranges = shift(ir, start(gal)[i] -1),
        strand = strand(gal)[i],
        read_id = mcols(gal)$qname[i]
      )

    })
  )
  tictoc::toc()

  all_juncs <- unlist(grl_juncs, use.names = FALSE)
  junction_table <- as.data.frame(all_juncs) %>%
    mutate(junction = paste0(seqnames, ":", start, "-", end)) %>%
    group_by(junction) %>%
    summarise(
      reads = list(read_id),
      n_reads = n()
    )
  remove(all_juncs);remove(juncs);remove(grl_juncs);remove(juncs_by_read)

  ### All data juncs
  all_data_juncs <- merge(juncs_df, junction_table, by = "junction")
  remove(juncs_df);remove(junction_table)
  }




