intron_addition <- function(gtf_df) {
  
  gtf_df_to_bind <- gtf_df %>% 
    dplyr::filter(type != "transcript")
  
  list_of_intron_additions <- list()
  
  for (i in 1:length(unique(gtf_df$transcript_id))) {
    
    transcript_now <- unique(gtf_df$transcript_id)[i]
    
    gtf_df_now <- gtf_df %>% 
      dplyr::filter(transcript_id == transcript_now) %>% 
      dplyr::mutate(exon_number = as.numeric(exon_number)) %>% 
      dplyr::filter(type != "transcript")
    
    strand_now <- unique(as.character(gtf_df_now$strand))
    
    # Skip transcripts with ambiguous/missing strand or only 1 exon
    exons_now <- gtf_df_now %>% dplyr::filter(type == "exon")
    if (!strand_now %in% c("+", "-") || nrow(exons_now) < 2) next
    
    if (strand_now == "+") {
      intron_additions <- exons_now %>% 
        dplyr::arrange(start) %>% 
        dplyr::mutate(
          next_exon_start = dplyr::lead(start),
          intron_number   = exon_number,
          intron_start    = end + 1,
          intron_end      = next_exon_start - 1
        ) %>% 
        dplyr::filter(!is.na(next_exon_start))
      
    } else {
      intron_additions <- exons_now %>% 
        dplyr::arrange(dplyr::desc(start)) %>% 
        dplyr::mutate(
          next_exon_end = dplyr::lead(end),
          intron_number = exon_number,
          intron_start  = next_exon_end + 1,
          intron_end    = start - 1
        ) %>% 
        dplyr::filter(!is.na(next_exon_end))
    }
    
    # Compute width AFTER resolving intron_start/intron_end
    intron_additions <- intron_additions %>%
      dplyr::mutate(width_new = intron_end - intron_start + 1) %>%
      dplyr::transmute(
        seqnames      = seqnames,
        start         = intron_start,
        end           = intron_end,
        width         = width_new,
        strand        = strand,
        type          = "intron",
        gene_id       = gene_id,
        transcript_id = transcript_id,
        exon_id       = exon_id,
        exon_number   = intron_number
      )
    
    list_of_intron_additions[[i]] <- intron_additions
  }
  
  # Remove NULL entries (skipped transcripts)
  list_of_intron_additions <- Filter(Negate(is.null), list_of_intron_additions)
  
  list_of_intron_additions <- do.call(rbind, list_of_intron_additions)
  list_of_intron_additions <- list_of_intron_additions %>% 
    dplyr::select(-width)
  
  gtf_w_introns <- rbind(list_of_intron_additions, gtf_df_to_bind) %>%
    dplyr::arrange(gene_id, transcript_id, start) %>% 
    dplyr::mutate(exon_number = as.numeric(exon_number))
  
  return(gtf_w_introns)
}