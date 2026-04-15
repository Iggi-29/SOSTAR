######################################
### Make the sostar Run per sample ###
######################################

### libraries
library(yaml)
library(tidyverse)
library(dplyr)

### pick the files
config_file <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR/config/config.yaml"
snake_file <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR/workflow/Snakefile"
sh_file <- "/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR/workflow/Sostar_in_cluster.sh"

config_content <- yaml::read_yaml(file = config_file)
for (i in 1:length(config_content$samples)) {
# for (i in 1:length(1)) {
  ### config_generation
  config_content_now <- config_content
  sample_now <- config_content$samples[i]
  
  config_content_now$samples <- config_content_now$samples[config_content_now$samples == sample_now]
  config_content_now$dirs$outdir <- paste0("/imppc/labs/eclab/ijarne/0_Recerca/pipelines/sostar/SOSTAR_output_EN_all_samples_individual_",
                                           sample_now)
  
  config_content_file_now <- paste0(x = c(
    gsub(x = config_file,
         pattern = "\\/config\\.yaml$",
         replacement = ""),
    "/",
    "config_",
    sample_now,
    ".yaml"),
    collapse = "")
  
  yaml::write_yaml(x = config_content_now, 
                   file = config_content_file_now)
  
  ### Generate the snakefile
  snake_file_content <- readLines(con = snake_file, warn = FALSE)
  snake_file_content[3] <- paste0(x = c(
    'configfile: \"../config/',
    basename(config_content_file_now),'"'), 
    collapse = "")
  
  snake_file_now <- paste0(x = c(
    gsub(x = snake_file,
         pattern = "\\/Snakefile$",
         replacement = ""),
    "/",
    "Snakefile_",
    sample_now,
    ".smk"
  ),
  collapse = "")
  writeLines(text = snake_file_content, con = snake_file_now)
  
  ### Generate the sh file
  sh_file_content <- readLines(con = sh_file, warn = FALSE)
  sh_file_content[51] <- paste0(x = c(
    "--rerun-incomplet"," -s ",basename(snake_file_now)),
    collapse = "")
  sh_file_content[12] <- paste0(x = c("#$ -e SOSTAR_EN_individual_",sample_now,".err"), collapse = "")
  sh_file_content[14] <- paste0(x = c("#$ -o SOSTAR_EN_individual_",sample_now,".log"), collapse = "")
  sh_file_content[15] <- "#$ -pe smp 4"                                                                            
  sh_file_content[48] <- "snakemake -j4 -p --skip-script-cleanup --resources mem_mb=160000 \\"                     
  
  
  sh_file_now <- paste0(x = c(
  gsub(x = sh_file,
       pattern = "\\/Sostar_in_cluster\\.sh$",
       replacement = ""),
  "/",
  "Sostar_in_cluster_",
  sample_now,
  ".sh"), 
  collapse = "")
  
  writeLines(text = sh_file_content, con = sh_file_now)
  
  ### execute the code!
  cmd_qsub <- paste0("qsub ",sh_file_now)
  system(command = cmd_qsub)
}



