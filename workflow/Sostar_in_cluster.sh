#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
# Execute from the current workig dir
#$ -cwd
# Name for the script in the queuing system
#$ -N SOSTAR_IGTP
# In order to load the environment variables and your path
# You can either use this or do a : source /etc/profile
#$ -V
# You can redirect the error output to a specific file
#$ -e SOSTAR_IGTP.err
# You can redirect the output to a specific file
#$ -o SOSTAR_IGTP.log
#$ -pe smp 10
#$ -q d10imppcv3  # los nodos nuevos!
#$ -l h_vmem=9G
# Avoid node 16 (low memory at the moment)
# -l h=!sge-exec-16

# Get the enviroment name
HOST=`hostname -s`

echo the parameters
echo "################################################################################"
echo ""
echo "Run at: $HOST"
date
echo "###############################################################################"
echo ""

# Start working
echo "Starting Analysis Pipeline..."
echo "Nanopore Analysis Pipeline started working at: $(date)"

# Get Variables
WORKDIR=$(pwd)
Myuser=$(echo $USER)
echo "Running at dir: $WORKDIR"

# Activate anacodna
source "/imppc/labs/eclab/${Myuser}/miniconda3/etc/profile.d/conda.sh"
conda activate snakemake

# Snakemake run
snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
snakemake --unlock
snakemake -j20 -p --skip-script-cleanup --resources mem_mb=160000 \
--use-singularity --singularity-args "--cleanenv --bind /imppc/labs/eclab/" \
--use-conda --conda-frontend conda \
--rerun-incomplet

echo "Done"
echo "Nanopore Analysis Pipeline finished at: $(date)"