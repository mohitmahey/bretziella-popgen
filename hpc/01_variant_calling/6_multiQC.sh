#!/bin/bash --login

# running multiQC on final product

DIR="/mnt/gs21/scratch/maheymoh/karan_project_jgi"
cd ${DIR}

## loading multiqc
module purge
module load MultiQC/1.14-foss-2022b

multiqc . -n gatk_multiqc_complete

echo "done with analysis and produced final raw counts file"

