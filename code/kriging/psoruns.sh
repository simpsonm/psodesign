#! /bin/bash

#SBATCH -A stsn
#SBATCH -p Lewis
#SBATCH -J psorunskriging
#SBATCH -o psoruns.o
#SBATCH -e psoruns.e
#SBATCH --mem 10000
#SBATCH -t 7-00:00
#SBATCH -n 20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=themattsimpson@gmail.com

VERSION="3.2.3"
module load R/R-${VERSION}

R --no-restore --no-save CMD BATCH psoruns.R
