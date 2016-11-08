#! /bin/bash

#SBATCH --qos long
#SBATCH -A stsn
#SBATCH -p Lewis
#SBATCH -J psorunskriging
#SBATCH -o psoruns.o
#SBATCH -e psoruns.e
#SBATCH --mem-per-cpu=3G
#SBATCH -t 7-00:00
#SBATCH -n 28  # cores
#SBATCH -N 1  # nodes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=themattsimpson@gmail.com

VERSION="3.2.3"
module load R/R-${VERSION}

R --no-restore --no-save CMD BATCH psoruns.R
