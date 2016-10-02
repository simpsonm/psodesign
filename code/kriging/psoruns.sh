#! /bin/bash

#SBATCH -A stsn
#SBATCH -J psorunskriging
#SBATCH -o psoruns.o
#SBATCH -e psoruns.e
#SBATCH --mem 2000
#SBATCH -t 5-00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=themattsimpson@gmail.com

R --no-restore --no-save CMD BATCH psoruns.R
