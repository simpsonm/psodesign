#! /bin/bash

#SBATCH -A stsn
#SBATCH -J accpsorunselection
#SBATCH -o accpsoruns.o
#SBATCH -e accpsoruns.e
#SBATCH --mem 10000
#SBATCH -t 5-00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=themattsimpson@gmail.com

R --no-restore --no-save CMD BATCH accpsoruns.R
