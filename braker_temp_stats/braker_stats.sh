#!/bin/bash
# set -e

gff_file="$1"
gtf_file="$2"
id="$3"


python braker_stats.py "$gff_file"  "$gtf_file" "$id" > "$id"_summary.stats

./braker_plots.R  "$id"_perc_support.txt



#Notes on usage:
# You need some packages in R...do:  Install.packages("ggplot2", "gridExtra", "ggpubr")
# if you want to change id's to get per chr support change the "SUPER_" here:
# chroms_subset<-perc[grep("SUPER_", perc$V1),]
