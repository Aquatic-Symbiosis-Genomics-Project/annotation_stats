#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("ggplot2")
library("gridExtra")
library("ggpubr")


perc_support<-args[1]

#Read in perc support from augustus.
perc<-data.frame(read.table(perc_support, header = FALSE))
attach(perc)

#Plot all support/ note this includes all sources e.g exons (which will be zero if only using rna-seq input) 
#all_support<-ggplot() + geom_histogram(aes(perc$V3), fill = "blue", bins = 100) + geom_density(alpha = 0.4) + theme_classic() + labs(x ="% support", y = "Count")#title="Distribution of transcript support (all sources) ",

# Export summary stat plots pdf file
ggexport(res, filename = "Summary_stats.pdf")

# subset only localised scaffolds based on SUPER_id 
chroms_subset<-perc[grep("SUPER_", perc$V1),]

# removed unlocalised (with SUPER_id)
chroms2<-data.frame(chroms_subset[!grepl("unloc",chroms_subset$V1), ])

#Get unique list of chromosome ids
chroms3<-unique(chroms2[1])

# Plot intron support 
chroms_support_intron<-ggplot(chroms2, aes(chroms2$V4)) + geom_histogram(fill = "blue", bins = 30) + geom_density(alpha = 0.4) + theme_classic() + labs(title="Transcript support (introns) ", x ="% support", y = "Count") + facet_wrap(~ chroms2$V1)

# Export chromosome stats to pdf file
ggexport(res2, filename = "Summary_stats_chromosomes.pdf")







