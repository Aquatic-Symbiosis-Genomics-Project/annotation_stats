#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("ggplot2")
library("gridExtra")
library("ggpubr")

intron_lens<-args[1]
exon_lens<-args[2]
gene_lens<-args[3]
perc_support<-args[4]

#Plot intron length dist.
intL<-data.frame(read.table(intron_lens, header = FALSE))
attach(intL)
Log_intLen<-ggplot() + geom_histogram(aes(log(intL$V1)), bins = 1000, fill = "blue") + geom_density(alpha = 0.4) + theme_classic() + labs(x ="Intron length(bp)", y = "Count") #title="Distribution of intron lengths (log)"
Intron_len<-ggplot() + geom_histogram(aes(intL$V1), bins = 1000, fill = "blue") + geom_density(alpha = 0.4) + theme_classic() + labs(x ="Intron length(bp)", y = "Count") # title="Distribution of intron lengths"
detach(intL)

# Plot exon length dist.
exL<-data.frame(read.table(exon_lens, header = FALSE))
attach(exL)
Log_exonLen<-ggplot() + geom_histogram(aes(log(exL$V1)), bins = 1000, fill = "blue") + geom_density(alpha = 0.4) + theme_classic() + labs(x ="Exon length(bp)", y = "Count")#title="Distribution of exon lengths (log)"
Exon_len<-ggplot() + geom_histogram(aes(exL$V1), bins = 1000, fill = "blue") + geom_density(alpha = 0.4) + theme_classic() + labs(x ="Exon length(bp)", y = "Count") #title="Distribution of exon lengths"
detach(exL)

#Plot gene length dist.
genL<-data.frame(read.table(gene_lens, header = FALSE))

attach(genL)
Log_geneLen<-ggplot() + geom_histogram(aes(log(genL$V1)), bins = 1000, fill = "blue") + geom_density(alpha = 0.4) + theme_classic() + labs(x ="Gene length(bp)", y = "Count") #title="Distribution of gene lengths (log)"
Gene_len<-ggplot() + geom_histogram(aes(genL$V1), bins = 1000, fill = "blue") + geom_density(alpha = 0.4) + theme_classic() + labs(x ="Gene length(bp)", y = "Count") #title="Distribution of gene lengths"
detach(genL)

#Plot transcript support from augustus.
perc<-data.frame(read.table(perc_support, header = FALSE))
attach(perc)

#Plot all support.
all_support<-ggplot() + geom_histogram(aes(perc$V3), fill = "blue", bins = 100) + geom_density(alpha = 0.4) + theme_classic() + labs(x ="% support", y = "Count")#title="Distribution of transcript support (all sources) ",

#Plot intron support.
intron_support<-ggplot() + geom_histogram(aes(perc$V4), fill = "blue", bins = 100) + geom_density(alpha = 0.4) + theme_classic() + labs(x ="% support", y = "Count")#title="Distribution of transcript support - Introns ",

# arrange plots
res <- marrangeGrob(list(Gene_len, Exon_len, Log_geneLen, Log_exonLen, Intron_len, all_support, Log_intLen, intron_support), nrow = 2, ncol = 2)

# Export summary stat plots pdf file
ggexport(res, filename = "Summary_stats.pdf")

# subset only localised scaffolds based on SUPER_id 
chroms_subset<-perc[grep("SUPER_", perc$V1),]

# removed unlocalised (with SUPER_id)
chroms2<-data.frame(chroms_subset[!grepl("unloc",chroms_subset$V1), ])

#Get unique list of chromosome ids
chroms3<-unique(chroms2[1])

# PLot perc support
chroms_support_all<-ggplot(chroms2, aes(chroms2$V3)) + geom_histogram(fill = "blue", bins = 30) + geom_density(alpha = 0.4) + theme_classic() + labs(title="Transcript support (all sources) ", x ="% support", y = "Count") + facet_wrap(~ chroms2$V1)

# Plot intron support 
chroms_support_intron<-ggplot(chroms2, aes(chroms2$V4)) + geom_histogram(fill = "blue", bins = 30) + geom_density(alpha = 0.4) + theme_classic() + labs(title="Transcript support (introns) ", x ="% support", y = "Count") + facet_wrap(~ chroms2$V1)

# arrange plots
res2 <- marrangeGrob(list(chroms_support_all, chroms_support_intron), nrow = 1, ncol = 1)

# Export chromosome stats to pdf file
ggexport(res2, filename = "Summary_stats_chromosomes.pdf")







