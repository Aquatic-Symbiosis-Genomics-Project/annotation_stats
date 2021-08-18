# annotation_stats

1. Install R libraries (will need to be in personal library on cluster, lib not writable....)


Install.packages("ggplot2", "gridExtra", "ggpubr")


3. Generate braker stats:        
        python braker_stats.py gff_file  gtf_file ID
        
4. Plot stats: 
        ./braker_plots.R intron_lens exon_lens gene_lens perc_support


