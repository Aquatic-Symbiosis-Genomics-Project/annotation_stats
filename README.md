# annotation_stats

1. Install R libraries (will need to be in personal library on cluster, lib not writable....)


Install.packages("ggplot2", "gridExtra", "ggpubr")


3. Generate braker stats:        
        python braker_stats.py gff_file  gtf_file ID
        
4. Plot stats: 
        ./braker_plots.R perc_support
        
        
Note chromosome id's need to be in format: e.g SUPER_1  and SUPER_1_unloc (to edit later on....)


