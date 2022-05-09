#!/usr/bin/env Rscript

### Exploring the repeat content of Marasmius oreades
#############################################################################
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-10-30
# Version 1
# =======================================

library(ggplot2)
library(dplyr)

# ============================
# Functions
# ============================
readgff <- function(gff_file){
  gff <- read.table(gff_file) # a gtf actually
  names(gff) <- c("sequence", "source", "feature", "start", "end", "identity", "strand", "phase", "Target", "Motif", "start_cons", "end_cons")
  # Remove some irrelevant columns
  gff <- select(gff, -c("source", "feature", "phase", "Target"))
  
  # Add length of alignment
  gff <- gff %>% mutate(len_ref = end - start + 1) %>%
    mutate(Motif= gsub("Motif:", "", Motif)) %>% 
    filter(!grepl(")n", Motif)) %>% filter(!grepl("-rich", Motif))
  
  # Make it a factor, not a character
  gff$Motif <- factor(gff$Motif)
  
  return(gff)
}

# ============================
# Data
# ============================
gff_file <- snakemake@input$gtf

# Read gff
gff <- readgff(gff_file)

# ============================
# Plot
# ============================

## Make a table with counts of occurance per Motif in the gff
motifcounts <- gff$Motif %>% table() %>% data.frame
names(motifcounts) <- c("Motif", "Count")
# Reorder by abundance
motifcounts <- motifcounts[order(motifcounts$Count, decreasing = TRUE),]

write.table(motifcounts, file = snakemake@output$tbl, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# How common are the motifs
TEhist <- ggplot(motifcounts, aes(Count)) + geom_histogram() +
  xlab("Count of repeat copies") +
  ggtitle(snakemake@wildcards$sample) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = TEhist, snakemake@output$hist, width = 3, height = 4)
