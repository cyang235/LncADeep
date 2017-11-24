rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

sampleFreqFile <- args[1]       # the input file of sample frequency
output_table <- args[2]         # the output file of KEGG pathway
#sampleFreqFile <- "pair.part.000000_ENST00000629393.1.ko.freq"

sample.df <- read.table(file=sampleFreqFile, sep ="\t", quote="", header = F)
sample.df$m <- sample.df$V3
sample.df$n <- sample.df$V5 - sample.df$V3
sample.df$k <- sample.df$V4
sample.df$x <- sample.df$V2 - 1

sample.df$p_value <- phyper(sample.df$x, sample.df$m, sample.df$n, sample.df$k, lower.tail = F)

## subset those pathway with at least 5 node and p_value < 0.05
sigSample <- subset(sample.df, V2 > 4 & p_value < 0.05)    

sigSample$adj_p_value <- p.adjust(sigSample$p_value, "BH")        ## for KEGG pathway
sigSample <- sigSample[ order(sigSample$adj_p_value), ]

output.df <- subset(sigSample[,c(1,6,11,12)], adj_p_value < 0.05)
colnames(output.df) <- c("KEGG_path_ID", "KEGG_pathway", "p_value", "adj_p_value") 
str(output.df)
#output_table = "tmp.results"
write.table(format(output.df, digits = 2, nsmall = 5), output_table, row.names = FALSE)
