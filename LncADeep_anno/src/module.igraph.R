rm(list=ls())

library(igraph)

args <- commandArgs(trailingOnly = TRUE)
interactingPairs <- args[1]     # interacting protein pairs of module
out_figure <- args[2]       # output figure
#interactingPairs <- "/data/cyang/analysis/dnn_LPDIL_test/gencode.v24.lncRNA-Uniprot.protein_pairs_results/pair.part.000000_LPIDL_Test/test.igraph.ppi.0"
#out_figure <- "test.pdf"

el <- read.table(interactingPairs)
el <- as.matrix(el)

g <- graph.edgelist(el[,1:2], directed = F)
g <- simplify(g, remove.loops=T)

E(g)$weight <- as.numeric(el[,3])

pdf(out_figure,width=8,height=6)
#l <- layout.circle(g)
#l <- layout.fruchterman.reingold(g)
l <- layout.kamada.kawai(g)
plot(g, layout = l, vertex.size = 6.18, edge.color="orange", vertex.frame.color="#ffffff", vertex.label.cex=0.382, edge.width=E(g)$weight)
dev.off()
