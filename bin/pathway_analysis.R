#!/usr/bin/env Rscript

library(org.Hs.eg.db)
library(fgsea)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
deseq_results <- args[1] # deseq-results-tidy.csv
hallmark      <- args[2] # h.all.v6.2.symbols.gmt
kegg          <- args[3] # c2.cp.kegg.v6.2.symbols.gmt
mir           <- args[4] # c3.mir.v6.2.symbols.gmt
go            <- args[5] # c5.all.v6.2.symbols.gmt

res <- read.csv(deseq_results, 
                header           = TRUE, 
                stringsAsFactors = FALSE, 
                check.names      = FALSE);

# Map Ensembl gene IDs to symbol. First create a mapping table.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)

# join them
res <- inner_join(res, ens2symbol, by=c("row"="ENSEMBL"))

# Further, all you’ll care about later on is the gene symbol and the test statistic. Get just those, and remove the NAs. Finally, if you have multiple test statistics for the same symbol, you’ll want to deal with that in some way. Here I’m just averaging them.
res2 <- res %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

# Using the fgsea package
ranks <- deframe(res2)

# Load the pathways into a named list
pathways.hallmark <- gmtPathways(hallmark)

# Run the fgsea algorithm with 1000 permutations:
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

# Tidy the results
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Plot the normalized enrichment scores. Color the bar indicating whether or not the pathway was significant:
# TODO: save plot
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
ggsave('hallmark_pathways_gsea.png')

# What genes are in each of these pathways? First, get a tibble with all the pathways and the genes in them. Continue to join that back to the original data to pull out genes in the pathways. Optionally, filter the list to include only those that are significant, etc. 
hallmark_pathways <- pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest() %>% 
  inner_join(res, by="SYMBOL")

# Let's try a different set of pathways. Let's look at KEGG pathways
kegg_pathways <- fgsea(pathways=gmtPathways(kegg), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)

# Or miR targets
mir_pathways <- fgsea(pathways=gmtPathways(mir), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)

# Or GO annotations
go_pathways <- fgsea(pathways=gmtPathways(go), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)

# TODO: save the following:
# fgseaResTidy
# hallmark_pathways
# kegg_pathways
# mir_pathways
# go_pathways