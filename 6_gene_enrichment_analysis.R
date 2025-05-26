library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)

### annotation format and genes informations
gene2loc <- read.csv("N:/3_P3-Projets/1_PROJETS AVEC RECETTES/BLUEREMEDIOMICS/SCIENTIFIQUE/Stage Florian Petrilli/3_analyses/transcripto/analyse_degs_banc/SC/fnl_analysis/geneid_to_locus.tsv", sep="\t", row.names=1)
gene2loc$GeneID <- rownames(gene2loc)
eggnog <- read.csv("N:/3_P3-Projets/1_PROJETS AVEC RECETTES/BLUEREMEDIOMICS/SCIENTIFIQUE/Stage Florian Petrilli/3_analyses/transcripto/analyse_degs_banc/SC/fnl_analysis/all_genomes_eggnog.tsv", sep="\t")
eggnog <- eggnog %>% rename(LocusTag=Label) 

gen_an <- merge(gene2loc,eggnog,by="LocusTag")
#write.csv(gen_an, "gene_anotated.csv", row.names = FALSE)


df <- DESeq_results
df$GeneID <- rownames(DESeq_results)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GeneID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

degs_an <- merge(df,gen_an, by = "GeneID")


### functionnal enrichment

signif_degs <- degs_an %>% filter(padj < 0.05)

gene_list <- signif_degs$LocusTag

TERM2GENE <- degs_an %>% select(COGclassID, LocusTag) %>% filter(COGclassID != "") %>% distinct()
TERM2GENE <- TERM2GENE %>% rename(gene=LocusTag ,term=COGclassID)

TERM2NAME <- gen_an %>% select(COGclassID, COGClassDescription) %>% filter(COGclassID != "") %>% unique()

universe <- unique(TERM2GENE$gene)

enrich_result <- enricher(
  gene = gene_list,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  universe = universe,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

lfc <- -degs_an$log2FoldChange
names(lfc) <- degs_an$LocusTag

# centplot visualisation

cnet_nvb <- cnetplot(enrich_result, categorySize="p.ajust", foldChange = lfc, circular = FALSE, colorEdge = TRUE, node_label = 'category')
cnet_pvb <- cnetplot(enrich_result, categorySize="p.ajust", foldChange = lfc, circular = FALSE, colorEdge = TRUE, node_label = 'category')
cnet_nvp <- cnetplot(enrich_result, categorySize="p.ajust", foldChange = lfc, circular = FALSE, colorEdge = TRUE, node_label = 'category')

plot_grid(cnet_nvb, cnet_pvb, cnet_nvp, labels=c("A", "B", "C"), ncol = 3)
