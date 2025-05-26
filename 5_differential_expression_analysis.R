library(Rsubread)
library(DESeq2)
library(tidyverse)
library(airway)
library(dplyr)
library(apeglm)
library(cowplot)

#############################################
######### données et metadonnées ############
#############################################

colData <- read.csv("colData.csv", sep = ",", row.names = 1)

# gene id_list
species <- read.csv("gene_to_species.csv", sep=";", row.names = 1) # gène_id espèce

# cob gene and species
data_cob <- read.csv("cob_genes.csv", sep=";", row.names = 1) # cob_gene_id cob_gene_name Species

## matrice de comptages ##

#matrice globale
counts <- read.csv("counts.csv", sep=",", row.names=1)
cob_global <- merge(counts, data_cob, by="row.names", all=FALSE) # comptages bruts pur gènes cob

# sous_selection
counts_sub <- merge(counts, species, by = 'row.names', all = FALSE)
counts_sub <- counts_sub[counts_sub$Species == 'Thalassospira sp.',]
rownames(counts_sub) <- counts_sub$Row.names
counts_sub <- counts_sub[,2:10]


#############################################
######### count table analysis ##############
#############################################

### gene on-off

BS <- ifelse(rowMeans(counts[, c("BS2", "BS3", "BS4")]) > 15, 1, 0)
NS <- ifelse(rowMeans(counts[, c("NS2", "NS3", "NS4")]) > 15, 1, 0)
PS <- ifelse(rowMeans(counts[, c("PS2", "PS3", "PS4")]) > 15, 1, 0)

oo_table <- data.frame(BS, NS, PS)

### diagramme de venn

# number of genes active in one limitation
area_BS <- sum(oo_table$BS == 1)
area_NS <- sum(oo_table$NS == 1)
area_PS <- sum(oo_table$PS == 1)

# number of genes active in two limitation
intersect_BS_NS <- sum(oo_table$BS == 1 & oo_table$NS == 1)
intersect_BS_PS <- sum(oo_table$BS == 1 & oo_table$PS == 1)
intersect_NS_PS <- sum(oo_table$NS == 1 & oo_table$PS == 1)

# number of genes active the three limitations
intersect_all <- sum(oo_table$BS == 1 & oo_table$NS == 1 & oo_table$PS == 1)


library(VennDiagram)
library(grid)

venn.plot <- draw.triple.venn(
  area1 = area_BS,
  area2 = area_NS,
  area3 = area_PS,
  n12 = intersect_BS_NS,
  n23 = intersect_NS_PS,
  n13 = intersect_BS_PS,
  n123 = intersect_all,
  category = c("B12", "N", "P"),
  fill = c("red", "green", "blue"),
  cex = 1.5,
  cat.cex = 1.5,
  ind = TRUE
)

grid.draw(venn.plot)

### species visulaisation (propotion)

oo_table <- merge(oo_table, species, by='row.names', all=F)
rownames(oo_table)<-oo_table$Row.names
oo_table$Row.names<- NULL

oo_table$type_activite <- with(oo_table, ifelse(BS == 1 & NS == 0 & PS == 0, "B12",
                                    ifelse(BS == 1 & NS == 1 & PS == 0, "B12&N",
                                           ifelse(BS == 1 & NS == 0 & PS == 1, "B12&P",
                                                  ifelse(BS == 1 & NS == 1 & PS == 1, "B12&N&P",
                                                         ifelse(BS == 0 & NS == 1 & PS == 0, "N",
                                                                ifelse(BS == 0 & NS == 1 & PS == 1, "N&P",
                                                                       ifelse(BS == 0 & NS == 0 & PS == 1, "P",
                                                                              "NA"))))))))

oo_table$GeneID <- rownames(oo_table)
library(dplyr)



prop_onoff <- oo_table %>%
  mutate(
    Category = case_when(
      BS == 1 & NS == 1 & PS == 0 ~ "B12&N",
      BS == 1 & NS == 0 & PS == 1 ~ "B12&P",
      BS == 0 & NS == 1 & PS == 1 ~ "N&P",
      BS == 1 & NS == 0 & PS == 0 ~ "B12",
      BS == 0 & NS == 1 & PS == 0 ~ "N",
      BS == 0 & NS == 0 & PS == 1 ~ "P",
      BS == 1 & NS == 1 & PS == 1 ~ "B12&N&P",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Category)) %>%
  group_by(Species, Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Proportion = Count / sum(Count))  # Normaliser pour que toutes barres soient de même taille

all_genes <- prop_onoff %>%
  group_by(Species) %>%
  summarise(Total_Genes = sum(Count), .groups = "drop")

prop_onoff_annot <- prop_onoff %>%
  left_join(all_genes, by = "Species")

order_species <- c("Marinovum sp.", "Cyclobacterium qasimii", "Marivita cryptomonadis",
                   "Marinobacter sp.", "Muricauda sp.", "Thalassospira sp.", "Alteromonas sp.", "Maribacter sp.")
prop_onoff$Species <- factor(prop_onoff$Species, levels = order_species)
all_genes$Species <- factor(all_genes$Species, levels = order_species)  # pour le geom_text aussi

order_act <- c("B12", "B12&N", "B12&P", "N", "N&P", "P", "B12&N&P")
prop_onoff$Category <- factor(prop_onoff$Category, levels = order_act)



ggplot(prop_onoff, aes(x = Species, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") + 
  geom_text(data = distinct(prop_onoff_annot, Species, Total_Genes), 
            aes(x = Species, y = 1.05, label = paste(Total_Genes)), 
            inherit.aes = FALSE, size = 3.5) +
  scale_fill_manual(values = c("B12" = "#FF6E6E", "N" = "#80FF80", "P" = "#8080FF",
                               "B12&N" = "#8AC66C", "B12&P" = "#A080C0", 
                               "N&P" = "#5DADE2", "B12&N&P" = "#2E4A7D")) +
  theme_minimal() +
  labs(title = "Proportion des gènes actif par espèce",
       x = "Espèce",
       y = "Proportion",
       fill = "Activité")+
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.1))


####################################################
###### deseq2 analysis Thalassopspira sp. ##########
####################################################

dds_t <- DESeqDataSetFromMatrix(countData = counts_sub, colData = colData, design = ~Stress)

#log transformation for visualization
rld <- rlog(dds_t, blind=FALSE)
rld_plot <- plotPCA(rld, intgroup = "Stress")

# change factor level
dds_t$Stress <- relevel(dds_t$Stress, ref = "Nitrogen")

# Run de DESeq2
dds_t <- DESeq(dds_t)
resultsNames(dds_t)

### Azote vs B12 ###
res_tal_n <- results(dds_t, name="Stress_Nitrogen_vs_B12")
res_tal_n <- as.data.frame(res_tal_n)

#### Phosphate vs B12 ###
res_tal_p <- results(dds_t, name="Stress_Phosphorus_vs_B12")
res_tal_p <- as.data.frame(res_tal_p)

### Phosphore vs Nitrigen ####
res_tal_np <- results(dds_t, name="Stress_Phosphorus_vs_Nitrogen")
res_tal_np <- as.data.frame(res_tal_np)

###############################################
##### Plot and analysis of DESeq2 results #####
###############################################
 
# Azote vs B12
vthal_nvb <- ggplot(res_tal_n, aes(x = -log2FoldChange, y = -log10(padj)))+
  geom_point(alpha = 0.6, size = 1) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = 0, color = "darkgrey", size = 1) +
  labs(title = "N-lim vs B12-lim",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

#Phosphate vs B12
vgthal_pvb <- ggplot(res_tal_p, aes(x = -log2FoldChange, y = -log10(padj)))+
  geom_point(alpha = 0.6, size = 1) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = 0, color = "darkgrey", size = 1) +
  labs(title = "P-lim vs B12-lim",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

#Phosphate vs Nitrogen
vthal_nvp <- ggplot(res_tal_np, aes(x = -log2FoldChange, y = -log10(padj)))+
  geom_point(alpha = 0.6, size = 1) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = 0, color = "darkgrey", size = 1) +
  labs(title = "P-lim vs N-lim",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

plot_grid(rld_plot, vthal_nvb, vgthal_pvb, vthal_nvp,  labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)

## cob genes Thalasso ##

cob_thal_n <- subset(res_tal_n, rownames(res_tal_n) %in% rownames(data_cob))
cob_thal_n <- merge(cob_thal_n, data_cob, by='row.names', all=F)
cob_thal_p <- subset(res_tal_p, rownames(res_tal_p) %in% rownames(data_cob))
cob_thal_p <- merge(cob_thal_p, data_cob, by='row.names', all=F)

# Thalassospira sp. cob genes differentially expressed between N-lim and B12-lim

cob_thal_n$log2FoldChange <- -cob_thal_n$log2FoldChange

cobn_plot <- ggplot(subset(cob_thal_n, padj < 0.05), aes(x = log2FoldChange, y = reorder(Gene_name, log2FoldChange), fill = padj)) +
  geom_col() +
  scale_fill_gradient(low = "#D73027", high = "#FEE090", name = "p-adj") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  labs(
    x = "log2 Fold Change",
    y = "Gènes",
    title = "Expression of Vit-B12 biosynthesis genes: N-lim vs B12-lim") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

# Thalassospira sp. cob genes differentially expressed between P-lim and B12-lim

cob_thal_p$log2FoldChange <- -cob_thal_p$log2FoldChange

cobp_plot <- ggplot(subset(cob_thal_p, padj < 0.05), aes(x = log2FoldChange, y = reorder(Gene_name, log2FoldChange), fill = padj)) +
  geom_col() +
  scale_fill_gradient(low = "#D73027", high = "#FEE090", name = "p-adj") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  labs(
    x = "log2 Fold Change",
    y = "Gènes",
    title = "Expression of Vit-B12 biosynthesis genes: P-lim vs B12-lim") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

plot_grid(cobn_plot, cobp_plot, labels=c("A", "B"), ncol = 2, nrow=1)
