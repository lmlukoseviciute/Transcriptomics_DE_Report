#!/usr/bin/env Rscript
### Loading libraries
suppressMessages(library(DESeq2))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(VennDiagram))


# Loading count matrix
cts_colli <- read.csv("inputs/counts_colli.csv", sep="," , row.names="gene_id")
cts_kapa <- read.csv("inputs/counts_kapa.csv", sep="," , row.names="gene_id")


# Loading columns names info
coldata <- read.csv("inputs/colnames.csv", row.names=1)
coldata$condition <- factor(coldata$condition)
coldata <- coldata[,c("condition","type")]
coldata$type <- factor(coldata$type)

################################################
############### DESeq2 #########################
################################################

dds_colli <- DESeqDataSetFromMatrix(countData = cts_colli, colData = coldata, design = ~ condition)
dds_colli <- DESeq(dds_colli)

dds_kapa <- DESeqDataSetFromMatrix(countData = cts_kapa, colData = coldata, design = ~ condition)
dds_kapa <- DESeq(dds_kapa)

# getting expression results
rezults_colli = results(dds_colli, contrast = c('condition', 'normal', 'cancerous'))
rezults_kapa = results(dds_kapa, contrast = c('condition', 'normal', 'cancerous'))


######################################################################
################################################
################# Vulcano ######################
################################################

# Volcano plot
png( "./outputs/plots/volcano_colli.png",
     width = 20, height = 25,
     units = "cm",
     res = 600,
     pointsize = 2 )
EnhancedVolcano(rezults_colli, title = 'Vulcano plot for Collibri', lab = rownames(rezults_colli), 
                x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.00000005)

png( "./outputs/plots/volcano_kapa.png",
     width = 20, height = 25,
     units = "cm",
     res = 600,
     pointsize = 2 )
EnhancedVolcano(rezults_kapa, title = 'Vulcano plot for KAPA', lab = rownames(rezults_kapa),
                x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.00000005)
dev = dev.off()

################################################
################# Lists ########################
################################################

# list of DE genes
gene_id = rezults_colli@rownames
p_adj = rezults_colli$padj

list_colli = cbind(gene_id, p_adj)
write.csv(list_colli, "outputs/DE_genes/Collibri.csv")

gene_id = rezults_kapa@rownames
p_adj = rezults_kapa$padj

list_kapa = cbind(gene_id, p_adj)
write.csv(list_kapa, "outputs/DE_genes/KAPA.csv")


################################################
############## Venn diagram ####################
################################################
colli = subset(rezults_colli, pvalue<0.00000005)
kapa = subset(rezults_kapa, pvalue<0.00000005)

write.csv(colli, "outputs/DE_genes/Collibri_pval.csv")
write.csv(kapa, "outputs/DE_genes/KAPA_pval.csv")

colli = row.names(colli[1])
kapa = row.names(kapa[1])


png( "./outputs/plots/venn_diag.png",
     width = 20, height = 25,
     units = "cm",
     res = 600,
     pointsize = 2 )

v_diag <- venn.diagram(list(colli, kapa), NULL, fill=c("blue", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Collibri", "KAPA"))
grid.draw(v_diag)
dev.off()

