library(data.table)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(mltools)
options(scipen = 999)

# used genes
dt <- fread('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/panhandles_preprocessed_filtered.tsv')
dt$gene_id <- str_split_fixed(dt$gene_id, '\\.', 2)[, 1]
dt$gene_length <- dt$end_gene - dt$start_gene + 1
dt <- dt[, c('gene_id', 'gene_length'), with = F]
dt <- dt[!duplicated(dt), ]

gene <- dt$gene_id
gene.df <- bitr(gene, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
gene.df <- merge(gene.df, dt, by.x = 'ENSEMBL', by.y = 'gene_id', all.x = T)
gene.df <- gene.df[order(gene.df$gene_length), ]

# match genes w and w/o ph by length
anno <- fread('../conservative_features/gencode.v19.annotation.gtf', skip = 5)
anno <- subset(anno, V3 == 'gene')
anno$gene_id <- str_split_fixed(anno$V9, '; ', 2)[, 1]
anno$gene_id <- gsub('gene_id ', '', anno$gene_id)
anno$gene_id <- gsub('\\"', '', anno$gene_id)
anno$gene_id <- str_split_fixed(anno$gene_id, '\\.', 2)[, 1]
anno$gene_length <- anno$V5 - anno$V4
anno <- anno[, c('gene_id', 'gene_length'), with = F]
anno <- anno[order(anno$gene_length), ]
anno$ph <- F
anno[anno$gene_id %in% gene.df$ENSEMBL, ]$ph <- T
convert <- bitr(anno$gene_id, fromType = "ENSEMBL",
                         toType = c("ENTREZID", "SYMBOL"),
                         OrgDb = org.Hs.eg.db)
anno2 <- merge(anno, convert, by.x = 'gene_id', by.y = 'ENSEMBL')
anno2 <- anno2[order(anno2$gene_length), ]
anno2$bin <- bin_data(x = anno2$gene_length, bins = 300)

x <- table(anno2$bin, anno2$ph)
sum(x[, 'FALSE'] < x[, 'TRUE'])
x[x[, 'FALSE'] < x[, 'TRUE'], ]
gene.list <- lapply(rownames(x), function(bin.){
  take <- min(x[bin., ])
  tmp <- anno2[anno2$bin == bin., ]
  tmp_ph <- subset(tmp, ph)
  tmp_no_ph <- subset(tmp, !ph)
  take_ph <- sample(nrow(tmp_ph), size = take)
  tmp_ph <- tmp_ph[take_ph, ]
  take_no_ph <- sample(nrow(tmp_no_ph), size = take)
  tmp_no_ph <- tmp_no_ph[take_no_ph, ]
  tmp <- rbind(tmp_ph, tmp_no_ph)
  tmp
})
gene.population <- Reduce(rbind, gene.list)
dim(gene.population)/2

# GO
ego2 <- enrichGO(gene = subset(gene.population, ph)$gene_id,
                 universe = gene.population$gene_id,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
ego2@result$p.adjust <- as.numeric(format(ego2@result$p.adjust, digits = 2))
options(scipen = 1)
pdf('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/GO/GO_enrichment.pdf')
barplot(ego2, showCategory=20)
dev.off()

# Disease
edo <- enrichDO(gene = subset(gene.population, ph)$ENTREZID,
                universe = gene.population$ENTREZID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
edo@result$p.adjust <- as.numeric(format(edo@result$p.adjust, digits = 2))
pdf('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/GO/DO_enrichment.pdf')
barplot(edo, showCategory=20)
dev.off()

# cancer genes
ncg <- enrichNCG(subset(gene.population, ph)$ENTREZID, universe = gene.population$ENTREZID)
pdf('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/cancer_enrichment.pdf')
barplot(ncg, showCategory=20)
dev.off()

# Reactome
kegg <- enrichKEGG(gene = subset(gene.population, ph)$ENTREZID, 
                   universe = gene.population$ENTREZID, 
                   organism = 'hsa', pvalueCutoff = 0.05)
pdf('../python_scripts/folding_pretty_copy/out/hg19_ss_flanks/KEGG_enrichment.pdf')
barplot(kegg, showCategory=20)
dev.off()