# load library
library(FigR)
library(optmatch)
library(data.table)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(doParallel)
library(ggplot2)
library(Pando)
library(grr)
data(motifs)
data(motif2tf)

# read data
rna <- readRDS("./demo_data/pbmc_rna.rds")
atac <- readRDS("./scJoint_result/pbmc_atac_add_label.rds")
DefaultAssay(atac) <- "peaks"

# find peak markers
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac)
Idents(atac) <- atac$celltype
markers <- FindAllMarkers(atac, test.use = 'LR',only.pos = T, logfc.threshold=0.5)
write.csv(markers,'./pando_result/peak_markers.csv')

# pseudomulti
embeddings <- read.csv("./scJoint_result/combine_embedding.csv",row.names=1)

isATAC <- rownames(embeddings) %in% colnames(atac) 
isRNA <- rownames(embeddings) %in% colnames(rna)

ATACcells <- rownames(embeddings)[isATAC]
RNAcells <- rownames(embeddings)[isRNA]

ATAC_PCs <- embeddings[isATAC,]
ATAC_PCs <- as.matrix(ATAC_PCs)
RNA_PCs <- embeddings[isRNA,]
RNA_PCs <- as.matrix(RNA_PCs)

pairing <- pairCells(ATAC = ATAC_PCs,
                     RNA = RNA_PCs,
                     keepUnique = TRUE
                     )

# Read umap information
umap <- read.csv("./scJoint_result/umap_embedding.csv",row.names=1)
meta <- data.frame(cell_name = rownames(umap),
                   celltype = umap$celltype,
                   batch=umap$batch)
ATAC_meta <- meta[meta$batch=="ATAC",c(1,2)]
RNA_meta <- meta[meta$batch=="RNA",c(1,2)]
colnames(ATAC_meta) <- c("cell_name","ATAC_celltype")
colnames(RNA_meta) <- c("cell_name","RNA_celltype")
pairing <- left_join(pairing, ATAC_meta, by = c("ATAC" = "cell_name"))
pairing <- left_join(pairing, RNA_meta, by = c("RNA" = "cell_name"))

pairing <- pairing[pairing$ATAC_celltype == pairing$RNA_celltype,]

## deduplicate
if (length(unique(pairing$RNA)) > length(unique(pairing$ATAC))) {
    dupli_ATAC <- unique(pairing$ATAC[duplicated(pairing$ATAC)])
    dupli_pairing <- pairing[which(pairing$ATAC %in% dupli_ATAC),]
    #dupli_pairing <- data.frame(ATAC = dupli_pairing$ATAC,RNA = dupli_pairing$RNA, dist = dupli_pairing$dist)
    ## deduplicate
    dedupli_pairing <- dupli_pairing %>% group_by(ATAC) %>% top_n(n=-1, wt=dist)
    unique_pairing <- pairing[which(! pairing$ATAC %in% dupli_ATAC),]
    total_uniq_pairing <- rbind(unique_pairing,dedupli_pairing)
}else {
    dupli_RNA <- unique(pairing$RNA[duplicated(pairing$RNA)])
    dupli_pairing <- pairing[which(pairing$RNA %in% dupli_RNA),]
    #dupli_pairing <- data.frame(ATAC = dupli_pairing$ATAC,RNA = dupli_pairing$RNA, dist = dupli_pairing$dist)
    ## deduplicate
    dedupli_pairing <- dupli_pairing %>% group_by(RNA) %>% top_n(n=-1, wt=dist)
    unique_pairing <- pairing[which(! pairing$RNA %in% dupli_RNA),]
    total_uniq_pairing <- rbind(unique_pairing,dedupli_pairing)
}

paired <- data.frame(ATAC = total_uniq_pairing$ATAC,
                     RNA = total_uniq_pairing$RNA,
                     multi = paste0("multi_cell",c(1:length(total_uniq_pairing$ATAC))),
                     multi_celltype=total_uniq_pairing$ATAC_celltype)
#gc()

# Filter cell types with more than 5 cells
paired <- paired %>%
  add_count(multi_celltype) %>%
  filter(n > 5)

print(table(paired$multi_celltype))
multi <- subset(rna, cells = paired$RNA)
multi <- RenameCells(multi,old.names = colnames(multi),new.names = paired$multi)
atac <- subset(atac, cells = paired$ATAC)
atac <- RenameCells(atac,old.names = colnames(atac),new.names = paired$multi)
multi[["peaks"]] <- atac[["peaks"]]
# multi$multi_celltype <- paired$multi_celltype
multi$multi_celltype <- paired[which(paired$multi %in% rownames(multi@meta.data)), "multi_celltype"]
saveRDS(multi, "./pando_result/multi.rds")

gc()

# difference peak
chr <- unlist(lapply(strsplit(as.character(markers$gene),"-"),"[",1))
start <- unlist(lapply(strsplit(as.character(markers$gene),"-"),"[",2))
end <- unlist(lapply(strsplit(as.character(markers$gene),"-"),"[",3))
start <- as.numeric(start)
end <- as.numeric(end)
gr <- GRanges(seqnames = Rle(chr),ranges = IRanges(start, end = end))


# Pando
seurat_object <- initiate_grn(multi,
                              rna_assay = 'RNA',
                              peak_assay = 'peaks',
                              regions = gr)

DefaultAssay(seurat_object) <- "peaks"
seurat_object <- find_motifs(
    seurat_object,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg19,
    motif_tfs = motif2tf
)

#saveRDS(seurat_object,"./pando_result/beforeGRN.rds")

# infer
# use highly variable genes for inference, which can also be customized
genesUse <- seurat_object@assays$RNA@var.features
registerDoParallel(8)

# inter_peak <- intersect(rownames(seurat_object@grn@regions@motifs@data),
#                         rownames(seurat_object@assays$peaks))

# seurat_object@grn@regions@motifs@data <- seurat_object@grn@regions@motifs@data[inter_peak,]
# seurat_object@grn@regions@peaks <- match(inter_peak, rownames(seurat_object@assays$peaks))

seurat_object <- infer_grn(
    seurat_object,
    peak_to_gene_method = 'GREAT',
    genes = genesUse,
    parallel = T
)

seurat_object <- find_modules(
 seurat_object, 
 p_thresh = 0.1,
 nvar_thresh = 2, 
 min_genes_per_module = 1, 
 rsq_thresh = 0.05
)

modules <- NetworkModules(seurat_object) 
head(modules@meta)
saveRDS(seurat_object,"./pando_result/pbmc_GRN.rds")

p1 <- plot_gof(seurat_object, point_size=3)
p2 <- plot_module_metrics(seurat_object)
ggsave("./pando_result/goodness-of-fit_metrics_plot.pdf",p1)
ggsave("./pando_result/module_metrics_plot.pdf",p2)
seurat_object <- get_network_graph(seurat_object)
p3 <- plot_network_graph(seurat_object, text_size=6)
ggsave("./pando_result/GRN_plot.pdf",p3)
Idents(seurat_object) <- "celltype"
p4 <- DimPlot(seurat_object, reduction = "umap",label = T,label.size = 5)
ggsave("./pando_result/umap_plot.pdf",p4)