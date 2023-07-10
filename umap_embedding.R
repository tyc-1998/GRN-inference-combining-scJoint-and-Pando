# Package and functions
library(Signac)
library(Seurat)
library(ggplot2)
library(scattermore)
library(uwot)
library(grDevices)
library(dplyr)
library(tibble)

# Reading label reference
label_class <- read.delim("./data/label_to_idx.txt",header = FALSE)
label_class_num <- unlist(lapply(strsplit(label_class$V1, " "),
                                 function(x) x[length(x)]))
label_class_name <- unlist(lapply(strsplit(label_class$V1, " "),
                                  function(x) paste(x[-length(x)], collapse = " ")))
label_class <- data.frame(name = label_class_name,num = label_class_num)

# The folder with output
results_dir <- './output'
embedding_files <- list.files(results_dir, "embeddings.txt")

embedding <- list()
for (i in 1:length(embedding_files)) {
    embedding[[i]] <- read.delim(file.path(results_dir, embedding_files[i]),
                                 header = FALSE, sep = " ")
}

names(embedding) <- gsub("_embeddings.txt", "", embedding_files)

cat("Dimension of embedding: ")
print(lapply(embedding, dim))

# Reading KNN prediction
knn_prediction_files <- list.files(results_dir, pattern = "knn_predictions.txt")

knn_prediction <- list()
for (i in 1:length(knn_prediction_files)) {
    knn_prediction[[i]] <- read.delim(file.path(results_dir, knn_prediction_files[i]),
                                      header = FALSE, sep = " ")
    knn_prediction[[i]] <- label_class$name[knn_prediction[[i]]$V1 + 1]
}

names(knn_prediction) <- gsub("_knn_predictions.txt", "", knn_prediction_files)

rna_dataset <- setdiff(names(embedding), names(knn_prediction))
print(rna_dataset)
rna_prediction <- list()
for (i in 1:length(rna_dataset)) {
    rna_prediction[[i]] <- read.delim(file.path(results_dir, paste0(rna_dataset[i], "_predictions.txt")),
                                      header = FALSE, sep = " ")
    rna_prediction[[i]] <- label_class$name[apply(rna_prediction[[i]], 1, which.max)]
}

names(rna_prediction) <- rna_dataset

prediction_list <- append(rna_prediction, knn_prediction)
prediction_list <- prediction_list[names(embedding)]

batch <- rep(names(prediction_list), unlist(lapply(prediction_list, length)))
combine_embedding <- do.call(rbind, embedding)
prediction <- do.call(c, prediction_list)

idx <- sort(sample(length(batch), round(length(batch))))
combine_embedding <- combine_embedding[idx, ]
prediction <- prediction[idx]
batch <- batch[idx]

cat("Dimension to be visualised: ")
print(dim(combine_embedding))

set.seed(1234)
umap_res <- uwot::umap(combine_embedding,min_dist = 0.1)

df <- data.frame(UMAP1 = umap_res[, 1], UMAP2 = umap_res[, 2],
                 celltype = prediction,
                 batch = batch)

rna <- read.csv("./data/rna_cell.csv")
atac <- read.csv("./data/atac_cell.csv")
cell <- c(atac$x,rna$x)
rownames(combine_embedding) <- cell
write.csv(combine_embedding,"./scJoint_result/combine_embedding.csv")
rownames(df) <- cell
df[grep("atac$",df$batch),]$batch <- "ATAC"
df[grep("rna$",df$batch),]$batch <- "RNA"
write.csv(df, "./scJoint_result/umap_embedding.csv")

atac<- readRDS("./demo_data/pbmc_atac.rds")
atac$celltype <- df[rownames(atac@meta.data),]$celltype
saveRDS(atac,"./scJoint_result/pbmc_atac_add_label.rds")

colors = c(
  "#7FC97F", "#BEAED4", "#FDC086", "#386CB0", "#F0027F", "#E78AC3",
  "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#B3DE69",
  "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928",
  "#FBB4AE", "#B3CDE3", "#BC80BD", "#CCEBC5", "#DECBE4", "#FED9A6",
  "#E5D8BD", "#FDDAEC", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4",
  "#E6F5C9", "#FFF2AE", "#F1E2CC", "#E41A1C", "#377EB8", "#984EA3",
  "#A65628", "#F781BF", "#FFED6F", "#66C2A5", "#FC8D62", "#8DA0CB", 
  "#A6D854", "#FFD92F", "#E5C494", "#8DD3C7", "#FFFFB3", "#FCCDE5",
  "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")
names(colors) <- unique(df$celltype)
colors <- colors[1:length(unique(df$celltype))]

# draw RNA, ATAC integrated cell type UMAP
df$batch <- factor(df$batch,levels=c('RNA','ATAC'))
median_df <- df %>% group_by(celltype) %>% summarise(x = median(UMAP1), y = median(UMAP2))
p=ggplot(df)+geom_point(aes(UMAP1,UMAP2,color=celltype),size=0.2,alpha=0.8)+
  cowplot::theme_cowplot() + facet_wrap(batch~.)+ 
  ggrepel::geom_text_repel(data = median_df, aes(x,y, label = celltype), min.segment.length = 0,size=4) + 
  scale_color_manual(values=colors) + theme(strip.background = element_blank()) + 
  labs(title='', x='UMAP_1', y='UMAP_2') + theme(legend.position="right")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=14),
  strip.text = element_text(size = 20))+ guides(color=guide_legend(override.aes = list(size=5,alpha=1)))+ 
  theme(legend.title = element_text(size = 12,face="bold"), legend.text = element_text(size = 10),legend.key.size = unit(0.1, "inches"))
ggsave(p,file='./scJoint_result/RNA_ATAC.png',width=10,height=5,units='in')