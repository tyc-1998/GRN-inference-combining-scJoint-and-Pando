# load library
library(Seurat)
library(Signac)
library(dplyr)
library(Matrix)
library(rhdf5)
library(HDF5Array)

#' This function will generate h5 files for a list expression matrices as input of process_db.py file
write_h5_scJoint <- function(exprs_list, h5file_list) {
  for (i in seq_along(exprs_list)) {
    if (file.exists(h5file_list[i])) {
      warning("h5file exists! will rewrite it.")
      system(paste("rm", h5file_list[i]))
    }
    h5createFile(h5file_list[i])
    h5createGroup(h5file_list[i], "matrix")
    writeHDF5Array(t((exprs_list[[i]])), h5file_list[i], name = "matrix/data")
    h5write(rownames(exprs_list[[i]]), h5file_list[i], name = "matrix/features")
    h5write(colnames(exprs_list[[i]]), h5file_list[i], name = "matrix/barcodes")
    print(h5ls(h5file_list[i]))
  }
}

#' This function will generate csv files for a list of cell types as input of process_db.py file
write_csv_scJoint <- function(cellType_list, csv_list) {
  for (i in seq_along(cellType_list)) {
    if (file.exists(csv_list[i])) {
      warning("csv_list exists! will rewrite it.")
      system(paste("rm", csv_list[i]))
    }
    names(cellType_list[[i]]) <- NULL
    write.csv(cellType_list[[i]], file = csv_list[i])
  }
}

# read data
rna <- readRDS("./demo_data/pbmc_rna.rds")
atac <- readRDS("./demo_data/pbmc_atac.rds")
# save cell name
write.csv(colnames(rna@assays$RNA@counts),"./data/rna_cell.csv")
write.csv(colnames(atac@assays$RNA@counts),"./data/atac_cell.csv")

rna_gene <- rna@assays$RNA@var.features
atac_gene <- rownames(atac@assays$RNA@counts)
gene <- intersect(rna_gene,atac_gene)
# intersect gene num:1880
print(length(gene))          
expres_rna <- Matrix(rna@assays$RNA@counts[gene,],sparse=F)
expres_atac <- Matrix(atac@assays$RNA@counts[gene,],sparse=F)

# generate result file
write_h5_scJoint(exprs_list = list(rna = expres_rna,
                                   atac = expres_atac), 
                 h5file_list = c("./data/exprs_rna.h5", 
                                 "./data/exprs_atac.h5"))
write_csv_scJoint(cellType_list =  list(rna = rna$celltype),
                  csv_list = c("./data/celltype_rna.csv"))