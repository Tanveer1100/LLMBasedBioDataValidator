# ---- Step 1: Install necessary packages (skip if already installed) ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ArrayExpress", "affy", "hgu219.db", "WGCNA", "readr", "writexl"), ask = FALSE, update = FALSE)

# ---- Step 2: Load packages ----
library(ArrayExpress)
library(affy)
library(readr)
library(hgu219.db)
library(WGCNA)
library(writexl)

# ---- Step 3: Download and unzip CEL files ----
getAE(accession = "E-MTAB-3610")
for (f in list.files(pattern = "zip")) {
  unzip(f, junkpaths = TRUE)
}

# ---- Step 4: Normalize and get expression matrix ----
arrays <- read.affybatch(filenames = list.files(pattern = "cel"))
eset <- rma(arrays)
mtx <- exprs(eset)
colnames(mtx) <- sub(pattern = ".cel", "", colnames(mtx), fixed = TRUE)

# ---- Step 5: Load sample annotations ----
sample_annotations <- read_tsv("E-MTAB-3610.sdrf.txt")
common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)
mtx <- mtx[, colnames(mtx) %in% common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]
stopifnot(all.equal(colnames(mtx), sample_annotations$`Assay Name`))

# ---- Step 6: Annotate probes with gene names ----
k <- keys(hgu219.db, keytype = "PROBEID")
gene_annotations <- select(hgu219.db, keys = k, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
common_genes <- intersect(rownames(mtx), gene_annotations$PROBEID)
mtx <- mtx[rownames(mtx) %in% common_genes, ]
gene_annotations <- gene_annotations[gene_annotations$PROBEID %in% common_genes, ]
gene_annotations <- gene_annotations[match(rownames(mtx), gene_annotations$PROBEID), ]
stopifnot(all.equal(rownames(mtx), gene_annotations$PROBEID))

# ---- Step 7: Collapse probes to gene level ----
collapsed <- collapseRows(datET = mtx, rowGroup = gene_annotations$SYMBOL, rowID = rownames(mtx))
mtx_collapsed <- collapsed$datETcollapsed
colnames(mtx_collapsed) <- sample_annotations$`Characteristics[cell line]`

# ---- Step 8: Save data locally ----
mtx <- data.frame(PROBEID = rownames(mtx), mtx)
mtx_collapsed <- data.frame(GENE = rownames(mtx_collapsed), mtx_collapsed)

write_csv(mtx_collapsed, "E-MTAB-3610_matrix.csv")
write_csv(sample_annotations, "E-MTAB-3610_cell_annotations.csv")

# Write Excel
x <- list(Summarized = mtx_collapsed, Samples = sample_annotations, Original = mtx, Genes = gene_annotations)
write_xlsx(x, "E-MTAB-3610_processed.xlsx")

# ---- Step 9: Upload files to S3 ----
# Set your bucket and folder path
bucket_name <- "your-s3-bucket-name"
folder_path <- "gdsc/E-MTAB-3610/"

system(paste0("aws s3 cp E-MTAB-3610_matrix.csv s3://", bucket_name, "/", folder_path))
system(paste0("aws s3 cp E-MTAB-3610_cell_annotations.csv s3://", bucket_name, "/", folder_path))
system(paste0("aws s3 cp E-MTAB-3610_processed.xlsx s3://", bucket_name, "/", folder_path))
