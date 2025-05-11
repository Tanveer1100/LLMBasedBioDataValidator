# Set working directory to the SageMaker notebook instance's home directory
setwd("/home/ec2-user/SageMaker")

# Install required packages if not already installed
install.packages(c("ArrayExpress", "affy", "readr", "WGCNA", "writexl", "aws.s3"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("hgu219.db", update = FALSE)

# Load required libraries
library(ArrayExpress)
library(affy)
library(readr)
library(WGCNA)
library(writexl)
library(aws.s3)
library(hgu219.db)

# Download the data, ~5Gb
getAE(accession = "E-MTAB-3610")
# Unzip files
for (f in list.files(pattern = "zip")) {
  unzip(f, junkpaths = TRUE)
}

# Read and process the data
arrays <- read.affybatch(filenames = list.files(pattern = "cel"))
eset = rma(arrays)
mtx <- exprs(eset)
colnames(mtx) <- sub(pattern = ".cel", "", colnames(mtx), fixed = TRUE)

# Cell annotations
sample_annotations <- read_tsv("E-MTAB-3610.sdrf.txt")
common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)
mtx <- mtx[, colnames(mtx) %in% common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]

# Gene annotations
k <- keys(hgu219.db,keytype="PROBEID")
gene_annotations <- select(hgu219.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="PROBEID")
common_genes <- intersect(rownames(mtx), gene_annotations$PROBEID)
mtx <- mtx[rownames(mtx) %in% common_genes, ]
gene_annotations <- gene_annotations[gene_annotations$PROBEID %in% common_genes, ]
gene_annotations <- gene_annotations[match(rownames(mtx), gene_annotations$PROBEID), ]

# Aggregate gene IDs
mtx_collapsed <- collapseRows(datET = mtx, rowGroup = gene_annotations$SYMBOL, rowID = rownames(mtx))$datETcollapsed
colnames(mtx_collapsed) <- sample_annotations$`Characteristics[cell line]`

# Prepare data
mtx <- data.frame(PROBEID = rownames(mtx), mtx)
mtx_collapsed <- data.frame(GENE = rownames(mtx_collapsed), mtx_collapsed)

# Define S3 location
s3_bucket <- "your-bucket-name"
s3_prefix <- "your/path/E-MTAB-3610/"

# Save summarized matrix to S3
temp_matrix <- tempfile(fileext = ".csv")
write_csv(round_df(mtx_collapsed), temp_matrix)
aws.s3::put_object(
  file = temp_matrix,
  object = paste0(s3_prefix, "E-MTAB-3610_matrix.csv"),
  bucket = s3_bucket
)
unlink(temp_matrix)

# Save cell annotations to S3
temp_annotations <- tempfile(fileext = ".csv")
write_csv(sample_annotations, temp_annotations)
aws.s3::put_object(
  file = temp_annotations,
  object = paste0(s3_prefix, "E-MTAB-3610_cell_annotations.csv"),
  bucket = s3_bucket
)
unlink(temp_annotations)

# Save all data in Excel to S3
x <- list(Summarized = mtx_collapsed, 
          Samples = sample_annotations, 
          Original = mtx,  
          Genes = gene_annotations)
temp_xlsx <- tempfile(fileext = ".xlsx")
write_xlsx(x, temp_xlsx)
aws.s3::put_object(
  file = temp_xlsx,
  object = paste0(s3_prefix, "E-MTAB-3610_processed.xlsx"),
  bucket = s3_bucket
)
unlink(temp_xlsx)
