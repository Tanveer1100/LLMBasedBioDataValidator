# Script using affy instead of oligo for E-MTAB-3610 data
# This uses older packages but should work on more installations
# Original link and source code - https://github.com/mdozmorov/E-MTAB-3610/blob/main/Processing_E-MTAB-3610.R

# Install required packages if needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install base packages
install.packages(c("readr", "WGCNA", "writexl", "R.utils", "httr", "jsonlite"))

# Install Bioconductor packages
BiocManager::install(c("affy", "hgu219.db"))

# Load libraries
library(affy)       # Older but widely compatible package
library(readr)
library(hgu219.db)
library(WGCNA)
library(writexl)
library(R.utils)
library(httr)        # For direct API access
library(jsonlite)    # For JSON parsing
library(tools)       # For file handling

# Create a directory for the data
data_dir <- "E-MTAB-3610_data"
if(!dir.exists(data_dir)) dir.create(data_dir)
setwd(data_dir)

# Download the sample annotation file
sdrf_url <- "https://www.ebi.ac.uk/biostudies/files/E-MTAB-3610/E-MTAB-3610.sdrf.txt"
download.file(sdrf_url, "E-MTAB-3610.sdrf.txt")

# Read the sample annotation to get CEL files
sample_annotations <- read_tsv("E-MTAB-3610.sdrf.txt")
cel_files <- sample_annotations$`Array Data File`

# Create a function to download the raw data files
download_cel_files <- function(cel_files) {
  base_url <- "https://www.ebi.ac.uk/biostudies/files/E-MTAB-3610/"
  for(file in cel_files) {
    # Check if file already exists
    if(!file.exists(file)) {
      file_url <- paste0(base_url, file)
      message("Downloading ", file)
      tryCatch({
        download.file(file_url, file)
      }, error = function(e) {
        message("Failed to download ", file, ": ", e$message)
      })
    }
  }
}

# Download CEL files - Note: This might take a while (~5GB)
download_cel_files(cel_files)

# Read the data using affy (consistent with original script)
message("Reading CEL files...")
arrays <- ReadAffy(filenames = list.files(pattern = "\\.cel$", ignore.case = TRUE))

# Normalize using RMA
message("Normalizing data...")
eset <- rma(arrays)

# Get expression matrix
mtx <- exprs(eset)

# Clean up column names
colnames(mtx) <- sub(pattern = ".cel", "", colnames(mtx), fixed = TRUE)
colnames(mtx) <- sub(pattern = ".CEL", "", colnames(mtx), fixed = TRUE)

# Get common array names
common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)

# Subset and match both matrices
mtx <- mtx[, colnames(mtx) %in% common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]
all.equal(colnames(mtx), sample_annotations$`Assay Name`) # Check if matching worked

# Gene annotations
message("Processing gene annotations...")
k <- keys(hgu219.db, keytype="PROBEID")
gene_annotations <- select(hgu219.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="PROBEID")

# Get common gene names
common_genes <- intersect(rownames(mtx), gene_annotations$PROBEID)

# Subset and match both matrices
mtx <- mtx[rownames(mtx) %in% common_genes, ]
gene_annotations <- gene_annotations[gene_annotations$PROBEID %in% common_genes, ]
gene_annotations <- gene_annotations[match(rownames(mtx), gene_annotations$PROBEID), ]
all.equal(rownames(mtx), gene_annotations$PROBEID)  # Check if matching worked

# Aggregate multiple gene IDs (rows) using the MaxMean method
message("Collapsing rows by gene symbol...")
mtx_collapsed <- collapseRows(datET = mtx, rowGroup = gene_annotations$SYMBOL, rowID = rownames(mtx))$datETcollapsed

# Replace array names with cell line names
colnames(mtx_collapsed) <- sample_annotations$`Characteristics[cell line]`

# Utility function for rounding data frame values
round_df <- function(df, digits = 3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# Prepare data to save
message("Preparing data for export...")
# Add Probe IDs to the matrices
mtx <- data.frame(PROBEID = rownames(mtx), mtx)
mtx_collapsed <- data.frame(GENE = rownames(mtx_collapsed), mtx_collapsed)

# Save summarized matrix
write_csv(round_df(mtx_collapsed), "E-MTAB-3610_matrix.csv")
gzip("E-MTAB-3610_matrix.csv")

# Save cell annotations
write_csv(sample_annotations, "E-MTAB-3610_cell_annotations.csv")
gzip("E-MTAB-3610_cell_annotations.csv")

# Save all data in Excel
message("Creating Excel file...")
x <- list(Summarized = mtx_collapsed, Samples = sample_annotations, Original = mtx, Genes = gene_annotations)
write_xlsx(x, "E-MTAB-3610_processed.xlsx")

message("Processing complete!")