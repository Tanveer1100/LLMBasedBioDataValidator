# --------- STEP 0: Setup working directory ---------
dir.create("/home/ec2-user/SageMaker/E-MTAB-3610", recursive = TRUE, showWarnings = FALSE)
setwd("/home/ec2-user/SageMaker/E-MTAB-3610")

# --------- STEP 1: Download all CEL zip files manually ---------
for (i in 1:25) {
  url <- sprintf("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3610/E-MTAB-3610.raw.%d.zip", i)
  dest <- sprintf("E-MTAB-3610.raw.%d.zip", i)
  download.file(url, destfile = dest, mode = "wb")
}

# --------- STEP 2: Unzip all zip files ---------
for (f in list.files(pattern = "zip")) {
  unzip(f, junkpaths = TRUE)
}

# --------- STEP 3: Install and load required R packages ---------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("affy", "hgu219.db", "WGCNA", "readr", "writexl"), ask = FALSE, update = FALSE)

library(affy)
library(hgu219.db)
library(WGCNA)
library(readr)
library(writexl)

# --------- STEP 4: Normalize CEL files ---------
arrays <- read.affybatch(filenames = list.files(pattern = "cel", ignore.case = TRUE))
eset <- rma(arrays)
mtx <- exprs(eset)
colnames(mtx) <- sub(".cel", "", colnames(mtx), fixed = TRUE)

# --------- STEP 5: Download and load sample annotations ---------
download.file("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3610/E-MTAB-3610.sdrf.txt", "E-MTAB-3610.sdrf.txt")
sample_annotations <- read_tsv("E-MTAB-3610.sdrf.txt")

common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)
mtx <- mtx[, common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]
stopifnot(all.equal(colnames(mtx), sample_annotations$`Assay Name`))

# --------- STEP 6: Map probes to gene symbols ---------
k <- keys(hgu219.db, keytype = "PROBEID")
gene_annotations <- select(hgu219.db, keys = k, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")

common_genes <- intersect(rownames(mtx), gene_annotations$PROBEID)
mtx <- mtx[common_genes, ]
gene_annotations <- gene_annotations[match(common_genes, gene_annotations$PROBEID), ]
stopifnot(all.equal(rownames(mtx), gene_annotations$PROBEID))

# --------- STEP 7: Collapse probes to gene-level expression ---------
collapsed <- collapseRows(datET = mtx, rowGroup = gene_annotations$SYMBOL, rowID = rownames(mtx))
mtx_collapsed <- collapsed$datETcollapsed
colnames(mtx_collapsed) <- sample_annotations$`Characteristics[cell line]`

# --------- STEP 8: Save output locally ---------
write_csv(data.frame(GENE = rownames(mtx_collapsed), mtx_collapsed), "E-MTAB-3610_matrix.csv")
write_csv(sample_annotations, "E-MTAB-3610_cell_annotations.csv")
write_xlsx(list(Summarized = mtx_collapsed, Samples = sample_annotations, Original = mtx, Genes = gene_annotations),
           "E-MTAB-3610_processed.xlsx")

# --------- STEP 9: Upload to S3 (set your bucket name) ---------
bucket <- "your-s3-bucket-name"   # <== REPLACE THIS
folder <- "gdsc/E-MTAB-3610/"

system(paste0("aws s3 cp E-MTAB-3610_matrix.csv s3://", bucket, "/", folder))
system(paste0("aws s3 cp E-MTAB-3610_cell_annotations.csv s3://", bucket, "/", folder))
system(paste0("aws s3 cp E-MTAB-3610_processed.xlsx s3://", bucket, "/", folder))
