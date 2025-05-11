# -------- STEP 0: Setup Working Directory --------
dir.create("/home/ec2-user/SageMaker/E-MTAB-3610", recursive = TRUE, showWarnings = FALSE)
setwd("/home/ec2-user/SageMaker/E-MTAB-3610")

# -------- STEP 1: Download all CEL zip files --------
for (i in 1:25) {
  url <- sprintf("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3610/E-MTAB-3610.raw.%d.zip", i)
  dest <- sprintf("E-MTAB-3610.raw.%d.zip", i)
  download.file(url, destfile = dest, mode = "wb")
}

# -------- STEP 2: Unzip all CEL files --------
for (f in list.files(pattern = "zip")) {
  unzip(f, junkpaths = TRUE)
}

# -------- STEP 3: Install Required Packages --------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("affy", "hgu219.db", "readr", "writexl"), ask = FALSE, update = FALSE, force= TRUE )

library(affy)
library(hgu219.db)
library(readr)
library(writexl)

# -------- STEP 4: Limit threads to avoid pthread errors --------
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

# -------- STEP 5: Normalize CEL files in safe batches --------
cel_files <- list.files(pattern = "cel", ignore.case = TRUE)
batch_size <- 300
all_exprs <- list()

for (i in seq(1, length(cel_files), by = batch_size)) {
  message(sprintf("Processing CEL files %d to %d", i, min(i + batch_size - 1, length(cel_files))))
  batch_files <- cel_files[i:min(i + batch_size - 1, length(cel_files))]
  arrays <- read.affybatch(filenames = batch_files)
  eset <- rma(arrays)
  all_exprs[[length(all_exprs) + 1]] <- exprs(eset)
}

# Merge all batches
mtx <- do.call(cbind, all_exprs)
colnames(mtx) <- sub(".cel", "", colnames(mtx), fixed = TRUE)

# -------- STEP 6: Download and Match Sample Annotations --------
download.file("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3610/E-MTAB-3610.sdrf.txt", "E-MTAB-3610.sdrf.txt")
sample_annotations <- read_tsv("E-MTAB-3610.sdrf.txt")

common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)
mtx <- mtx[, common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]
stopifnot(all.equal(colnames(mtx), sample_annotations$`Assay Name`))

# -------- STEP 7: Save Results --------
write_csv(data.frame(PROBEID = rownames(mtx), mtx), "E-MTAB-3610_matrix_probe.csv")
write_csv(sample_annotations, "E-MTAB-3610_cell_annotations.csv")

# Optional: save Excel
write_xlsx(list(ProbeMatrix = mtx, Samples = sample_annotations), "E-MTAB-3610_probe_data.xlsx")

# -------- STEP 8: Upload to S3 --------
bucket <- "your-s3-bucket-name"  # REPLACE THIS
folder <- "gdsc/E-MTAB-3610/"

system(paste0("aws s3 cp E-MTAB-3610_matrix_probe.csv s3://", bucket, "/", folder))
system(paste0("aws s3 cp E-MTAB-3610_cell_annotations.csv s3://", bucket, "/", folder))
system(paste0("aws s3 cp E-MTAB-3610_probe_data.xlsx s3://", bucket, "/", folder))
