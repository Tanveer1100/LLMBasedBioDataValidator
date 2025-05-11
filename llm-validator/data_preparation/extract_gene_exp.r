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
BiocManager::install(c("affy", "readr", "writexl", "hgu219.db"), ask = FALSE, update = FALSE)

library(affy)
library(readr)
library(writexl)
library(hgu219.db)

# -------- STEP 4: Limit threads to avoid errors --------
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# -------- STEP 5: Normalize CEL files in safe batches --------
cel_files <- list.files(pattern = "cel", ignore.case = TRUE)
batch_size <- 100
all_batches <- list()

for (i in seq(1, length(cel_files), by = batch_size)) {
  message(sprintf("🟢 Processing CEL files %d to %d", i, min(i + batch_size - 1, length(cel_files))))
  
  batch_files <- cel_files[i:min(i + batch_size - 1, length(cel_files))]
  
  tryCatch({
    arrays <- read.affybatch(filenames = batch_files)
    eset <- rma(arrays)
    expr_matrix <- exprs(eset)
    all_batches[[length(all_batches) + 1]] <- expr_matrix
    
    # Save intermediate batch to disk
    write_csv(data.frame(PROBEID = rownames(expr_matrix), expr_matrix),
              sprintf("E-MTAB-3610_batch_%03d.csv", i))
    
  }, error = function(e) {
    message(sprintf("❌ Batch %d–%d failed: %s", i, i + batch_size - 1, e$message))
  })
}

# -------- STEP 6: Merge all batches into full matrix --------
mtx <- do.call(cbind, all_batches)
colnames(mtx) <- sub(".cel", "", colnames(mtx), fixed = TRUE)

# -------- STEP 7: Download and apply sample annotations --------
download.file("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3610/E-MTAB-3610.sdrf.txt", "E-MTAB-3610.sdrf.txt")
sample_annotations <- read_tsv("E-MTAB-3610.sdrf.txt")

common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)
mtx <- mtx[, common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]
stopifnot(all.equal(colnames(mtx), sample_annotations$`Assay Name`))

# -------- STEP 8: Save Final Outputs --------
write_csv(data.frame(PROBEID = rownames(mtx), mtx), "E-MTAB-3610_matrix_probe.csv")
write_csv(sample_annotations, "E-MTAB-3610_cell_annotations.csv")
write_xlsx(list(ProbeMatrix = mtx, Samples = sample_annotations), "E-MTAB-3610_probe_data.xlsx")

# -------- STEP 9: Upload to S3 --------
bucket <- "your-s3-bucket-name"  # 🔁 REPLACE this with your bucket
folder <- "gdsc/E-MTAB-3610/"

system(paste0("aws s3 cp E-MTAB-3610_matrix_probe.csv s3://", bucket, "/", folder))
system(paste0("aws s3 cp E-MTAB-3610_cell_annotations.csv s3://", bucket, "/", folder))
system(paste0("aws s3 cp E-MTAB-3610_probe_data.xlsx s3://", bucket, "/", folder))
