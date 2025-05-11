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

BiocManager::install(c("affy", "hgu219.db", "readr", "writexl"), ask = FALSE, update = FALSE)

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
download.file("ftp://ftp.
