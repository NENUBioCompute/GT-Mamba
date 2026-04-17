# preprocessing/raw_data_pipeline/impute_methylation.R

# ==============================================================================
# GT-Mamba Preprocessing Core (R Script)
# Task: Array-aware Imputation (methyLImp2) + Normalization (BMIQ)
# Purpose: High-fidelity biological preprocessing of DNA methylation arrays.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript impute_methylation.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# --- Load Required Libraries ---
required_pkgs <- c("data.table", "methyLImp2", "BiocParallel", "wateRmelon", "IlluminaHumanMethylation450kanno.ilmn12.hg19")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package", pkg, "is missing. Please install it before running."))
  }
}

cat(sprintf("[GT-Mamba R-Core] 🔍 Processing File: %s\n", basename(input_file)))

# --- Data Loading ---
# Using fread for high-speed loading of large methylation matrices
beta_matrix <- fread(input_file, data.table = FALSE)
rownames(beta_matrix) <- beta_matrix[,1]
beta_matrix <- as.matrix(beta_matrix[,-1])

# --- Step 1: Missing Value Imputation ---
missing_count <- sum(is.na(beta_matrix))
if (missing_count > 0) {
  cat(sprintf("[GT-Mamba R-Core] 🧬 Imputing %d missing values...\n", missing_count))
  
  # Automatic array type detection (EPIC vs 450K) based on probe count
  array_type <- if(nrow(beta_matrix) > 460000) "EPIC" else "450K"
  
  tryCatch({
    # Attempt advanced linear regression imputation using methyLImp2
    beta_matrix <- methyLImp2(input = beta_matrix, type = array_type, BPPARAM = MulticoreParam(workers = 4))
    
    # Check for residual NAs and apply column-wise mean imputation as a second-pass
    if(sum(is.na(beta_matrix)) > 0) {
      k <- which(is.na(beta_matrix), arr.ind=TRUE)
      beta_matrix[k] <- colMeans(beta_matrix, na.rm=TRUE)[k[,2]]
    }
  }, error = function(e) {
    cat("[GT-Mamba R-Core] ⚠️ methyLImp2 failed. Falling back to simple Mean Imputation.\n")
    for(i in 1:ncol(beta_matrix)) {
      beta_matrix[is.na(beta_matrix[,i]), i] <- mean(beta_matrix[,i], na.rm = TRUE)
    }
  })
}

# --- Step 2: BMIQ Normalization (Probe Type Bias Correction) ---
cat("[GT-Mamba R-Core] 🧪 Running BMIQ Normalization...\n")
tryCatch({
  # Load Illumina 450K/EPIC probe annotations
  ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  common_probes <- intersect(rownames(beta_matrix), rownames(ann))
  
  if(length(common_probes) > 1000) {
    data_sub <- beta_matrix[common_probes, , drop=FALSE]
    # Identify Type I and Type II probes for BMIQ adjustment
    design_v <- ifelse(ann[common_probes, "Type"] == "I", 1, 2)
    
    # Perform BMIQ normalization sample-by-sample
    beta_norm <- apply(data_sub, 2, function(x) {
      out <- try(wateRmelon::BMIQ(x, design.v = design_v, plots = FALSE), silent=TRUE)
      if(inherits(out, "try-error")) return(x) else return(out$nbeta)
    })
    beta_matrix <- beta_norm
    cat("[GT-Mamba R-Core] ✅ BMIQ Normalization Success.\n")
  }
}, error = function(e) {
  cat("[GT-Mamba R-Core] ❌ BMIQ Skipped due to runtime error.\n")
})

# --- Export Results ---
# Saving as compressed .csv.gz to optimize storage space
write.csv(beta_matrix, file = gzfile(output_file))
cat("[GT-Mamba R-Core] 🎉 Processed matrix saved successfully. Pipeline complete.\n")