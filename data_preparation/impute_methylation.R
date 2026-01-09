# ==============================================================
# GT-Mamba Preprocessing Core (R Script)
# Task: Imputation (methyLImp2) + Normalization (BMIQ)
# ==============================================================

# 1. 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript impute_methylation.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# 2. 加载依赖包
# 注意：BMIQ 需要 wateRmelon，区分探针 Type I/II 需要 450K 注释包
required_pkgs <- c("data.table", "methyLImp2", "BiocParallel", "wateRmelon", "IlluminaHumanMethylation450kanno.ilmn12.hg19")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package", pkg, "is missing. Please install via BiocManager::install()"))
  }
}

cat(sprintf("[R-Core] Reading Data: %s\n", input_file))

# 3. 快速读取数据
beta_matrix <- fread(input_file, data.table = FALSE)
rownames(beta_matrix) <- beta_matrix[,1]
beta_matrix <- beta_matrix[,-1]
data_matrix <- as.matrix(beta_matrix)

# ==============================================================
# Step 1: 缺失值插补 (Imputation)
# Method: methyLImp2 (Array-aware) with Robust Fallback
# ==============================================================
missing_count <- sum(is.na(data_matrix))

if (missing_count > 0) {
  cat(sprintf("[R-Core] Found %d missing values. Running methyLImp2...\n", missing_count))
  
  # 自动检测芯片类型 (EPIC vs 450K)
  array_type <- if(nrow(data_matrix) > 460000) "EPIC" else "450K"
  
  tryCatch({
    data_matrix <- methyLImp2(
      input = data_matrix, 
      type = array_type,
      BPPARAM = MulticoreParam(workers = 4)
    )
    
    # 兜底策略：如果算法未收敛导致仍有 NA，使用列均值填充
    if(sum(is.na(data_matrix)) > 0) {
      k <- which(is.na(data_matrix), arr.ind=TRUE)
      data_matrix[k] <- colMeans(data_matrix, na.rm=TRUE)[k[,2]]
    }
  }, error = function(e) {
    cat("[R-Core] methyLImp2 failed, switching to Mean Imputation.\n")
    # 终极兜底：列均值
    for(i in 1:ncol(data_matrix)) {
      data_matrix[is.na(data_matrix[,i]), i] <- mean(data_matrix[,i], na.rm = TRUE)
    }
  })
} else {
  cat("[R-Core] No missing values. Skipping imputation.\n")
}

# ==============================================================
# Step 2: 数据归一化 (Normalization)
# Method: BMIQ (Beta-mixture Quantile Normalization)
# ==============================================================
cat("[R-Core] Running BMIQ Normalization...\n")

tryCatch({
  # 2.1 获取 450K 注释信息 (区分 Type I 和 Type II 探针)
  ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  common_probes <- intersect(rownames(data_matrix), rownames(ann))
  
  if(length(common_probes) > 1000) {
    # 对齐数据
    data_sub <- data_matrix[common_probes, , drop=FALSE]
    
    # 获取设计向量: Type I = 1, Type II = 2
    design_info <- ann[common_probes, "Type"]
    design_v <- ifelse(design_info == "I", 1, 2)
    
    # 2.2 定义单样本 BMIQ 函数 (含错误捕获)
    run_bmiq <- function(x, des) {
      # tryCatch 防止单个样本报错导致整个程序崩溃
      out <- try(wateRmelon::BMIQ(x, design.v = des, plots = FALSE), silent=TRUE)
      if(inherits(out, "try-error")) return(x) else return(out$nbeta)
    }
    
    # 2.3 对每列(样本)执行 BMIQ
    # 注意：BMIQ 计算较慢，这里使用 apply 逐列处理
    data_norm <- apply(data_sub, 2, run_bmiq, des = design_v)
    
    data_matrix <- data_norm
    cat("[R-Core] BMIQ Normalization completed.\n")
    
  } else {
    cat("[R-Core] Warning: Too few matching probes for annotation. Skipping BMIQ.\n")
  }
  
}, error = function(e) {
  cat("[R-Core] BMIQ Error:", e$message, "\n")
  cat("[R-Core] Skipping BMIQ to ensure pipeline continuity.\n")
})

# ==============================================================
# Step 3: 保存结果
# ==============================================================
cat(sprintf("[R-Core] Saving processed matrix to: %s\n", output_file))
write.csv(data_matrix, file = gzfile(output_file))