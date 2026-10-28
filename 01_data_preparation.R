# 01 Data preparation (sources: 拆分脚本.R / 分析数据.R / 分析数据个数.R)


# =====================================================
# Source: 拆分脚本.R
# =====================================================

# =====================================================
# 完整脚本：使用230bp序列建树，生成ASV表和TAX表
# 输出格式：.xls 和 .txt 两个版本
# =====================================================

rm(list=ls())
library(ape)
library(phangorn)

base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
setwd(base_dir)

# 创建输出目录
output_dir <- file.path(base_dir, "建树结果_230bp")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==================== 1. 读取数据 ====================
cat("========== 读取数据 ==========\n")

feature_tax <- read.table("feature_tax.xls", header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE, check.names = FALSE)

rep_seqs_tax <- read.table("rep_seqs_tax.xls", header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE, check.names = FALSE, quote = "")

cat("feature_tax维度:", dim(feature_tax), "\n")
cat("rep_seqs_tax维度:", dim(rep_seqs_tax), "\n")

# ==================== 2. 只保留230bp的序列 ====================
cat("\n========== 筛选230bp序列 ==========\n")

seq_lengths <- nchar(rep_seqs_tax$represent_seq)
keep_idx <- which(seq_lengths == 230)
keep_asv <- rep_seqs_tax$ASV_ID[keep_idx]
keep_sequences <- rep_seqs_tax$represent_seq[keep_idx]

cat("总序列数:", nrow(rep_seqs_tax), "\n")
cat("230bp序列数:", length(keep_asv), "\n")
cat("排除序列数:", nrow(rep_seqs_tax) - length(keep_asv), "\n")

# ==================== 3. 创建FASTA并建树 ====================
cat("\n========== 建树 ==========\n")

# 创建FASTA
fasta_file <- file.path(output_dir, "sequences_230bp.fasta")
fasta_lines <- c()
for(i in 1:length(keep_asv)) {
  fasta_lines <- c(fasta_lines,
                   paste0(">", keep_asv[i]),
                   toupper(keep_sequences[i]))
}
writeLines(fasta_lines, fasta_file)

# 读取序列
dna <- read.dna(fasta_file, format = "fasta")
cat("序列数:", length(labels(dna)), "\n")
cat("序列长度:", ncol(dna), "bp\n")

# 计算距离矩阵
cat("计算距离矩阵...\n")
dist_mat <- dist.dna(dna, model = "TN93", pairwise.deletion = TRUE)

# 构建NJ树
cat("构建NJ树...\n")
tree <- njs(dist_mat)

# 生根
tree_rooted <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)

# 保存树文件（两个版本）
write.tree(tree_rooted, file.path(output_dir, "tree.nwk"))
write.tree(tree_rooted, file.path(output_dir, "tree.txt"))
cat("✓ 树文件已保存\n")

# 清理临时文件
file.remove(fasta_file)

# ==================== 4. 解析TAXONOMY函数 ====================
parse_taxonomy <- function(tax_str) {
  result <- c(Kingdom = NA, Phylum = NA, Class = NA, Order = NA, 
              Family = NA, Genus = NA, Species = NA)
  
  if(is.na(tax_str) || tax_str == "") return(result)
  
  if(grepl("k__", tax_str)) {
    k <- sub(".*(k__[^;]+).*", "\\1", tax_str)
    result["Kingdom"] <- ifelse(k == tax_str, NA, k)
  }
  if(grepl("p__", tax_str)) {
    p <- sub(".*(p__[^;]+).*", "\\1", tax_str)
    result["Phylum"] <- ifelse(p == tax_str, NA, p)
  }
  if(grepl("c__", tax_str)) {
    c <- sub(".*(c__[^;]+).*", "\\1", tax_str)
    result["Class"] <- ifelse(c == tax_str, NA, c)
  }
  if(grepl("o__", tax_str)) {
    o <- sub(".*(o__[^;]+).*", "\\1", tax_str)
    result["Order"] <- ifelse(o == tax_str, NA, o)
  }
  if(grepl("f__", tax_str)) {
    f <- sub(".*(f__[^;]+).*", "\\1", tax_str)
    result["Family"] <- ifelse(f == tax_str, NA, f)
  }
  if(grepl("g__", tax_str)) {
    g <- sub(".*(g__[^;]+).*", "\\1", tax_str)
    result["Genus"] <- ifelse(g == tax_str, NA, g)
  }
  if(grepl("s__", tax_str)) {
    s <- sub(".*(s__[^;]+).*", "\\1", tax_str)
    result["Species"] <- ifelse(s == tax_str, NA, s)
  }
  
  return(result)
}

# ==================== 5. 生成TAX表 ====================
cat("\n========== 生成TAX表 ==========\n")

# 获取保留ASV的分类信息
tax_list <- list()
for(i in 1:length(keep_asv)) {
  asv <- keep_asv[i]
  tax_row <- rep_seqs_tax[rep_seqs_tax$ASV_ID == asv, ]
  if(nrow(tax_row) > 0) {
    tax_parsed <- parse_taxonomy(tax_row$taxonomy[1])
    tax_list[[i]] <- data.frame(
      ASV_ID = asv,
      Kingdom = tax_parsed["Kingdom"],
      Phylum = tax_parsed["Phylum"],
      Class = tax_parsed["Class"],
      Order = tax_parsed["Order"],
      Family = tax_parsed["Family"],
      Genus = tax_parsed["Genus"],
      Species = tax_parsed["Species"],
      stringsAsFactors = FALSE
    )
  }
}

tax_table <- do.call(rbind, tax_list)
cat("TAX表行数:", nrow(tax_table), "\n")

# 保存TAX表（两个版本）
write.table(tax_table, file.path(output_dir, "tax.xls"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tax_table, file.path(output_dir, "tax.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("✓ TAX表已保存\n")

# ==================== 6. 拆分时间点并生成ASV表 ====================
cat("\n========== 拆分时间点生成ASV表 ==========\n")

all_cols <- colnames(feature_tax)
asv_col <- "ASV_ID"

# 定义7个时间点
timepoints <- list(
  DH04 = list(pattern = "^C4DH", name = "DH04"),
  DH07 = list(pattern = "^X7DH", name = "DH07"),
  DH09 = list(pattern = "^Q9DH", name = "DH09"),
  DH10 = list(pattern = "^Q10DH", name = "DH10"),
  WLSH04 = list(pattern = "^C4WLSH", name = "WLSH04"),
  WLSH07 = list(pattern = "^X7WLSH", name = "WLSH07"),
  WLSH10 = list(pattern = "^Q10WLSH", name = "WLSH10")
)

# 创建时间点子目录
for(tp in names(timepoints)) {
  tp_dir <- file.path(output_dir, timepoints[[tp]]$name)
  if(!dir.exists(tp_dir)) dir.create(tp_dir, recursive = TRUE)
}

results_summary <- data.frame(
  时间点 = character(),
  总ASV数 = integer(),
  保留ASV数 = integer(),
  覆盖率 = numeric(),
  总reads = integer(),
  保留reads = integer(),
  reads保留率 = numeric(),
  stringsAsFactors = FALSE
)

for(tp in names(timepoints)) {
  cat("\n处理:", tp, "\n")
  
  pattern <- timepoints[[tp]]$pattern
  sample_cols <- grep(pattern, all_cols, value = TRUE)
  
  if(length(sample_cols) == 0) {
    cat("  警告：没有找到匹配的列，跳过\n")
    next
  }
  
  # 提取该时间点的数据
  extract_cols <- c(asv_col, sample_cols)
  subset_data <- feature_tax[, extract_cols, drop = FALSE]
  
  # 过滤全0行
  row_sums <- rowSums(subset_data[, sample_cols, drop = FALSE], na.rm = TRUE)
  subset_data_all <- subset_data[row_sums > 0, ]
  total_reads_all <- sum(row_sums)
  total_asv_all <- nrow(subset_data_all)
  
  # 只保留230bp的ASV
  asv_ids_all <- subset_data_all[[asv_col]]
  keep_in_time <- asv_ids_all[asv_ids_all %in% keep_asv]
  subset_data_filtered <- subset_data_all[subset_data_all[[asv_col]] %in% keep_in_time, ]
  
  # 计算保留的reads
  row_sums_filtered <- rowSums(subset_data_filtered[, sample_cols, drop = FALSE], na.rm = TRUE)
  total_reads_filtered <- sum(row_sums_filtered)
  total_asv_filtered <- nrow(subset_data_filtered)
  
  cat("  原始ASV:", total_asv_all, "\n")
  cat("  保留ASV:", total_asv_filtered, sprintf("(%.1f%%)", total_asv_filtered/total_asv_all*100), "\n")
  cat("  原始reads:", total_reads_all, "\n")
  cat("  保留reads:", total_reads_filtered, sprintf("(%.1f%%)", total_reads_filtered/total_reads_all*100), "\n")
  
  # 保存ASV表（两个版本）
  tp_name <- timepoints[[tp]]$name
  tp_dir <- file.path(output_dir, tp_name)
  
  write.table(subset_data_filtered, file.path(tp_dir, "ASV.xls"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(subset_data_filtered, file.path(tp_dir, "ASV.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 同时保存该时间点的TAX表（只包含该时间点的ASV）
  tax_time <- tax_table[tax_table$ASV_ID %in% keep_in_time, ]
  write.table(tax_time, file.path(tp_dir, "tax.xls"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(tax_time, file.path(tp_dir, "tax.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 保存该时间点的树（从完整树提取子树）
  if(length(keep_in_time) >= 2) {
    subtree <- keep.tip(tree_rooted, keep_in_time)
    write.tree(subtree, file.path(tp_dir, "tree.nwk"))
    write.tree(subtree, file.path(tp_dir, "tree.txt"))
  }
  
  cat("  ✓ 已保存到:", tp_dir, "\n")
  
  # 记录结果
  results_summary <- rbind(results_summary, data.frame(
    时间点 = tp_name,
    总ASV数 = total_asv_all,
    保留ASV数 = total_asv_filtered,
    覆盖率 = round(total_asv_filtered/total_asv_all*100, 2),
    总reads = total_reads_all,
    保留reads = total_reads_filtered,
    reads保留率 = round(total_reads_filtered/total_reads_all*100, 2),
    stringsAsFactors = FALSE
  ))
}

# ==================== 7. 保存汇总报告 ====================
cat("\n\n========================================\n")
cat("========== 处理完成 ==========\n")
cat("结果保存在:", output_dir, "\n\n")

print(results_summary)

# 保存汇总报告
write.table(results_summary, file.path(output_dir, "汇总报告.xls"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(results_summary, file.path(output_dir, "汇总报告.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 保存完整TAX表
write.table(tax_table, file.path(output_dir, "tax_complete.xls"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tax_table, file.path(output_dir, "tax_complete.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==================== 8. 打印目录结构 ====================
cat("\n\n输出目录结构:\n")
cat(output_dir, "/\n")
cat("  ├─ tree.nwk\n")
cat("  ├─ tree.txt\n")
cat("  ├─ tax_complete.xls\n")
cat("  ├─ tax_complete.txt\n")
cat("  ├─ 汇总报告.xls\n")
cat("  ├─ 汇总报告.txt\n")

for(tp in names(timepoints)) {
  tp_name <- timepoints[[tp]]$name
  tp_dir <- file.path(output_dir, tp_name)
  if(dir.exists(tp_dir)) {
    cat("\n  ", tp_name, "/\n")
    cat("    ├─ ASV.xls\n")
    cat("    ├─ ASV.txt\n")
    cat("    ├─ tax.xls\n")
    cat("    ├─ tax.txt\n")
    cat("    ├─ tree.nwk\n")
    cat("    ├─ tree.txt\n")
  }
}

cat("\n========== 全部完成 ==========\n")


# =====================================================
# 生成完整的metadata.xls文件（修正Group2格式）
# Group2格式：季节代码+月份+类型，如 C4DH, X7DH, Q9DH, Q10DH, C4WLSH, X7WLSH, Q10WLSH
# =====================================================

rm(list=ls())

base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
setwd(base_dir)

# ==================== 1. 定义采样点坐标 ====================

# DH点坐标（9个）
dh_coords <- data.frame(
  SiteID = paste0("DH", sprintf("%02d", 1:9)),
  Longitude = c(112.65891991, 112.66657677, 112.66770020, 112.67759617, 
                112.69162299, 112.71829538, 112.70754482, 112.70849324, 112.73368777),
  Latitude = c(40.57169287, 40.55624468, 40.59140684, 40.57981620, 
               40.56555673, 40.57062358, 40.58574064, 40.60106757, 40.58459969),
  Type = "DH",
  stringsAsFactors = FALSE
)

# WLSH点坐标（21个）
wlsh_coords <- data.frame(
  SiteID = paste0("WLSH", 1:21),
  Longitude = c(108.86489868, 108.90403748, 108.91880035, 108.91777039, 
                108.92910004, 108.93510818, 108.92978668, 108.90678406, 
                108.87508517, 108.87073517, 108.87382507, 108.86713028, 
                108.84086609, 108.80584717, 108.80687714, 108.79829407, 
                108.78318787, 108.76258850, 108.74130249, 108.71658325, 
                108.84120941),
  Latitude = c(41.11039943, 41.09875804, 41.07288091, 41.02705341, 
               41.00503451, 40.98132359, 40.95527062, 40.93685883, 
               40.91173496, 40.99246749, 40.94775093, 40.92635404, 
               40.91792306, 40.91850677, 40.89898185, 40.87951606, 
               40.85874624, 40.84082704, 40.81926564, 40.79379857, 
               40.87354541),
  Type = "WLSH",
  stringsAsFactors = FALSE
)

# 合并所有点
all_sites <- rbind(dh_coords, wlsh_coords)

# ==================== 2. 定义季节代码和Group2映射 ====================
# 季节映射：前缀 -> 季节代码 + 月份 + 类型
season_mapping <- list(
  C4 = list(prefix = "C4", month = "04", season = "Spring", group2_suffix = "DH", group2_wlsh = "WLSH"),
  X7 = list(prefix = "X7", month = "07", season = "Summer", group2_suffix = "DH", group2_wlsh = "WLSH"),
  Q9 = list(prefix = "Q9", month = "09", season = NA, group2_suffix = "DH", group2_wlsh = "WLSH"),
  Q10 = list(prefix = "Q10", month = "10", season = "Autumn", group2_suffix = "DH", group2_wlsh = "WLSH")
)

# 定义每个时间点有哪些采样点
# DH样本
dh_samples <- list(
  C4 = c("DH01", "DH02", "DH03", "DH04", "DH05", "DH06", "DH07", "DH09"),
  X7 = c("DH04", "DH05", "DH06"),
  Q9 = c("DH01", "DH02", "DH03", "DH04", "DH05", "DH06", "DH07", "DH09"),
  Q10 = c("DH02", "DH03", "DH04", "DH05", "DH06", "DH07", "DH08", "DH09")
)

# WLSH样本
wlsh_samples <- list(
  C4 = paste0("WLSH", c(1, 3:21)),
  X7 = paste0("WLSH", c(2, 4:21)),
  Q10 = paste0("WLSH", c(1:9, 12:21))
)

# ==================== 3. 生成metadata ====================

metadata <- data.frame(
  SampleID = character(),
  SiteID = character(),
  Group1 = character(),
  Group2 = character(),
  Season = character(),
  Time = character(),
  Longitude = numeric(),
  Latitude = numeric(),
  stringsAsFactors = FALSE
)

# 处理DH样本
for(prefix in names(dh_samples)) {
  season_info <- season_mapping[[prefix]]
  sites <- dh_samples[[prefix]]
  
  for(site in sites) {
    sample_id <- paste0(prefix, site)
    group2 <- paste0(prefix, "DH")  # C4DH, X7DH, Q9DH, Q10DH
    
    # 查找坐标
    coord <- all_sites[all_sites$SiteID == site, ]
    
    metadata <- rbind(metadata, data.frame(
      SampleID = sample_id,
      SiteID = site,
      Group1 = "DH",
      Group2 = group2,
      Season = season_info$season,
      Time = season_info$month,
      Longitude = if(nrow(coord) > 0) coord$Longitude else NA,
      Latitude = if(nrow(coord) > 0) coord$Latitude else NA,
      stringsAsFactors = FALSE
    ))
  }
}

# 处理WLSH样本
for(prefix in names(wlsh_samples)) {
  season_info <- season_mapping[[prefix]]
  sites <- wlsh_samples[[prefix]]
  
  for(site in sites) {
    sample_id <- paste0(prefix, site)
    group2 <- paste0(prefix, "WLSH")  # C4WLSH, X7WLSH, Q10WLSH
    
    # 查找坐标
    coord <- all_sites[all_sites$SiteID == site, ]
    
    metadata <- rbind(metadata, data.frame(
      SampleID = sample_id,
      SiteID = site,
      Group1 = "WLSH",
      Group2 = group2,
      Season = season_info$season,
      Time = season_info$month,
      Longitude = if(nrow(coord) > 0) coord$Longitude else NA,
      Latitude = if(nrow(coord) > 0) coord$Latitude else NA,
      stringsAsFactors = FALSE
    ))
  }
}

# 按Group2和SampleID排序
metadata <- metadata[order(metadata$Group2, metadata$SampleID), ]

# ==================== 4. 查看结果 ====================
cat("========== Metadata统计 ==========\n")
cat("总样本数:", nrow(metadata), "\n")
cat("\n各Group2样本数:\n")
print(table(metadata$Group2))

cat("\n各季节样本数:\n")
print(table(metadata$Season, useNA = "ifany"))

cat("\n前20行:\n")
print(head(metadata, 20))

# 检查缺失坐标
missing_coord <- metadata[is.na(metadata$Longitude), ]
if(nrow(missing_coord) > 0) {
  cat("\n警告：以下样本缺少坐标:\n")
  print(missing_coord[, c("SampleID", "SiteID")])
}

# ==================== 5. 保存文件 ====================
write.table(metadata, file.path(base_dir, "metadata.xls"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(metadata, file.path(base_dir, "metadata.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n✓ metadata.xls 和 metadata.txt 已保存到:", base_dir, "\n")

# ==================== 6. 输出预览 ====================
cat("\n========== metadata.xls 预览 ==========\n")
cat("SampleID\tSiteID\tGroup1\tGroup2\tSeason\tTime\tLongitude\tLatitude\n")
for(i in 1:min(10, nrow(metadata))) {
  cat(metadata$SampleID[i], "\t",
      metadata$SiteID[i], "\t",
      metadata$Group1[i], "\t",
      metadata$Group2[i], "\t",
      metadata$Season[i], "\t",
      metadata$Time[i], "\t",
      metadata$Longitude[i], "\t",
      metadata$Latitude[i], "\n")
}

cat("\n========== 完成 ==========\n")

# =====================================================
# Source: 分析数据.R
# =====================================================

# =====================================================
# 为每个时间点生成两个版本的环境文件
# env.txt: 所有环境数据（原文件，保持不变）
# env_asv.txt: 只包含有ASV对应的环境数据（即使内容相同也生成）
# =====================================================

rm(list=ls())

base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
setwd(data_dir)

# 定义时间点
timepoints <- c("DH04", "DH07", "DH09", "DH10", "WLSH04", "WLSH07", "WLSH10")
has_env <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
names(has_env) <- timepoints

# 提取SiteID的函数
extract_site_id <- function(sample_name) {
  prefixes <- c("C4", "X7", "Q9", "Q10")
  for(prefix in prefixes) {
    if(startsWith(sample_name, prefix)) {
      return(substr(sample_name, nchar(prefix) + 1, nchar(sample_name)))
    }
  }
  return(sample_name)
}

# ==================== 为每个时间点生成文件 ====================
cat("========== 生成 env.txt 和 env_asv.txt ==========\n")

for(tp in timepoints) {
  cat("\n", tp, ":\n")
  
  tp_dir <- file.path(data_dir, tp)
  if(!dir.exists(tp_dir)) {
    cat("  目录不存在，跳过\n")
    next
  }
  
  # 读取原始env.txt
  env_file <- file.path(tp_dir, "env.txt")
  if(!file.exists(env_file)) {
    cat("  env.txt不存在，跳过\n")
    next
  }
  
  env <- read.table(env_file, header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE, check.names = FALSE)
  colnames(env)[1] <- "SiteID"
  
  # 统一SiteID格式（去掉连字符）
  env$SiteID <- gsub("-", "", env$SiteID)
  
  cat("  原始环境样本数:", nrow(env), "\n")
  
  # 读取ASV样本（获取有ASV的站点）
  asv_file <- file.path(tp_dir, "ASV.txt")
  asv_sites <- c()
  
  if(file.exists(asv_file)) {
    asv <- read.table(asv_file, header = TRUE, sep = "\t", nrows = 1)
    sample_cols <- colnames(asv)[!colnames(asv) %in% c("ASV_ID", "taxonomy")]
    asv_sites <- sapply(sample_cols, extract_site_id)
    cat("  ASV样本数:", length(asv_sites), "\n")
  } else {
    cat("  ASV.txt不存在\n")
  }
  
  # 生成env_asv.txt
  if(length(asv_sites) > 0) {
    env_asv <- env[env$SiteID %in% asv_sites, ]
  } else {
    # 如果没有ASV数据，创建空数据框（但保留列结构）
    env_asv <- env[0, ]
    cat("  警告：无ASV数据，生成空文件\n")
  }
  
  # 无论内容是否相同，都生成env_asv.txt
  env_asv_file <- file.path(tp_dir, "env_asv.txt")
  write.table(env_asv, env_asv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  ✓ env_asv.txt 已生成:", nrow(env_asv), "行\n")
  
  # 显示信息
  if(nrow(env) == nrow(env_asv)) {
    cat("  （与环境数据完全一致）\n")
  } else {
    excluded <- env[!env$SiteID %in% asv_sites, ]
    cat("  排除的样本:", paste(excluded$SiteID, collapse=", "), "\n")
  }
  
  # 确保原始env.txt存在（它本来就存在，不需要修改）
  cat("  ✓ env.txt 保持不变\n")
}

# ==================== 汇总 ====================
cat("\n========== 完成 ==========\n")
cat("每个时间点文件夹内现在有:\n")
cat("  - env.txt (所有环境数据)\n")
cat("  - env_asv.txt (ASV对应的环境数据)\n")

# =====================================================
# Source: 分析数据个数.R
# =====================================================

# =====================================================
# 诊断：找出所有ASV和环境不匹配的样本
# =====================================================

rm(list=ls())

base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
setwd(data_dir)

# 读取已生成的共同样本数据
common_samples <- read.table("common_env_data.txt", header = TRUE, sep = "\t", 
                             stringsAsFactors = FALSE, check.names = FALSE)

cat("========== 共同样本（有ASV+环境）==========\n")
cat("数量:", nrow(common_samples), "\n\n")

# 按时间点统计共同样本
cat("按时间点统计:\n")
print(table(common_samples$TimePoint))

# ==================== 读取各时间点的ASV样本 ====================
cat("\n\n========== 各时间点ASV样本 ==========\n")

timepoints <- c("DH04", "DH07", "DH09", "DH10", "WLSH04", "WLSH07", "WLSH10")
has_env <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
names(has_env) <- timepoints

asv_counts <- data.frame()

for(tp in timepoints) {
  tp_dir <- file.path(data_dir, tp)
  asv_file <- file.path(tp_dir, "ASV.txt")
  
  if(file.exists(asv_file)) {
    asv <- read.table(asv_file, header = TRUE, sep = "\t", nrows = 1)
    sample_cols <- colnames(asv)[!colnames(asv) %in% c("ASV_ID", "taxonomy")]
    
    asv_counts <- rbind(asv_counts, data.frame(
      TimePoint = tp,
      ASV样本数 = length(sample_cols),
      ASV样本 = paste(sample_cols, collapse=", "),
      stringsAsFactors = FALSE
    ))
  }
}

print(asv_counts)

# ==================== 读取各时间点的环境样本 ====================
cat("\n\n========== 各时间点环境样本 ==========\n")

env_counts <- data.frame()

for(tp in timepoints) {
  if(has_env[tp]) {
    env_file <- file.path(data_dir, tp, "env.txt")
    if(file.exists(env_file)) {
      env <- read.table(env_file, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, check.names = FALSE)
      colnames(env)[1] <- "SiteID"
      env$SiteID_clean <- gsub("-", "", env$SiteID)
      
      env_counts <- rbind(env_counts, data.frame(
        TimePoint = tp,
        环境样本数 = nrow(env),
        环境SiteID = paste(env$SiteID_clean, collapse=", "),
        stringsAsFactors = FALSE
      ))
    }
  } else {
    env_counts <- rbind(env_counts, data.frame(
      TimePoint = tp,
      环境样本数 = 0,
      环境SiteID = "无",
      stringsAsFactors = FALSE
    ))
  }
}

print(env_counts)

# ==================== 找出每个时间点的具体不匹配 ====================
cat("\n\n========== 各时间点匹配详情 ==========\n")

for(tp in timepoints) {
  if(tp == "DH09") {
    cat("\n", tp, ": 被排除分析\n")
    next
  }
  
  cat("\n", tp, ":\n")
  
  # 读取ASV样本
  tp_dir <- file.path(data_dir, tp)
  asv_file <- file.path(tp_dir, "ASV.txt")
  
  if(!file.exists(asv_file)) {
    cat("  ASV文件不存在\n")
    next
  }
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t", nrows = 1)
  asv_samples <- colnames(asv)[!colnames(asv) %in% c("ASV_ID", "taxonomy")]
  
  # 提取SiteID
  extract_site_id <- function(sample_name) {
    prefixes <- c("C4", "X7", "Q9", "Q10")
    for(prefix in prefixes) {
      if(startsWith(sample_name, prefix)) {
        return(substr(sample_name, nchar(prefix) + 1, nchar(sample_name)))
      }
    }
    return(sample_name)
  }
  
  asv_sites <- sapply(asv_samples, extract_site_id)
  
  # 读取环境数据
  if(has_env[tp]) {
    env_file <- file.path(data_dir, tp, "env.txt")
    if(file.exists(env_file)) {
      env <- read.table(env_file, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, check.names = FALSE)
      colnames(env)[1] <- "SiteID_raw"
      env_sites <- gsub("-", "", env$SiteID_raw)
      
      # 找出匹配和未匹配
      matched <- asv_sites[asv_sites %in% env_sites]
      asv_only <- asv_sites[!asv_sites %in% env_sites]
      env_only <- env_sites[!env_sites %in% asv_sites]
      
      cat("  ASV样本数:", length(asv_samples), "\n")
      cat("  环境样本数:", length(env_sites), "\n")
      cat("  匹配数:", length(matched), "\n")
      
      if(length(asv_only) > 0) {
        cat("  只有ASV没环境:", paste(asv_only, collapse=", "), "\n")
      }
      if(length(env_only) > 0) {
        cat("  只有环境没ASV:", paste(env_only, collapse=", "), "\n")
      }
    }
  } else {
    cat("  无环境数据\n")
  }
}

# ==================== 汇总所有不匹配 ====================
cat("\n\n========== 不匹配汇总 ==========\n")

all_asv_only <- c()
all_env_only <- c()

for(tp in timepoints) {
  if(tp == "DH09") next
  
  tp_dir <- file.path(data_dir, tp)
  asv_file <- file.path(tp_dir, "ASV.txt")
  
  if(!file.exists(asv_file)) next
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t", nrows = 1)
  asv_samples <- colnames(asv)[!colnames(asv) %in% c("ASV_ID", "taxonomy")]
  asv_sites <- sapply(asv_samples, extract_site_id)
  
  if(has_env[tp]) {
    env_file <- file.path(data_dir, tp, "env.txt")
    if(file.exists(env_file)) {
      env <- read.table(env_file, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, check.names = FALSE)
      colnames(env)[1] <- "SiteID_raw"
      env_sites <- gsub("-", "", env$SiteID_raw)
      
      asv_only <- asv_sites[!asv_sites %in% env_sites]
      env_only <- env_sites[!env_sites %in% asv_sites]
      
      if(length(asv_only) > 0) {
        all_asv_only <- c(all_asv_only, paste0(tp, ":", asv_only))
      }
      if(length(env_only) > 0) {
        all_env_only <- c(all_env_only, paste0(tp, ":", env_only))
      }
    }
  }
}

cat("\n只有ASV没有环境的样本:\n")
if(length(all_asv_only) > 0) {
  cat(paste(all_asv_only, collapse="\n"), "\n")
} else {
  cat("无\n")
}

cat("\n只有环境没有ASV的样本:\n")
if(length(all_env_only) > 0) {
  cat(paste(all_env_only, collapse="\n"), "\n")
} else {
  cat("无\n")
}
cat("\n========== 完成 ==========\n")
