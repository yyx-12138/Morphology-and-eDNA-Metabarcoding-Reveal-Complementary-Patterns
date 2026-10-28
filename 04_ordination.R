# 04 Ordination (PCoA / RDA)


# =====================================================
# Source: 图4C_H.R
# =====================================================

# =====================================================
# PCoA分析 - 形态学和eDNA（修复索引错误）
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(ape)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir_morph <- file.path(base_dir, "输出图片/Fig4CEG")
output_dir_edna <- file.path(base_dir, "输出图片/Fig4DFH")
if(!dir.exists(output_dir_morph)) dir.create(output_dir_morph, recursive = TRUE)
if(!dir.exists(output_dir_edna)) dir.create(output_dir_edna, recursive = TRUE)

# 配色
lake_colors <- c("DH" = "#C82423", "WLSH" = "#5A7D9A")

# =====================================================
# 通用函数：按季节独立计算PCoA，合并后用facet_wrap
# =====================================================

plot_pcoa_facet <- function(data_list, metadata, title_name, output_dir, file_prefix) {
  
  seasons <- c("Spring", "Summer", "Autumn")
  all_coords <- data.frame()
  var_info <- data.frame(Season = character(), PC1_var = numeric(), PC2_var = numeric())
  
  for(season in seasons) {
    cat("  处理季节:", season, "\n")
    
    # 筛选该季节的样本
    season_samples <- metadata %>% filter(Season == season) %>% pull(Sample)
    season_data <- data_list[names(data_list) %in% season_samples]
    
    if(length(season_data) < 3) {
      cat("    样本数不足，跳过\n")
      next
    }
    
    # 获取所有属名
    all_genus <- unique(unlist(lapply(season_data, names)))
    
    # 构建矩阵
    season_matrix <- matrix(0, nrow = length(season_data), ncol = length(all_genus))
    rownames(season_matrix) <- names(season_data)
    colnames(season_matrix) <- all_genus
    
    for(sample_name in names(season_data)) {
      abu <- season_data[[sample_name]]
      for(genus in names(abu)) {
        season_matrix[sample_name, genus] <- abu[genus]
      }
    }
    
    # 移除全0列
    col_sums <- colSums(season_matrix)
    season_matrix <- season_matrix[, col_sums > 0, drop = FALSE]
    
    if(nrow(season_matrix) < 3 || ncol(season_matrix) < 2) {
      cat("    数据不足，跳过\n")
      next
    }
    
    # 计算Bray-Curtis距离
    dist_matrix <- vegdist(season_matrix, method = "bray")
    
    # PCoA
    pcoa_result <- pcoa(dist_matrix)
    
    # 提取坐标
    coords <- as.data.frame(pcoa_result$vectors[, 1:2])
    colnames(coords) <- c("PC1", "PC2")
    coords$Sample <- rownames(coords)
    coords <- merge(coords, metadata, by = "Sample")
    
    # 添加季节列
    coords$Season <- season
    
    # 解释方差
    var1 <- round(pcoa_result$values$Relative_eig[1] * 100, 2)
    var2 <- round(pcoa_result$values$Relative_eig[2] * 100, 2)
    
    cat("    样本数:", nrow(coords), ", PC1:", var1, "%, PC2:", var2, "%\n")
    
    var_info <- rbind(var_info, data.frame(Season = season, PC1_var = var1, PC2_var = var2))
    all_coords <- rbind(all_coords, coords)
  }
  
  if(nrow(all_coords) == 0) {
    cat("  没有有效的季节数据\n")
    return(NULL)
  }
  
  # 设置因子顺序
  all_coords$Season <- factor(all_coords$Season, levels = c("Spring", "Summer", "Autumn"))
  
  # 创建分面标签（包含季节和解释方差）
  season_labels <- c(
    Spring = paste0("Spring\n(PC1: ", var_info$PC1_var[var_info$Season == "Spring"], 
                    "%, PC2: ", var_info$PC2_var[var_info$Season == "Spring"], "%)"),
    Summer = paste0("Summer\n(PC1: ", var_info$PC1_var[var_info$Season == "Summer"], 
                    "%, PC2: ", var_info$PC2_var[var_info$Season == "Summer"], "%)"),
    Autumn = paste0("Autumn\n(PC1: ", var_info$PC1_var[var_info$Season == "Autumn"], 
                    "%, PC2: ", var_info$PC2_var[var_info$Season == "Autumn"], "%)")
  )
  
  # 绘图 - 使用真正的facet_wrap
  p <- ggplot(all_coords, aes(x = PC1, y = PC2, color = Lake, fill = Lake)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 1, linewidth = 0.8) +
    scale_color_manual(values = lake_colors, name = "Lake") +
    scale_fill_manual(values = lake_colors, name = "Lake") +
    facet_wrap(~Season, nrow = 1, ncol = 3, scales = "free",
               labeller = labeller(Season = season_labels)) +
    labs(x = NULL, y = NULL, title = title_name) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 10, margin = margin(5, 0, 5, 0)),
      panel.grid = element_blank(),
      axis.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # 保存
  pdf(file.path(output_dir, paste0(file_prefix, "_PCoA.pdf")), 
      width = 12, height = 5.5, useDingbats = FALSE)
  print(p)
  dev.off()
  
  ggsave(file.path(output_dir, paste0(file_prefix, "_PCoA.png")), 
         p, width = 12, height = 5.5, dpi = 300, bg = "white")
  
  # 保存坐标
  for(season in unique(all_coords$Season)) {
    season_coords <- all_coords %>% filter(Season == season)
    write.csv(season_coords, 
              file.path(output_dir, paste0(file_prefix, "_PCoA_", season, ".csv")), 
              row.names = FALSE)
  }
  
  cat("  ✅ 已保存到:", output_dir, "\n")
  
  return(all_coords)
}

# =====================================================
# 第一部分：形态学数据
# =====================================================

cat("\n", rep("=", 60), "\n", sep="")
cat("【形态学数据】\n")
cat(rep("=", 60), "\n", sep="")

timepoints <- c("DH04", "DH07", "DH10", "WLSH04", "WLSH07", "WLSH10")
season_map <- c("Spring", "Summer", "Autumn", "Spring", "Summer", "Autumn")
names(season_map) <- timepoints

# 存储形态学数据
morpho_data <- list()
morpho_metadata <- data.frame()

for(tp in timepoints) {
  file_path <- file.path(data_dir, tp, "morpho.txt")
  if(!file.exists(file_path)) next
  
  morpho <- read.table(file_path, 
                       header = TRUE, 
                       sep = "\t", 
                       check.names = FALSE,
                       row.names = 1,
                       stringsAsFactors = FALSE)
  
  morpho <- as.data.frame(lapply(morpho, function(x) as.numeric(as.character(x))))
  morpho[is.na(morpho)] <- 0
  
  lake <- ifelse(grepl("^DH", tp), "DH", "WLSH")
  season <- season_map[tp]
  
  for(sample_name in colnames(morpho)) {
    abu <- morpho[[sample_name]]
    names(abu) <- rownames(morpho)
    abu <- abu[abu > 0]
    
    if(length(abu) > 0) {
      morpho_metadata <- rbind(morpho_metadata, data.frame(
        Sample = sample_name,
        Lake = lake,
        Season = season,
        TimePoint = tp,
        stringsAsFactors = FALSE
      ))
      morpho_data[[sample_name]] <- abu
    }
  }
}

cat("形态学总样本数:", length(morpho_data), "\n")

# 计算形态学PCoA
morpho_coords <- plot_pcoa_facet(morpho_data, morpho_metadata, 
                                 "Morphology", output_dir_morph, "Fig4C_Morphology")

# =====================================================
# 第二部分：eDNA数据
# =====================================================

cat("\n", rep("=", 60), "\n", sep="")
cat("【eDNA数据】\n")
cat(rep("=", 60), "\n", sep="")

timepoints <- c("DH04","DH07","DH10","WLSH04","WLSH07","WLSH10")

edna_data <- list()
edna_metadata <- data.frame()

for(tp in timepoints) {
  asv_file <- file.path(data_dir, tp, "ASV.txt")
  tax_file <- file.path(data_dir, tp, "tax.txt")
  
  if(!file.exists(asv_file) | !file.exists(tax_file)) next
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE)
  sample_cols <- setdiff(colnames(asv), c("ASV_ID", "taxonomy"))
  if(length(sample_cols) == 0) next
  
  tax <- read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    na.strings = c("NA", "", "uncultured"))
  
  # 提取属水平，排除Unclassified
  tax$Genus <- gsub("^g__", "", tax$Genus)
  tax$Genus[is.na(tax$Genus) | tax$Genus == ""] <- "Unclassified"
  
  # 排除Unclassified
  tax_clean <- tax[tax$Genus != "Unclassified", ]
  
  if(nrow(tax_clean) == 0) next
  
  # 创建ASV到属的映射
  asv_to_genus <- tax_clean$Genus
  names(asv_to_genus) <- tax_clean$ASV_ID
  
  # 只保留有分类信息的ASV
  asv_subset <- asv[asv$ASV_ID %in% names(asv_to_genus), ]
  
  if(nrow(asv_subset) == 0) next
  
  asv_matrix <- as.matrix(asv_subset[, sample_cols])
  rownames(asv_matrix) <- asv_subset$ASV_ID
  
  lake <- ifelse(grepl("DH", tp), "DH", "WLSH")
  season <- ifelse(grepl("04$", tp), "Spring",
                   ifelse(grepl("07$", tp), "Summer", "Autumn"))
  
  for(sample_name in sample_cols) {
    sample_abu <- asv_matrix[, sample_name]
    names(sample_abu) <- rownames(asv_matrix)
    
    # 按属聚合
    genus_abu <- tapply(sample_abu, asv_to_genus[names(sample_abu)], sum)
    genus_abu <- genus_abu[genus_abu > 0]
    
    if(length(genus_abu) > 0) {
      edna_metadata <- rbind(edna_metadata, data.frame(
        Sample = sample_name,
        Lake = lake,
        Season = season,
        TimePoint = tp,
        stringsAsFactors = FALSE
      ))
      edna_data[[sample_name]] <- genus_abu
    }
  }
}

cat("eDNA总样本数:", length(edna_data), "\n")

if(length(edna_data) > 0) {
  # 计算eDNA PCoA
  edna_coords <- plot_pcoa_facet(edna_data, edna_metadata, 
                                 "eDNA", output_dir_edna, "Fig4D_eDNA")
} else {
  cat("⚠️ 没有有效的eDNA数据\n")
}

# =====================================================
# 完成
# =====================================================

cat("\n", rep("=", 60), "\n", sep="")
cat("【PCoA 分析完成】\n")
cat(rep("=", 60), "\n", sep="")

cat("\n形态学 PCoA:\n")
cat("  输出文件夹: Fig4CEG\n")
cat("  - Fig4C_Morphology_PCoA.pdf\n")
cat("  - Fig4C_Morphology_PCoA.png\n")

cat("\neDNA PCoA:\n")
cat("  输出文件夹: Fig4DFH\n")
cat("  - Fig4D_eDNA_PCoA.pdf\n")
cat("  - Fig4D_eDNA_PCoA.png\n")

cat("\n图形特点:\n")
cat("  - ✅ 真正的facet_wrap分面图\n")
cat("  - ✅ 每个季节独立计算PCoA\n")
cat("  - ✅ 灰色背景的strip标签\n")
cat("  - ✅ 每个分面独立的坐标轴（scales = free）\n")
cat("  - ✅ 距离算法: Bray-Curtis\n")
cat("  - ✅ 数据水平: 属水平\n")
cat("  - ✅ 置信椭圆: 95%\n")
cat("  - ✅ eDNA已排除Unclassified\n")

# =====================================================
# Source: 图4IJK.R
# =====================================================

# =====================================================
# 形态学 RDA 分面图 —— 输出两张图
#  图1：只用原始环境因子（去掉 WQI、TLI）
#  图2：只用 WQI 和 TLI
# 两张图绘图设置完全一致
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ========== 路径 & 基础参数 ==========
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir_morph <- file.path(base_dir, "输出图片/Fig4IJK_Morph_分WQI版本")
if(!dir.exists(output_dir_morph)) dir.create(output_dir_morph, recursive = TRUE)

lake_colors <- c("DH" = "#C82423", "WLSH" = "#5A7D9A")

timepoints <- c("DH04", "DH07", "DH10", "WLSH04", "WLSH07", "WLSH10")
season_map <- c("Spring", "Summer", "Autumn", "Spring", "Summer", "Autumn")
names(season_map) <- timepoints
season_levels <- c("Spring", "Summer", "Autumn")

min_species <- 1

# ========== WQI 单项评分函数 ==========
get_Ci_TN <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.1) return(100)
  else if(value < 0.2) return(90)
  else if(value < 0.35) return(80)
  else if(value < 0.5) return(70)
  else if(value < 0.75) return(60)
  else if(value < 1) return(50)
  else if(value < 1.25) return(40)
  else if(value < 1.5) return(30)
  else if(value < 1.75) return(20)
  else if(value <= 2) return(10)
  else return(0)
}
get_Ci_TP <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.01) return(100)
  else if(value < 0.02) return(90)
  else if(value < 0.05) return(80)
  else if(value < 0.1) return(70)
  else if(value < 0.15) return(60)
  else if(value < 0.2) return(50)
  else if(value < 0.25) return(40)
  else if(value < 0.3) return(30)
  else if(value < 0.35) return(20)
  else if(value <= 0.4) return(10)
  else return(0)
}
get_Ci_NO3 <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.5) return(100)
  else if(value < 1) return(90)
  else if(value < 1.5) return(80)
  else if(value < 2) return(70)
  else if(value < 3) return(60)
  else if(value < 4) return(50)
  else if(value < 5) return(40)
  else if(value < 6) return(30)
  else if(value < 8) return(20)
  else if(value <= 10) return(10)
  else return(0)
}
get_Ci_NH4 <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.01) return(100)
  else if(value < 0.05) return(90)
  else if(value < 0.1) return(80)
  else if(value < 0.2) return(70)
  else if(value < 0.3) return(60)
  else if(value < 0.4) return(50)
  else if(value < 0.5) return(40)
  else if(value < 0.75) return(30)
  else if(value < 1) return(20)
  else if(value <= 1.25) return(10)
  else return(0)
}
get_Ci_Chla <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 1) return(100)
  else if(value < 4) return(90)
  else if(value < 7) return(80)
  else if(value < 10) return(70)
  else if(value < 15) return(60)
  else if(value < 20) return(50)
  else if(value < 30) return(40)
  else if(value < 40) return(30)
  else if(value < 50) return(20)
  else if(value <= 65) return(10)
  else return(0)
}
get_Ci_CODcr <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 15) return(100)
  else if(value < 16) return(90)
  else if(value < 18) return(80)
  else if(value < 19) return(70)
  else if(value < 20) return(60)
  else if(value < 25) return(50)
  else if(value < 30) return(40)
  else if(value < 35) return(30)
  else if(value < 37) return(20)
  else if(value <= 40) return(10)
  else return(0)
}
get_Ci_DO <- function(value) {
  if(is.na(value)) return(NA)
  if(value >= 9) return(100)
  else if(value > 8) return(90)
  else if(value > 7) return(80)
  else if(value > 6) return(70)
  else if(value > 5) return(60)
  else if(value > 4) return(50)
  else if(value > 3.5) return(40)
  else if(value > 3) return(30)
  else if(value > 2) return(20)
  else if(value >= 1) return(10)
  else return(0)
}
get_Ci_pH <- function(value) {
  if(is.na(value)) return(NA)
  if(value == 7) return(100)
  else if(value >= 7 & value <= 8) return(90)
  else if(value >= 7 & value <= 8.5) return(80)
  else if(value >= 6.5 & value <= 9) return(70)
  else if(value >= 6 & value <= 9.5) return(60)
  else if(value >= 5 & value <= 10) return(50)
  else if(value >= 4 & value <= 11) return(40)
  else if(value >= 3 & value <= 12) return(30)
  else if(value >= 2 & value <= 13) return(20)
  else if(value >= 1 & value <= 14) return(10)
  else return(0)
}
get_Ci_Temp <- function(value) {
  if(is.na(value)) return(NA)
  if(value >= 16 & value <= 21) return(100)
  else if((value >= 15 & value < 16) | (value > 21 & value <= 22)) return(90)
  else if((value >= 14 & value < 15) | (value > 22 & value <= 24)) return(80)
  else if((value >= 12 & value < 14) | (value > 24 & value <= 26)) return(70)
  else if((value >= 10 & value < 12) | (value > 26 & value <= 28)) return(60)
  else if((value >= 5 & value < 10) | (value > 28 & value <= 30)) return(50)
  else if((value >= 0 & value < 5) | (value > 30 & value <= 32)) return(40)
  else if((value >= -2 & value < 0) | (value > 32 & value <= 36)) return(30)
  else if((value >= -4 & value < -2) | (value > 36 & value <= 40)) return(20)
  else if((value >= -6 & value < -4) | (value > 40 & value <= 45)) return(10)
  else return(0)
}

# ========== WQI & TLI 综合计算函数 ==========
calculate_WQI <- function(row) {
  total_weight <- 0; weighted_sum <- 0
  tn_val <- as.numeric(row["TN"])
  if(!is.na(tn_val)) { ci <- get_Ci_TN(tn_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  tp_val <- as.numeric(row["TP"])
  if(!is.na(tp_val)) { ci <- get_Ci_TP(tp_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 1; total_weight <- total_weight + 1 } }
  no3_val <- as.numeric(row["NO3--N"])
  if(!is.na(no3_val)) { ci <- get_Ci_NO3(no3_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 2; total_weight <- total_weight + 2 } }
  nh4_val <- as.numeric(row["NH4+-N"])
  if(!is.na(nh4_val)) { ci <- get_Ci_NH4(nh4_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  chla_val <- as.numeric(row["Chl-a"])
  if(!is.na(chla_val)) { ci <- get_Ci_Chla(chla_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  cod_val <- as.numeric(row["CODcr"])
  if(!is.na(cod_val)) { ci <- get_Ci_CODcr(cod_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  do_val <- as.numeric(row["DO"])
  if(!is.na(do_val)) { ci <- get_Ci_DO(do_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 4; total_weight <- total_weight + 4 } }
  ph_val <- as.numeric(row["Ph"])
  if(!is.na(ph_val)) { ci <- get_Ci_pH(ph_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 1; total_weight <- total_weight + 1 } }
  temp_val <- as.numeric(row["Temp"])
  if(!is.na(temp_val)) { ci <- get_Ci_Temp(temp_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 1; total_weight <- total_weight + 1 } }
  if(total_weight > 0) return(weighted_sum / total_weight) else return(NA)
}

calculate_TLI <- function(row) {
  chla_val <- as.numeric(row["Chl-a"])
  if(!is.na(chla_val) && chla_val > 0) {
    TLI_Chla <- 10 * (2.5 + 1.086 * log(chla_val))
  } else { TLI_Chla <- NA }
  tn_val <- as.numeric(row["TN"])
  if(!is.na(tn_val) && tn_val > 0) {
    TLI_TN <- 10 * (5.453 + 1.694 * log(tn_val))
  } else { TLI_TN <- NA }
  tp_val <- as.numeric(row["TP"])
  if(!is.na(tp_val) && tp_val > 0) {
    TLI_TP <- 10 * (9.436 + 1.624 * log(tp_val))
  } else { TLI_TP <- NA }
  if(!is.na(TLI_Chla) && !is.na(TLI_TP) && !is.na(TLI_TN)) {
    return(0.4205 * TLI_Chla + 0.2967 * TLI_TP + 0.2828 * TLI_TN)
  } else { return(NA) }
}

# ========== 样本名映射函数 ==========
map_env_to_morph <- function(env_names, tp) {
  if(tp == "DH04") {
    nums <- gsub("DH", "", env_names)
    nums_formatted <- sprintf("%02d", as.numeric(nums))
    return(paste0("C4DH", nums_formatted))
  } 
  else if(tp == "DH07") {
    nums <- gsub("DH", "", env_names)
    nums_formatted <- sprintf("%02d", as.numeric(nums))
    return(paste0("X7DH", nums_formatted))
  } 
  else if(tp == "DH10") {
    nums <- gsub("DH", "", env_names)
    nums_formatted <- sprintf("%02d", as.numeric(nums))
    return(paste0("Q10DH", nums_formatted))
  } 
  else if(tp == "WLSH04") {
    nums <- gsub("WLSH", "", env_names)
    return(paste0("C4WLSH", nums))
  } 
  else if(tp == "WLSH07") {
    nums <- gsub("WLSH", "", env_names)
    return(paste0("X7WLSH", nums))
  } 
  else if(tp == "WLSH10") {
    nums <- gsub("WLSH", "", env_names)
    return(paste0("Q10WLSH", nums))
  }
  return(env_names)
}

read_morph_matrix <- function(morph_file, sample_cols) {
  morph <- read.table(morph_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  morph_matrix <- as.matrix(morph[, sample_cols, drop = FALSE])
  return(morph_matrix)
}

merge_matrices <- function(mat1, mat2) {
  all_rows <- union(rownames(mat1), rownames(mat2))
  new_mat1 <- matrix(0, nrow = length(all_rows), ncol = ncol(mat1))
  new_mat2 <- matrix(0, nrow = length(all_rows), ncol = ncol(mat2))
  rownames(new_mat1) <- all_rows
  rownames(new_mat2) <- all_rows
  colnames(new_mat1) <- colnames(mat1)
  colnames(new_mat2) <- colnames(mat2)
  new_mat1[rownames(mat1), ] <- mat1
  new_mat2[rownames(mat2), ] <- mat2
  return(cbind(new_mat1, new_mat2))
}

# ========== 读取环境因子（保留 WQI/TLI 的完整版，后续再按需筛选）==========
read_env_data_full <- function(env_file) {
  env <- read.table(env_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  site_ids <- as.character(env[, 1])
  if(any(grepl("-", site_ids))) {
    site_ids <- gsub("-", "", site_ids)
  }
  env <- env[, -1, drop = FALSE]
  env <- as.data.frame(lapply(env, as.numeric))
  rownames(env) <- site_ids
  env$WQI <- apply(env, 1, calculate_WQI)
  env$TLI <- apply(env, 1, calculate_TLI)
  return(env)
}

# =====================================================
# 核心函数：按给定环境因子列集合执行 RDA 并绘图
#   env_keep = NULL  -> 使用全部原始环境因子（去掉 WQI/TLI）
#   env_keep = c("WQI","TLI") -> 只用 WQI/TLI
# =====================================================
run_rda_plot <- function(env_keep, out_prefix, title_suffix = "") {
  
  plot_data_list <- list()
  arrow_data_list <- list()
  season_info <- list()
  
  for(season in season_levels) {
    cat("\n===== 正在计算：", season, "(", out_prefix, ") =====\n")
    tps <- names(season_map[season_map == season])
    all_morph <- NULL
    all_env <- NULL
    all_metadata <- data.frame()
    
    for(tp in tps) {
      morph_file <- file.path(data_dir, tp, "morpho.txt")
      env_file <- file.path(data_dir, tp, "env.txt")
      if(!file.exists(morph_file) || !file.exists(env_file)) next
      
      morph_header <- read.table(morph_file, header = TRUE, sep = "\t", check.names = FALSE, nrow = 0)
      morph_samples <- colnames(morph_header)[-1]
      morph_samples <- morph_samples[morph_samples != ""]
      if(length(morph_samples) == 0) next
      
      morph_matrix <- read_morph_matrix(morph_file, morph_samples)
      morph_matrix <- morph_matrix[rowSums(morph_matrix) > 0, , drop = FALSE]
      if(nrow(morph_matrix) < min_species) next
      
      env <- read_env_data_full(env_file)
      
      # ★ 关键：按 env_keep 筛选环境因子
      if(is.null(env_keep)) {
        # 去掉 WQI/TLI
        env <- env[, setdiff(colnames(env), c("WQI","TLI")), drop = FALSE]
      } else {
        keep_cols <- intersect(env_keep, colnames(env))
        if(length(keep_cols) == 0) { next }
        env <- env[, keep_cols, drop = FALSE]
      }
      
      rownames_mapped <- map_env_to_morph(rownames(env), tp)
      common <- intersect(morph_samples, rownames_mapped)
      if(length(common) < 3) next
      
      morph_sub <- morph_matrix[, common, drop = FALSE]
      env_sub <- env[rownames_mapped %in% common, , drop = FALSE]
      rownames(env_sub) <- rownames_mapped[rownames_mapped %in% common]
      
      if(is.null(all_morph)) {
        all_morph <- morph_sub
        all_env <- env_sub
      } else {
        all_morph <- merge_matrices(all_morph, morph_sub)
        all_env <- rbind(all_env, env_sub)
      }
      
      lake <- ifelse(grepl("^DH", tp), "DH", "WLSH")
      for(s in common) {
        all_metadata <- rbind(all_metadata, data.frame(
          Sample = s, Lake = lake, Season = season, stringsAsFactors = FALSE
        ))
      }
    }
    
    if(is.null(all_morph) || ncol(all_morph) < 3) next
    
    env_var <- apply(all_env, 2, var, na.rm = TRUE)
    all_env <- all_env[, env_var > 0 & !is.na(env_var), drop = FALSE]
    if(ncol(all_env) == 0) next
    
    morph_hell <- decostand(t(all_morph), method = "hellinger")
    rda_result <- rda(morph_hell ~ ., data = all_env)
    rda_summ <- summary(rda_result)
    
    R2_adj <- RsquareAdj(rda_result)$adj.r.squared
    if(is.na(R2_adj)) R2_adj <- 0
    r2_value <- round(R2_adj, 3)
    
    rda1_expl <- round(rda_summ$cont$importance[2,1] * 100, 2)
    rda2_expl <- round(rda_summ$cont$importance[2,2] * 100, 2)
    
    season_info[[season]] <- list(
      R2 = r2_value,
      title = paste0(season, " RDA (R² = ", r2_value, ")"),
      explain_text = paste0("RDA1: ", rda1_expl, "%    RDA2: ", rda2_expl, "%")
    )
    
    site_scores <- as.data.frame(scores(rda_result, display = "site", choices = 1:2))
    colnames(site_scores) <- c("RDA1", "RDA2")
    site_scores$Sample <- rownames(site_scores)
    site_scores <- merge(site_scores, all_metadata, by = "Sample", all.x = TRUE)
    site_scores$Season <- season
    
    bp_scores <- as.data.frame(scores(rda_result, display = "bp", choices = 1:2))
    if(nrow(bp_scores) > 0) {
      colnames(bp_scores) <- c("RDA1", "RDA2")
      bp_scores$Env <- rownames(bp_scores)
      bp_scores <- bp_scores[!is.na(bp_scores$RDA1), ]
      bp_scores$Season <- season
    }
    
    plot_data_list[[season]] <- site_scores
    arrow_data_list[[season]] <- bp_scores
  }
  
  successful_seasons <- names(plot_data_list)
  if(length(successful_seasons) == 0) {
    cat("⚠️", out_prefix, "没有成功计算的季节，跳过绘图\n")
    return(NULL)
  }
  
  all_plot_data <- do.call(rbind, plot_data_list)
  all_arrow_data <- do.call(rbind, arrow_data_list)
  
  all_plot_data$Season <- factor(all_plot_data$Season, levels = season_levels)
  if(!is.null(all_arrow_data) && nrow(all_arrow_data) > 0) {
    all_arrow_data$Season <- factor(all_arrow_data$Season, levels = season_levels)
    
    all_arrow_data$RDA1_scaled <- NA
    all_arrow_data$RDA2_scaled <- NA
    for(season in successful_seasons) {
      idx <- which(all_arrow_data$Season == season)
      if(length(idx) > 0) {
        plot_data <- plot_data_list[[season]]
        arrow_data <- arrow_data_list[[season]]
        if(!is.null(arrow_data) && nrow(arrow_data) > 0) {
          x_rng <- range(plot_data$RDA1)
          y_rng <- range(plot_data$RDA2)
          data_scale <- min(diff(x_rng), diff(y_rng)) * 0.8
          max_arrow <- max(abs(arrow_data$RDA1), abs(arrow_data$RDA2))
          arrow_scale <- ifelse(max_arrow > 0, data_scale / max_arrow, 1)
          all_arrow_data$RDA1_scaled[idx] <- arrow_data$RDA1 * arrow_scale
          all_arrow_data$RDA2_scaled[idx] <- arrow_data$RDA2 * arrow_scale
        }
      }
    }
    all_arrow_data <- all_arrow_data[!is.na(all_arrow_data$RDA1_scaled), ]
  }
  
  explain_labels <- data.frame(
    Season = factor(successful_seasons, levels = successful_seasons),
    x = -Inf, y = -Inf,
    label = unlist(lapply(successful_seasons, function(s) season_info[[s]]$explain_text))
  )
  
  p <- ggplot() +
    geom_point(data = all_plot_data,
               aes(x = RDA1, y = RDA2, color = Lake),
               shape = 16, size = 3.5, alpha = 0.8, stroke = 0.8) +
    stat_ellipse(data = all_plot_data,
                 aes(x = RDA1, y = RDA2, color = Lake),
                 level = 0.95, linetype = 1, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_text(data = explain_labels,
              aes(x = x, y = y, label = label),
              size = 3.5, hjust = -0.05, vjust = -0.5,
              color = "black", fontface = "plain") +
    scale_color_manual(values = lake_colors, name = "Lake") +
    facet_wrap(~Season, nrow = 1, scales = "free",
               labeller = as_labeller(function(x) {
                 paste0(x, " RDA (R² = ",
                        sapply(x, function(s) season_info[[s]]$R2), ")")
               })) +
    labs(x = "RDA1", y = "RDA2") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(10, 15, 10, 15)
    )
  
  if(!is.null(all_arrow_data) && nrow(all_arrow_data) > 0) {
    p <- p +
      geom_segment(data = all_arrow_data,
                   aes(x = 0, y = 0, xend = RDA1_scaled, yend = RDA2_scaled),
                   color = "black", linewidth = 0.7,
                   arrow = arrow(length = unit(0.25, "cm"))) +
      geom_text_repel(data = all_arrow_data,
                      aes(x = RDA1_scaled, y = RDA2_scaled, label = Env),
                      size = 3.5, fontface = "bold", color = "darkblue",
                      box.padding = 0.4, point.padding = 0.2,
                      segment.color = "grey50", segment.size = 0.3,
                      max.overlaps = 15)
  }
  
  pdf(file.path(output_dir_morph, paste0(out_prefix, ".pdf")), width = 14, height = 5.5, useDingbats = FALSE)
  print(p)
  dev.off()
  
  ggsave(file.path(output_dir_morph, paste0(out_prefix, ".png")),
         p, width = 14, height = 5.5, dpi = 300, bg = "white")
  
  cat("\n✅ 已保存:", out_prefix, "\n")
  for(season in successful_seasons) {
    cat("  ", season, ": R² =", sprintf("%.3f", season_info[[season]]$R2),
        "|", season_info[[season]]$explain_text, "\n")
  }
  invisible(p)
}

# =====================================================
# 图1：去掉 WQI 和 TLI（只用原始环境因子）
# =====================================================
p_env_only <- run_rda_plot(env_keep = NULL,
                           out_prefix = "Morph_RDA_EnvOnly")

# =====================================================
# 图2：只用 WQI 和 TLI
# =====================================================
p_wqi_tli <- run_rda_plot(env_keep = c("WQI", "TLI"),
                          out_prefix = "Morph_RDA_WQI_TLI")

cat("\n===== 全部完成 =====\n")
cat("输出目录:", output_dir_morph, "\n")
cat("图1（原始环境因子，无 WQI/TLI）: Morph_RDA_EnvOnly.pdf/png\n")
cat("图2（仅 WQI+TLI）: Morph_RDA_WQI_TLI.pdf/png\n")

# =====================================================
# Source: 图4LMN.R
# =====================================================

# =====================================================
# eDNA RDA 分面图 —— 输出两张图
#  图1：只用原始环境因子（去掉 WQI、TLI）
#  图2：只用 WQI 和 TLI
# 两张图绘图设置完全一致（facet_wrap、实线椭圆、左下角解释率）
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ========== 路径 & 基础参数 ==========
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir_edna <- file.path(base_dir, "输出图片/Fig4LMN_分WQI版本")
if(!dir.exists(output_dir_edna)) dir.create(output_dir_edna, recursive = TRUE)

lake_colors <- c("DH" = "#C82423", "WLSH" = "#5A7D9A")

timepoints <- c("DH04", "DH07", "DH10", "WLSH04", "WLSH07", "WLSH10")
season_map <- c("Spring", "Summer", "Autumn", "Spring", "Summer", "Autumn")
names(season_map) <- timepoints
season_levels <- c("Spring", "Summer", "Autumn")

min_rel_abundance <- 0.0001

# ========== WQI 单项评分函数 ==========
get_Ci_TN <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.1) return(100)
  else if(value < 0.2) return(90)
  else if(value < 0.35) return(80)
  else if(value < 0.5) return(70)
  else if(value < 0.75) return(60)
  else if(value < 1) return(50)
  else if(value < 1.25) return(40)
  else if(value < 1.5) return(30)
  else if(value < 1.75) return(20)
  else if(value <= 2) return(10)
  else return(0)
}
get_Ci_TP <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.01) return(100)
  else if(value < 0.02) return(90)
  else if(value < 0.05) return(80)
  else if(value < 0.1) return(70)
  else if(value < 0.15) return(60)
  else if(value < 0.2) return(50)
  else if(value < 0.25) return(40)
  else if(value < 0.3) return(30)
  else if(value < 0.35) return(20)
  else if(value <= 0.4) return(10)
  else return(0)
}
get_Ci_NO3 <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.5) return(100)
  else if(value < 1) return(90)
  else if(value < 1.5) return(80)
  else if(value < 2) return(70)
  else if(value < 3) return(60)
  else if(value < 4) return(50)
  else if(value < 5) return(40)
  else if(value < 6) return(30)
  else if(value < 8) return(20)
  else if(value <= 10) return(10)
  else return(0)
}
get_Ci_NH4 <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 0.01) return(100)
  else if(value < 0.05) return(90)
  else if(value < 0.1) return(80)
  else if(value < 0.2) return(70)
  else if(value < 0.3) return(60)
  else if(value < 0.4) return(50)
  else if(value < 0.5) return(40)
  else if(value < 0.75) return(30)
  else if(value < 1) return(20)
  else if(value <= 1.25) return(10)
  else return(0)
}
get_Ci_Chla <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 1) return(100)
  else if(value < 4) return(90)
  else if(value < 7) return(80)
  else if(value < 10) return(70)
  else if(value < 15) return(60)
  else if(value < 20) return(50)
  else if(value < 30) return(40)
  else if(value < 40) return(30)
  else if(value < 50) return(20)
  else if(value <= 65) return(10)
  else return(0)
}
get_Ci_CODcr <- function(value) {
  if(is.na(value)) return(NA)
  if(value < 15) return(100)
  else if(value < 16) return(90)
  else if(value < 18) return(80)
  else if(value < 19) return(70)
  else if(value < 20) return(60)
  else if(value < 25) return(50)
  else if(value < 30) return(40)
  else if(value < 35) return(30)
  else if(value < 37) return(20)
  else if(value <= 40) return(10)
  else return(0)
}
get_Ci_DO <- function(value) {
  if(is.na(value)) return(NA)
  if(value >= 9) return(100)
  else if(value > 8) return(90)
  else if(value > 7) return(80)
  else if(value > 6) return(70)
  else if(value > 5) return(60)
  else if(value > 4) return(50)
  else if(value > 3.5) return(40)
  else if(value > 3) return(30)
  else if(value > 2) return(20)
  else if(value >= 1) return(10)
  else return(0)
}
get_Ci_pH <- function(value) {
  if(is.na(value)) return(NA)
  if(value == 7) return(100)
  else if(value >= 7 & value <= 8) return(90)
  else if(value >= 7 & value <= 8.5) return(80)
  else if(value >= 6.5 & value <= 9) return(70)
  else if(value >= 6 & value <= 9.5) return(60)
  else if(value >= 5 & value <= 10) return(50)
  else if(value >= 4 & value <= 11) return(40)
  else if(value >= 3 & value <= 12) return(30)
  else if(value >= 2 & value <= 13) return(20)
  else if(value >= 1 & value <= 14) return(10)
  else return(0)
}
get_Ci_Temp <- function(value) {
  if(is.na(value)) return(NA)
  if(value >= 16 & value <= 21) return(100)
  else if((value >= 15 & value < 16) | (value > 21 & value <= 22)) return(90)
  else if((value >= 14 & value < 15) | (value > 22 & value <= 24)) return(80)
  else if((value >= 12 & value < 14) | (value > 24 & value <= 26)) return(70)
  else if((value >= 10 & value < 12) | (value > 26 & value <= 28)) return(60)
  else if((value >= 5 & value < 10) | (value > 28 & value <= 30)) return(50)
  else if((value >= 0 & value < 5) | (value > 30 & value <= 32)) return(40)
  else if((value >= -2 & value < 0) | (value > 32 & value <= 36)) return(30)
  else if((value >= -4 & value < -2) | (value > 36 & value <= 40)) return(20)
  else if((value >= -6 & value < -4) | (value > 40 & value <= 45)) return(10)
  else return(0)
}

# ========== WQI & TLI 综合计算函数（eDNA 列名用 . 形式）==========
calculate_WQI <- function(row) {
  total_weight <- 0; weighted_sum <- 0
  tn_val <- as.numeric(row["TN"])
  if(!is.na(tn_val)) { ci <- get_Ci_TN(tn_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  tp_val <- as.numeric(row["TP"])
  if(!is.na(tp_val)) { ci <- get_Ci_TP(tp_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 1; total_weight <- total_weight + 1 } }
  no3_val <- as.numeric(row["NO3..N"])
  if(!is.na(no3_val)) { ci <- get_Ci_NO3(no3_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 2; total_weight <- total_weight + 2 } }
  nh4_val <- as.numeric(row["NH4..N"])
  if(!is.na(nh4_val)) { ci <- get_Ci_NH4(nh4_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  chla_val <- as.numeric(row["Chl.a"])
  if(!is.na(chla_val)) { ci <- get_Ci_Chla(chla_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  cod_val <- as.numeric(row["CODcr"])
  if(!is.na(cod_val)) { ci <- get_Ci_CODcr(cod_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 3; total_weight <- total_weight + 3 } }
  do_val <- as.numeric(row["DO"])
  if(!is.na(do_val)) { ci <- get_Ci_DO(do_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 4; total_weight <- total_weight + 4 } }
  ph_val <- as.numeric(row["Ph"])
  if(!is.na(ph_val)) { ci <- get_Ci_pH(ph_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 1; total_weight <- total_weight + 1 } }
  temp_val <- as.numeric(row["Temp"])
  if(!is.na(temp_val)) { ci <- get_Ci_Temp(temp_val); if(!is.na(ci)) { weighted_sum <- weighted_sum + ci * 1; total_weight <- total_weight + 1 } }
  if(total_weight > 0) return(weighted_sum / total_weight) else return(NA)
}

calculate_TLI <- function(row) {
  chla_val <- as.numeric(row["Chl.a"])
  if(!is.na(chla_val) && chla_val > 0) {
    TLI_Chla <- 10 * (2.5 + 1.086 * log(chla_val))
  } else { TLI_Chla <- NA }
  tn_val <- as.numeric(row["TN"])
  if(!is.na(tn_val) && tn_val > 0) {
    TLI_TN <- 10 * (5.453 + 1.694 * log(tn_val))
  } else { TLI_TN <- NA }
  tp_val <- as.numeric(row["TP"])
  if(!is.na(tp_val) && tp_val > 0) {
    TLI_TP <- 10 * (9.436 + 1.624 * log(tp_val))
  } else { TLI_TP <- NA }
  if(!is.na(TLI_Chla) && !is.na(TLI_TP) && !is.na(TLI_TN)) {
    return(0.4205 * TLI_Chla + 0.2967 * TLI_TP + 0.2828 * TLI_TN)
  } else { return(NA) }
}

# ========== 数据读取与预处理工具函数 ==========
map_env_to_otu <- function(env_names, tp) {
  nums <- gsub("DH", "", env_names)
  nums <- gsub("WLSH", "", nums)
  nums_int <- as.numeric(nums)
  nums_formatted <- sprintf("%02d", nums_int)
  if(tp == "DH04") return(paste0("C4DH", nums_formatted))
  else if(tp == "DH07") return(paste0("X7DH", nums_formatted))
  else if(tp == "DH10") return(paste0("Q10DH", nums_formatted))
  else if(tp == "WLSH04") return(paste0("C4WLSH", nums))
  else if(tp == "WLSH07") return(paste0("X7WLSH", nums))
  else if(tp == "WLSH10") return(paste0("Q10WLSH", nums))
  return(env_names)
}

read_asv_matrix <- function(asv_file, sample_cols) {
  asv <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE)
  asv_matrix <- as.matrix(asv[, sample_cols, drop = FALSE])
  rownames(asv_matrix) <- asv$ASV_ID
  return(asv_matrix)
}

filter_low_abundance <- function(asv_matrix, min_rel_abund = 0.0001) {
  if(min_rel_abund <= 0) return(asv_matrix)
  total_sum <- sum(asv_matrix)
  rel_abund <- rowSums(asv_matrix) / total_sum
  keep <- rel_abund >= min_rel_abund
  return(asv_matrix[keep, , drop = FALSE])
}

merge_asv_matrices <- function(mat1, mat2) {
  all_asvs <- union(rownames(mat1), rownames(mat2))
  new_mat1 <- matrix(0, nrow = length(all_asvs), ncol = ncol(mat1))
  new_mat2 <- matrix(0, nrow = length(all_asvs), ncol = ncol(mat2))
  rownames(new_mat1) <- all_asvs
  rownames(new_mat2) <- all_asvs
  colnames(new_mat1) <- colnames(mat1)
  colnames(new_mat2) <- colnames(mat2)
  new_mat1[rownames(mat1), ] <- mat1
  new_mat2[rownames(mat2), ] <- mat2
  return(cbind(new_mat1, new_mat2))
}

# 读取环境因子（保留 WQI/TLI 的完整版，后续按需筛选）
read_env_data_full <- function(env_file) {
  env <- read.table(env_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  site_ids <- env$SiteID
  env <- env[, -1, drop = FALSE]
  env <- as.data.frame(lapply(env, as.numeric))
  rownames(env) <- site_ids
  env$WQI <- apply(env, 1, calculate_WQI)
  env$TLI <- apply(env, 1, calculate_TLI)
  return(env)
}

# =====================================================
# 核心函数：按给定环境因子列集合执行 RDA 并绘图
#   env_keep = NULL            -> 只用原始环境因子（去掉 WQI/TLI）
#   env_keep = c("WQI","TLI")  -> 只用 WQI/TLI
# =====================================================
run_rda_plot <- function(env_keep, out_prefix) {
  
  plot_data_list <- list()
  arrow_data_list <- list()
  season_info <- list()
  
  for(season in season_levels) {
    cat("\n===== 正在计算：", season, "(", out_prefix, ") =====\n")
    tps <- names(season_map[season_map == season])
    all_asv <- NULL
    all_env <- NULL
    all_metadata <- data.frame()
    
    for(tp in tps) {
      asv_file <- file.path(data_dir, tp, "ASV.txt")
      env_file <- file.path(data_dir, tp, "env_asv.txt")
      if(!file.exists(asv_file) || !file.exists(env_file)) next
      
      asv_header <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE, nrow = 0)
      sample_cols <- setdiff(colnames(asv_header), c("ASV_ID", "taxonomy"))
      if(length(sample_cols) == 0) next
      
      asv_matrix <- read_asv_matrix(asv_file, sample_cols)
      asv_matrix <- filter_low_abundance(asv_matrix, min_rel_abundance)
      if(nrow(asv_matrix) < 3) next
      
      env <- read_env_data_full(env_file)
      
      # ★ 按 env_keep 筛选环境因子
      if(is.null(env_keep)) {
        env <- env[, setdiff(colnames(env), c("WQI","TLI")), drop = FALSE]
      } else {
        keep_cols <- intersect(env_keep, colnames(env))
        if(length(keep_cols) == 0) next
        env <- env[, keep_cols, drop = FALSE]
      }
      
      rownames(env) <- map_env_to_otu(rownames(env), tp)
      common <- intersect(colnames(asv_matrix), rownames(env))
      if(length(common) < 3) next
      
      asv_sub <- asv_matrix[, common, drop = FALSE]
      env_sub <- env[common, , drop = FALSE]
      asv_sub <- asv_sub[rowSums(asv_sub) > 0, , drop = FALSE]
      if(nrow(asv_sub) < 3) next
      
      if(is.null(all_asv)) {
        all_asv <- asv_sub
        all_env <- env_sub
      } else {
        all_asv <- merge_asv_matrices(all_asv, asv_sub)
        all_env <- rbind(all_env, env_sub)
      }
      
      lake <- ifelse(grepl("^DH", tp), "DH", "WLSH")
      for(s in common) {
        all_metadata <- rbind(all_metadata, data.frame(
          Sample = s, Lake = lake, Season = season, stringsAsFactors = FALSE
        ))
      }
    }
    
    if(is.null(all_asv) || ncol(all_asv) < 3) next
    
    env_var <- apply(all_env, 2, var, na.rm = TRUE)
    all_env <- all_env[, env_var > 0 & !is.na(env_var), drop = FALSE]
    if(ncol(all_env) == 0) next
    
    asv_hell <- decostand(t(all_asv), method = "hellinger")
    rda_result <- rda(asv_hell ~ ., data = all_env)
    rda_summ <- summary(rda_result)
    
    R2_adj <- RsquareAdj(rda_result)$adj.r.squared
    if(is.na(R2_adj)) R2_adj <- 0
    r2_value <- round(R2_adj, 3)
    
    rda1_expl <- round(rda_summ$cont$importance[2,1] * 100, 2)
    rda2_expl <- round(rda_summ$cont$importance[2,2] * 100, 2)
    
    facet_title <- paste0(season, " RDA (R² = ", r2_value, ")")
    
    season_info[[season]] <- list(
      R2 = r2_value,
      rda1 = rda1_expl,
      rda2 = rda2_expl,
      title = facet_title,
      explain_text = paste0("RDA1: ", rda1_expl, "%    RDA2: ", rda2_expl, "%")
    )
    
    site_scores <- as.data.frame(scores(rda_result, display = "site", choices = 1:2))
    colnames(site_scores) <- c("RDA1", "RDA2")
    site_scores$Sample <- rownames(site_scores)
    site_scores <- merge(site_scores, all_metadata, by = "Sample", all.x = TRUE)
    site_scores$Season <- season
    
    bp_scores <- as.data.frame(scores(rda_result, display = "bp", choices = 1:2))
    if(nrow(bp_scores) > 0) {
      colnames(bp_scores) <- c("RDA1", "RDA2")
      bp_scores$Env <- rownames(bp_scores)
      bp_scores <- bp_scores[!is.na(bp_scores$RDA1), ]
      bp_scores$Season <- season
    }
    
    plot_data_list[[season]] <- site_scores
    arrow_data_list[[season]] <- bp_scores
    
    write.csv(site_scores, file.path(output_dir_edna, paste0(out_prefix, "_", season, "_Coordinates.csv")), row.names = FALSE)
    if(nrow(bp_scores) > 0) write.csv(bp_scores, file.path(output_dir_edna, paste0(out_prefix, "_", season, "_Arrows.csv")), row.names = FALSE)
    saveRDS(rda_result, file.path(output_dir_edna, paste0(out_prefix, "_", season, "_result.rds")))
  }
  
  successful_seasons <- names(plot_data_list)
  if(length(successful_seasons) == 0) {
    cat("⚠️", out_prefix, "没有成功计算的季节，跳过绘图\n")
    return(NULL)
  }
  
  all_plot_data <- do.call(rbind, plot_data_list)
  all_arrow_data <- do.call(rbind, arrow_data_list)
  
  all_plot_data$Season <- factor(all_plot_data$Season, levels = season_levels)
  if(!is.null(all_arrow_data) && nrow(all_arrow_data) > 0) {
    all_arrow_data$Season <- factor(all_arrow_data$Season, levels = season_levels)
    all_arrow_data$RDA1_scaled <- NA
    all_arrow_data$RDA2_scaled <- NA
    for(season in successful_seasons) {
      idx <- which(all_arrow_data$Season == season)
      if(length(idx) > 0) {
        plot_data <- plot_data_list[[season]]
        arrow_data <- arrow_data_list[[season]]
        if(!is.null(arrow_data) && nrow(arrow_data) > 0) {
          x_rng <- range(plot_data$RDA1)
          y_rng <- range(plot_data$RDA2)
          data_scale <- min(diff(x_rng), diff(y_rng)) * 0.8
          max_arrow <- max(abs(arrow_data$RDA1), abs(arrow_data$RDA2))
          arrow_scale <- ifelse(max_arrow > 0, data_scale / max_arrow, 1)
          all_arrow_data$RDA1_scaled[idx] <- arrow_data$RDA1 * arrow_scale
          all_arrow_data$RDA2_scaled[idx] <- arrow_data$RDA2 * arrow_scale
        }
      }
    }
    all_arrow_data <- all_arrow_data[!is.na(all_arrow_data$RDA1_scaled), ]
  }
  
  facet_labels <- sapply(successful_seasons, function(s) season_info[[s]]$title)
  
  explain_labels <- data.frame(
    Season = factor(successful_seasons, levels = successful_seasons),
    x = -Inf, y = -Inf,
    label = sapply(successful_seasons, function(s) season_info[[s]]$explain_text)
  )
  
  p <- ggplot() +
    geom_point(data = all_plot_data,
               aes(x = RDA1, y = RDA2, color = Lake),
               shape = 16, size = 3.5, alpha = 0.8) +
    stat_ellipse(data = all_plot_data,
                 aes(x = RDA1, y = RDA2, color = Lake),
                 level = 0.95, linetype = 1, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_text(data = explain_labels,
              aes(x = x, y = y, label = label),
              size = 3.5, hjust = -0.05, vjust = -0.5,
              color = "black", fontface = "plain") +
    scale_color_manual(values = lake_colors, name = "Lake") +
    facet_wrap(~Season, nrow = 1, scales = "free",
               labeller = labeller(Season = facet_labels)) +
    labs(x = "RDA1", y = "RDA2") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(10, 15, 10, 15)
    )
  
  if(!is.null(all_arrow_data) && nrow(all_arrow_data) > 0) {
    p <- p +
      geom_segment(data = all_arrow_data,
                   aes(x = 0, y = 0, xend = RDA1_scaled, yend = RDA2_scaled),
                   color = "black", linewidth = 0.7,
                   arrow = arrow(length = unit(0.25, "cm"))) +
      geom_text_repel(data = all_arrow_data,
                      aes(x = RDA1_scaled, y = RDA2_scaled, label = Env),
                      size = 3.5, fontface = "bold", color = "darkblue",
                      box.padding = 0.4, point.padding = 0.2,
                      segment.color = "grey50", segment.size = 0.3,
                      max.overlaps = 15)
  }
  
  pdf(file.path(output_dir_edna, paste0(out_prefix, ".pdf")), width = 14, height = 5.5, useDingbats = FALSE)
  print(p)
  dev.off()
  
  ggsave(file.path(output_dir_edna, paste0(out_prefix, ".png")),
         p, width = 14, height = 5.5, dpi = 300, bg = "white")
  
  cat("\n✅ 已保存:", out_prefix, "\n")
  for(season in successful_seasons) {
    cat("  ", season, ": R² =", sprintf("%.3f", season_info[[season]]$R2),
        "|", season_info[[season]]$explain_text, "\n")
  }
  invisible(p)
}

# =====================================================
# 图1：去掉 WQI 和 TLI（只用原始环境因子）
# =====================================================
p_env_only <- run_rda_plot(env_keep = NULL,
                           out_prefix = "eDNA_RDA_EnvOnly")

# =====================================================
# 图2：只用 WQI 和 TLI
# =====================================================
p_wqi_tli <- run_rda_plot(env_keep = c("WQI", "TLI"),
                          out_prefix = "eDNA_RDA_WQI_TLI")

cat("\n===== 全部完成 =====\n")
cat("输出目录:", output_dir_edna, "\n")
cat("图1（原始环境因子，无 WQI/TLI）: eDNA_RDA_EnvOnly.pdf/png\n")
cat("图2（仅 WQI+TLI）: eDNA_RDA_WQI_TLI.pdf/png\n")
