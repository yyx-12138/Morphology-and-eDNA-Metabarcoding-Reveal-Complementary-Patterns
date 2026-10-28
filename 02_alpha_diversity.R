# 02 Alpha diversity


# =====================================================
# Source: 图4A_Morph.R
# =====================================================

# =====================================================
# 形态学数据：多样性分析 + 出图（不分季节，分面图）
# 按照论文公式：
#    Shannon: H = -∑P_i × log₂(P_i)
#    Margalef: D = (S-1) / log₂(N)
#    Pielou: J = H / ln(S)
# 图结构：X轴为湖泊（DH, WLSH），分面为三个指数（1行3列）
# 每个分面独立计算显著性，纵坐标独立
# 保存文件夹：Fig4B
# 分面顺序：Shannon, Pielou, Margalef（强制从左到右）
# 
# 重要：每个样本使用自己的S（该样本的物种数）和N（该样本的总个体数）
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/Fig4A_Morph")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 配色
lake_colors <- c("DH" = "#D6857E", "WLSH" = "#5A7D9A")

# =====================================================
# 按照论文公式定义函数（每个样本用自己的S和N）
# =====================================================
# Shannon指数: log2
shannon_log2 <- function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(NA)
  p <- x / sum(x)
  -sum(p * log2(p))
}

# Pielou指数: H / ln(S) - H用log2计算的Shannon，分母用ln(S)
pielou_formula <- function(x) {
  x <- x[x > 0]
  if(length(x) <= 1) return(NA)
  H <- shannon_log2(x)  # 使用log2的Shannon
  S <- length(x)
  H / log(S)  # 自然对数 ln(S)
}

# Margalef指数: (S-1) / log2(N)
margalef_log2 <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  N <- sum(x)
  if(N > 0) (S - 1) / log2(N) else NA
}
# =====================================================
# 第一步：修复所有形态学文件（合并重复物种）
# =====================================================

cat("========== 第一步：修复数据文件 ==========\n")

timepoints <- c("DH04", "DH07", "DH10", "WLSH04", "WLSH07", "WLSH10")

fix_and_merge <- function(file_path) {
  morpho <- read.table(file_path, 
                       header = TRUE, 
                       sep = "\t", 
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       fill = TRUE)
  
  tp_name <- basename(dirname(file_path))
  cat("\n处理:", tp_name, "\n")
  cat("  原始行数:", nrow(morpho), "\n")
  
  species_names <- morpho[[1]]
  
  # 移除空物种名的行
  empty_rows <- which(species_names == "" | is.na(species_names))
  if(length(empty_rows) > 0) {
    cat("  移除空物种名行:", length(empty_rows), "行\n")
    morpho <- morpho[-empty_rows, ]
    species_names <- morpho[[1]]
  }
  
  # 合并重复的物种名
  dup_names <- unique(species_names[duplicated(species_names)])
  
  if(length(dup_names) > 0) {
    cat("  发现重复物种名:", paste(dup_names, collapse=", "), "\n")
    
    for(dup in dup_names) {
      dup_rows <- which(species_names == dup)
      
      if(length(dup_rows) > 1) {
        sample_cols <- 2:ncol(morpho)
        
        for(i in 2:length(dup_rows)) {
          morpho[dup_rows[1], sample_cols] <- as.numeric(morpho[dup_rows[1], sample_cols]) + 
            as.numeric(morpho[dup_rows[i], sample_cols])
        }
        
        morpho <- morpho[-dup_rows[-1], ]
        cat("    合并:", dup, "- 共", length(dup_rows), "行\n")
        species_names <- morpho[[1]]
      }
    }
  }
  
  write.table(morpho, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  保存成功，新行数:", nrow(morpho), "\n")
  
  return(morpho)
}

# 修复所有文件
for(tp in timepoints) {
  file_path <- file.path(data_dir, tp, "morpho.txt")
  if(file.exists(file_path)) {
    fix_and_merge(file_path)
  } else {
    cat("⚠️ 文件不存在:", file_path, "\n")
  }
}

# =====================================================
# 第二步：读取修复后的数据并计算多样性（每个样本用自己的S和N）
# =====================================================

cat("\n========== 第二步：计算多样性指数（每个样本用自己的S和N）==========\n")

get_lake <- function(tp) {
  if(grepl("^DH", tp)) return("DH")
  if(grepl("^WLSH", tp)) return("WLSH")
  return(NA)
}

all_diversity <- data.frame()

for(tp in timepoints) {
  file_path <- file.path(data_dir, tp, "morpho.txt")
  
  if(!file.exists(file_path)) {
    cat("⚠️ 缺失:", file_path, "\n")
    next
  }
  
  morpho <- read.table(file_path, 
                       header = TRUE, 
                       sep = "\t", 
                       check.names = FALSE,
                       row.names = 1,
                       stringsAsFactors = FALSE)
  
  morpho <- as.data.frame(lapply(morpho, function(x) as.numeric(as.character(x))))
  morpho[is.na(morpho)] <- 0
  
  lake <- get_lake(tp)
  
  cat("✓", tp, "- 湖泊:", lake, "- 样本数:", ncol(morpho), "\n")
  
  for(sample_name in colnames(morpho)) {
    abu <- morpho[[sample_name]]
    abu <- abu[abu > 0]
    
    # 每个样本用自己的S和N计算
    sha <- shannon_log2(abu)
    pie <- pielou_formula(abu)
    mar <- margalef_log2(abu)
    
    all_diversity <- rbind(all_diversity, data.frame(
      Sample = sample_name,
      Lake = lake,
      Shannon = sha,
      Pielou = pie,
      Margalef = mar,
      stringsAsFactors = FALSE
    ))
  }
}

# 检查数据
cat("\n数据统计:\n")
print(table(all_diversity$Lake))

if(nrow(all_diversity) == 0) {
  stop("没有找到任何数据文件！")
}

# =====================================================
# 第三步：转为长格式（强制顺序：Shannon, Pielou, Margalef）
# =====================================================

all_long <- all_diversity %>%
  pivot_longer(cols = c(Shannon, Pielou, Margalef), 
               names_to = "Index", 
               values_to = "Value") %>%
  filter(!is.na(Value))

all_long$Index <- factor(all_long$Index, 
                         levels = c("Shannon", "Pielou", "Margalef"))
all_long$Lake <- factor(all_long$Lake, levels = c("DH", "WLSH"))

# 验证顺序
cat("\n因子顺序验证:\n")
print(levels(all_long$Index))

# =====================================================
# 第四步：计算显著性（每个指数独立计算，不分季节）
# =====================================================

cat("\n========== 第三步：显著性检验（每个指数独立）==========\n")

# 为每个指数计算显著性标记的位置
sig_data <- data.frame()

for(idx in c("Shannon", "Pielou", "Margalef")) {
  
  # 提取当前指数的数据
  sub_data <- all_long %>% filter(Index == idx)
  
  # 提取两个湖泊的数据
  dh_values <- sub_data %>% filter(Lake == "DH") %>% pull(Value)
  wl_values <- sub_data %>% filter(Lake == "WLSH") %>% pull(Value)
  
  # 计算显著性
  if(length(dh_values) >= 2 & length(wl_values) >= 2) {
    pv <- wilcox.test(dh_values, wl_values)$p.value
    
    # 确定显著性标记
    if(pv < 0.001) {
      sig_label <- "***"
    } else if(pv < 0.01) {
      sig_label <- "**"
    } else if(pv < 0.05) {
      sig_label <- "*"
    } else {
      sig_label <- "ns"
    }
    
    # 计算Y轴位置（在当前指数的最大值上方）
    y_max <- max(sub_data$Value, na.rm = TRUE)
    y_min <- min(sub_data$Value, na.rm = TRUE)
    y_pos <- y_max + (y_max - y_min) * 0.15
    
    # 存储显著性数据
    sig_data <- rbind(sig_data, data.frame(
      Index = idx,
      x = 1.5,
      y = y_pos,
      label = sig_label,
      p_value = pv,
      y_max_orig = y_max,
      y_min_orig = y_min,
      y_max_new = y_max + (y_max - y_min) * 0.25
    ))
    
    cat("  ", idx, ": p =", format(pv, scientific = TRUE), sig_label, "\n")
  }
}

if(nrow(sig_data) == 0) {
  cat("  无显著差异 (p < 0.05)\n")
}

# 计算线段位置
sig_data <- sig_data %>%
  mutate(
    y_seg_end = y - (y_max_new - y_min_orig) * 0.05,
    y_text = y
  )

# =====================================================
# 第五步：创建分面图（1行3列，横版，强制顺序，纵坐标自由）
# =====================================================

cat("\n========== 第四步：绘图 ==========\n")

# 创建标签
index_labels <- c(
  Shannon = "Shannon",
  Pielou = "Pielou", 
  Margalef = "Margalef"
)

# 主图
p <- ggplot(all_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 21, size = 2.5, alpha = 0.7, width = 0.15, stroke = 0.3) +
  facet_wrap(~Index, 
             nrow = 1, 
             ncol = 3,
             scales = "free_y",
             labeller = labeller(Index = as_labeller(index_labels))) +
  scale_fill_manual(values = lake_colors, name = "Lake") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11, margin = margin(8, 0, 8, 0)),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 15)
  )

# 添加显著性标记
if(nrow(sig_data) > 0) {
  sig_plot_data <- sig_data %>% filter(label != "ns")
  
  if(nrow(sig_plot_data) > 0) {
    p <- p +
      geom_segment(data = sig_plot_data,
                   aes(x = 1, xend = 2, 
                       y = y_seg_end, yend = y_seg_end),
                   linewidth = 0.6, inherit.aes = FALSE) +
      geom_text(data = sig_plot_data,
                aes(x = x, y = y_text, label = label),
                size = 6, fontface = "bold", inherit.aes = FALSE)
  }
}

# =====================================================
# 第六步：保存为可编辑PDF和预览图
# =====================================================

cat("\n========== 第五步：保存文件 ==========\n")

# 保存PDF
pdf(file.path(output_dir, "Fig4B_Morphology_Diversity_Editable.pdf"), 
    width = 8, height = 5, 
    useDingbats = FALSE)
print(p)
dev.off()

# 保存PNG用于预览
ggsave(file.path(output_dir, "Fig4B_Morphology_Diversity_Preview.png"), 
       p, width = 8, height = 5, dpi = 300, bg = "white")

# 保存每个指数的单独图
for(idx in c("Shannon", "Pielou", "Margalef")) {
  p_single <- all_long %>%
    filter(Index == idx) %>%
    ggplot(aes(x = Lake, y = Value, fill = Lake)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(shape = 21, size = 2.5, alpha = 0.7, width = 0.15, stroke = 0.3) +
    scale_fill_manual(values = lake_colors) +
    labs(x = NULL, y = NULL, title = idx) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
    )
  
  # 添加显著性标记
  sig_sub <- sig_data %>% filter(Index == idx, label != "ns")
  if(nrow(sig_sub) > 0) {
    y_range <- range(all_long %>% filter(Index == idx) %>% pull(Value), na.rm = TRUE)
    y_offset <- (y_range[2] - y_range[1]) * 0.05
    
    p_single <- p_single +
      annotate("segment", x = 1, xend = 2, 
               y = y_range[2] + y_offset, yend = y_range[2] + y_offset,
               linewidth = 0.6) +
      annotate("text", x = 1.5, y = y_range[2] + y_offset * 1.5, 
               label = sig_sub$label[1], size = 6, fontface = "bold") +
      coord_cartesian(ylim = c(y_range[1], y_range[2] + (y_range[2] - y_range[1]) * 0.2))
  }
  
  ggsave(file.path(output_dir, paste0("Fig4B_Morphology_", idx, ".pdf")), 
         p_single, width = 2.5, height = 4, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0("Fig4B_Morphology_", idx, ".png")), 
         p_single, width = 2.5, height = 4, dpi = 300, bg = "white")
}

# =====================================================
# 完成
# =====================================================

cat("\n✅ 完成！生成横版分面图（1行3列，不分季节）：\n")
cat("   - 图片尺寸：宽度8英寸，高度5英寸\n")
cat("   - 分面顺序：Shannon, Pielou, Margalef（从左到右）\n")
cat("   - 纵坐标：每个分面独立（scales = free_y）\n")
cat("   - 每个样本使用自己的S和N计算多样性\n")
cat("   - 箱线图 + 抖动点\n")
cat("   - 每个指数独立计算显著性（*, **, ***）\n\n")
cat("📁 输出文件夹: Fig4B\n")
cat("   输出文件:\n")
cat("   - Fig4B_Morphology_Diversity_Editable.pdf\n")
cat("   - Fig4B_Morphology_Diversity_Preview.png\n")
cat("   - Fig4B_Morphology_Shannon.pdf/png\n")
cat("   - Fig4B_Morphology_Pielou.pdf/png\n")
cat("   - Fig4B_Morphology_Margalef.pdf/png\n\n")
cat("按照论文公式计算（每个样本用自己的S和N）：\n")
cat("  公式(2-3) Shannon: H = -∑P_i × log₂(P_i)\n")
cat("  公式(2-4) Margalef: D = (S-1) / log₂(N)\n")
cat("  公式(2-5) Pielou: J = H / ln(S)\n")

# =====================================================
# Source: 图4B_eDNA.R
# =====================================================

# =====================================================
# eDNA多样性分析 - 分面图（不分季节）- 属水平
# 排除Unclassified和NA
# 每个样本用自己的S
# 保存文件夹：Fig4B.eDNA
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/Fig4B.eDNA")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 配色
lake_colors <- c("DH" = "#D6857E", "WLSH" = "#5A7D9A")

# =====================================================
# 按照论文公式定义函数
# =====================================================

# Shannon指数: 使用 log2
shannon_log2 <- function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(NA)
  p <- x / sum(x)
  -sum(p * log2(p))
}

# Pielou指数: H / ln(S) - 用样本自己的S
pielou_formula <- function(x) {
  x <- x[x > 0]
  if(length(x) <= 1) return(NA)
  H <- shannon_log2(x)
  S <- length(x)
  H / log(S)
}

# Margalef指数: (S-1) / log₂(N) - 用样本自己的S和N
margalef_log2 <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  N <- sum(x)
  if(N > 0) (S - 1) / log2(N) else NA
}

# =====================================================
# 计算多样性指数（属水平，排除Unclassified和NA）
# =====================================================

timepoints <- c("DH04","DH07","DH10","WLSH04","WLSH07","WLSH10")

all_diversity <- data.frame()

for(tp in timepoints) {
  asv_file <- file.path(data_dir, tp, "ASV.txt")
  tax_file <- file.path(data_dir, tp, "tax.txt")
  
  if(!file.exists(asv_file)) {
    cat("⚠️ ASV文件不存在:", asv_file, "\n")
    next
  }
  if(!file.exists(tax_file)) {
    cat("⚠️ tax文件不存在:", tax_file, "\n")
    next
  }
  
  # 读取ASV表
  asv <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE)
  sample_cols <- setdiff(colnames(asv), c("ASV_ID", "taxonomy"))
  if(length(sample_cols) == 0) next
  
  # 读取分类文件
  tax <- read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 提取属水平，排除Unclassified和NA
  tax$Genus <- gsub("^g__", "", tax$Genus)
  # 只保留有明确属名的ASV（排除NA、空值、Unclassified、uncultured）
  tax_clean <- tax[!is.na(tax$Genus) & 
                     tax$Genus != "" & 
                     tax$Genus != "Unclassified" & 
                     tax$Genus != "uncultured", ]
  
  if(nrow(tax_clean) == 0) {
    cat("⚠️", tp, "没有有效属水平分类\n")
    next
  }
  
  # 创建ASV到属的映射
  asv_to_genus <- tax_clean$Genus
  names(asv_to_genus) <- tax_clean$ASV_ID
  
  # 提取ASV丰度矩阵
  asv_matrix <- as.matrix(asv[, sample_cols])
  rownames(asv_matrix) <- asv$ASV_ID
  
  # 只保留有分类信息的ASV
  valid_asvs <- intersect(rownames(asv_matrix), names(asv_to_genus))
  asv_matrix <- asv_matrix[valid_asvs, , drop = FALSE]
  
  if(nrow(asv_matrix) == 0) next
  
  # 获取所有唯一的属名
  unique_genus <- unique(asv_to_genus[valid_asvs])
  
  # 创建属水平矩阵（同一属的ASV合并）
  genus_matrix <- matrix(0, nrow = ncol(asv_matrix), ncol = length(unique_genus))
  rownames(genus_matrix) <- colnames(asv_matrix)
  colnames(genus_matrix) <- unique_genus
  
  # 将同一属的ASV丰度相加（聚合）
  for(asv_id in rownames(asv_matrix)) {
    genus <- asv_to_genus[asv_id]
    genus_matrix[, genus] <- genus_matrix[, genus] + asv_matrix[asv_id, ]
  }
  
  lake <- ifelse(grepl("DH", tp), "DH", "WLSH")
  
  cat("处理:", tp, "-", lake, "- 样本数:", nrow(genus_matrix), "- 有效属数:", ncol(genus_matrix), "\n")
  
  # 计算每个样本的多样性
  for(sample_name in rownames(genus_matrix)) {
    abu <- genus_matrix[sample_name, ]
    
    sha <- shannon_log2(abu)
    pie <- pielou_formula(abu)
    mar <- margalef_log2(abu)
    
    all_diversity <- rbind(all_diversity, data.frame(
      Sample = sample_name,
      Lake = lake,
      Shannon = sha,
      Pielou = pie,
      Margalef = mar,
      stringsAsFactors = FALSE
    ))
  }
}

# 查看结果统计
cat("\n多样性指数统计（排除Unclassified，属水平）:\n")
cat("Shannon: H = -∑P_i × log₂(P_i)\n")
cat("Margalef: D = (S-1) / log₂(N)\n")
cat("Pielou: J = H / ln(S)\n")
cat("注：已排除Unclassified和NA的ASV\n\n")

# 转换成长格式，便于分面
all_long <- all_diversity %>%
  pivot_longer(cols = c(Shannon, Pielou, Margalef), 
               names_to = "Index", 
               values_to = "Value") %>%
  filter(!is.na(Value))

# 设置因子顺序
all_long$Lake <- factor(all_long$Lake, levels = c("DH", "WLSH"))
all_long$Index <- factor(all_long$Index, levels = c("Shannon", "Pielou", "Margalef"))

# =====================================================
# 计算显著性（每个指数独立计算）
# =====================================================

sig_data <- data.frame()

for(idx in c("Shannon", "Pielou", "Margalef")) {
  
  sub_data <- all_long %>% filter(Index == idx)
  
  dh_values <- sub_data %>% filter(Lake == "DH") %>% pull(Value)
  wl_values <- sub_data %>% filter(Lake == "WLSH") %>% pull(Value)
  
  if(length(dh_values) >= 2 & length(wl_values) >= 2) {
    pv <- wilcox.test(dh_values, wl_values)$p.value
    
    if(pv < 0.001) {
      sig_label <- "***"
    } else if(pv < 0.01) {
      sig_label <- "**"
    } else if(pv < 0.05) {
      sig_label <- "*"
    } else {
      sig_label <- "ns"
    }
    
    y_max <- max(sub_data$Value, na.rm = TRUE)
    y_min <- min(sub_data$Value, na.rm = TRUE)
    y_pos <- y_max + (y_max - y_min) * 0.15
    
    sig_data <- rbind(sig_data, data.frame(
      Index = idx,
      x = 1.5,
      y = y_pos,
      label = sig_label,
      p_value = pv,
      y_max_orig = y_max,
      y_min_orig = y_min,
      y_max_new = y_max + (y_max - y_min) * 0.25
    ))
  }
}

cat("\n显著性检验结果（Wilcoxon检验，排除Unclassified）:\n")
for(i in 1:nrow(sig_data)) {
  cat(sig_data$Index[i], ": p =", sig_data$p_value[i], sig_data$label[i], "\n")
}

# 计算线段位置
sig_data <- sig_data %>%
  mutate(
    y_seg_end = y - (y_max_new - y_min_orig) * 0.05,
    y_text = y
  )

# =====================================================
# 创建分面图（1行3列，横版）
# =====================================================

index_labels <- c(
  Shannon = "Shannon",
  Pielou = "Pielou", 
  Margalef = "Margalef"
)

p <- ggplot(all_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 21, size = 2, alpha = 0.7, width = 0.15, stroke = 0.3) +
  facet_wrap(~Index, nrow = 1, ncol = 3, 
             scales = "free_y",
             labeller = labeller(Index = as_labeller(index_labels))) +
  scale_fill_manual(values = lake_colors, name = "Lake") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 10, margin = margin(6, 0, 6, 0)),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(8, 8, 8, 8)
  )

# 添加显著性标记
if(nrow(sig_data) > 0) {
  sig_plot_data <- sig_data %>% filter(label != "ns")
  
  if(nrow(sig_plot_data) > 0) {
    p <- p +
      geom_segment(data = sig_plot_data,
                   aes(x = 1, xend = 2, 
                       y = y_seg_end, yend = y_seg_end),
                   linewidth = 0.5, inherit.aes = FALSE) +
      geom_text(data = sig_plot_data,
                aes(x = x, y = y_text, label = label),
                size = 5, fontface = "bold", inherit.aes = FALSE)
  }
}

# =====================================================
# 保存为可编辑PDF和预览图
# =====================================================

# 保存PDF
pdf(file.path(output_dir, "Fig4B_eDNA_Diversity_Editable.pdf"), 
    width = 8, height = 5, 
    useDingbats = FALSE)
print(p)
dev.off()

# 保存PNG预览
ggsave(file.path(output_dir, "Fig4B_eDNA_Diversity_Preview.png"), 
       p, width = 8, height = 5, dpi = 300, bg = "white")

# 保存每个指数的单独图
for(idx in c("Shannon", "Pielou", "Margalef")) {
  p_single <- all_long %>%
    filter(Index == idx) %>%
    ggplot(aes(x = Lake, y = Value, fill = Lake)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(shape = 21, size = 2, alpha = 0.7, width = 0.15, stroke = 0.3) +
    scale_fill_manual(values = lake_colors) +
    labs(x = NULL, y = NULL, title = idx) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
    )
  
  # 添加显著性标记
  sig_sub <- sig_data %>% filter(Index == idx, label != "ns")
  if(nrow(sig_sub) > 0) {
    y_range <- range(all_long %>% filter(Index == idx) %>% pull(Value), na.rm = TRUE)
    y_offset <- (y_range[2] - y_range[1]) * 0.05
    
    p_single <- p_single +
      annotate("segment", x = 1, xend = 2, 
               y = y_range[2] + y_offset, yend = y_range[2] + y_offset,
               linewidth = 0.5) +
      annotate("text", x = 1.5, y = y_range[2] + y_offset * 1.5, 
               label = sig_sub$label[1], size = 5, fontface = "bold") +
      coord_cartesian(ylim = c(y_range[1], y_range[2] + (y_range[2] - y_range[1]) * 0.2))
  }
  
  ggsave(file.path(output_dir, paste0("Fig4B_eDNA_", idx, ".pdf")), 
         p_single, width = 2.5, height = 4, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0("Fig4B_eDNA_", idx, ".png")), 
         p_single, width = 2.5, height = 4, dpi = 300, bg = "white")
}

cat("\n✅ 完成！生成横版分面图（1行3列，不分季节，属水平）：\n")
cat("   - 图片尺寸：宽度8英寸，高度5英寸\n")
cat("   - 数据水平：属（Genus）- 排除Unclassified和NA\n")
cat("   - 同一属的ASV已合并\n")
cat("   - 每个样本用自己的S计算Pielou和Margalef\n")
cat("   - 分面：Shannon, Pielou, Margalef（并排横版）\n")
cat("   - 每个小图：左侧岱海(DH)，右侧永定河(WLSH)\n")
cat("   - 标题在顶部，灰色背景strip\n")
cat("   - 箱线图 + 抖动点\n")
cat("   - 每个指数独立计算显著性（*, **, ***）\n")
cat("   - 每个分面纵坐标独立（scales = free_y）\n\n")
cat("📁 输出文件夹: Fig4B.eDNA\n")
cat("   输出文件:\n")
cat("   - Fig4B_eDNA_Diversity_Editable.pdf (横版分面图，8×5英寸)\n")
cat("   - Fig4B_eDNA_Diversity_Preview.png (预览图)\n")
cat("   - Fig4B_eDNA_Shannon.pdf/png\n")
cat("   - Fig4B_eDNA_Pielou.pdf/png\n")
cat("   - Fig4B_eDNA_Margalef.pdf/png\n\n")
cat("按照论文公式计算（排除Unclassified）：\n")
cat("  公式(2-3) Shannon: H = -∑P_i × log₂(P_i)\n")
cat("  公式(2-4) Margalef: D = (S-1) / log₂(N)\n")
cat("  公式(2-5) Pielou: J = H / ln(S)\n")
cat("  注：S和N均为每个样本自己的值，已排除Unclassified\n")

# =====================================================
# Source: 图4.A.R
# =====================================================

# =====================================================
# eDNA多样性分析 - 一个图，三个指数作为分面
# 按照论文公式：
# Shannon: -∑P_i × log₂(P_i)
# Margalef: (S-1) / log₂(N)
# Pielou: H / ln(S)
# 图结构：X轴为湖泊（DH, WLSH），分面1：三个指数，分面2：三个季节
# 保存文件夹：S2A
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/S2A")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 配色
lake_colors <- c("DH" = "#D6857E", "WLSH" = "#5A7D9A")

# =====================================================
# 按照论文公式定义函数
# =====================================================

# Shannon指数: 使用 log2
shannon_log2 <- function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(NA)
  p <- x / sum(x)
  -sum(p * log2(p))
}

# Pielou指数: H / ln(S)
pielou_formula <- function(x) {
  x <- x[x > 0]
  if(length(x) <= 1) return(NA)
  H <- shannon_log2(x)
  S <- length(x)
  H / log(S)
}

# Margalef指数: (S-1) / log₂(N)
margalef_log2 <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  N <- sum(x)
  if(N > 0) (S - 1) / log2(N) else NA
}

# =====================================================
# 计算多样性指数
# =====================================================

timepoints <- c("DH04","DH07","DH10","WLSH04","WLSH07","WLSH10")
season_map <- c("Spring","Summer","Autumn","Spring","Summer","Autumn")
names(season_map) <- timepoints

all_diversity <- data.frame()

for(tp in timepoints) {
  asv_file <- file.path(data_dir, tp, "ASV.txt")
  if(!file.exists(asv_file)) next
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE)
  sample_cols <- setdiff(colnames(asv), c("ASV_ID","taxonomy"))
  if(length(sample_cols) == 0) next
  
  asv_t <- as.data.frame(t(asv[,sample_cols]))
  colnames(asv_t) <- asv$ASV_ID
  
  lake <- ifelse(grepl("DH", tp), "DH", "WLSH")
  season <- season_map[tp]
  
  cat("处理:", tp, "-", lake, "-", season, "- 样本数:", nrow(asv_t), "\n")
  
  for(s in rownames(asv_t)){
    abu <- as.numeric(asv_t[s,])
    
    sha <- shannon_log2(abu)
    pie <- pielou_formula(abu)
    mar <- margalef_log2(abu)
    
    all_diversity <- rbind(all_diversity, data.frame(
      Sample = s, Lake = lake, Season = season,
      Shannon = sha, Pielou = pie, Margalef = mar
    ))
  }
}

# 查看结果统计
cat("\n多样性指数统计（按论文公式）:\n")
cat("Shannon: H = -∑P_i × log₂(P_i)\n")
cat("Margalef: D = (S-1) / log₂(N)\n")
cat("Pielou: J = H / ln(S)\n\n")
print(summary(all_diversity[,c("Shannon","Pielou","Margalef")]))

# 转换成长格式，便于分面
all_long <- all_diversity %>%
  pivot_longer(cols = c(Shannon, Pielou, Margalef), 
               names_to = "Index", 
               values_to = "Value")

# 设置因子顺序
all_long$Season <- factor(all_long$Season, levels = c("Spring", "Summer", "Autumn"))
all_long$Lake <- factor(all_long$Lake, levels = c("DH", "WLSH"))
all_long$Index <- factor(all_long$Index, levels = c("Shannon", "Pielou", "Margalef"))

# =====================================================
# 计算显著性（按指数和季节）
# =====================================================

calculate_significance <- function(data, index_name, season_list) {
  sig_df <- data.frame()
  
  for(s in season_list) {
    sub <- data %>% filter(Index == index_name, Season == s)
    dh <- sub %>% filter(Lake == "DH") %>% pull(Value)
    wl <- sub %>% filter(Lake == "WLSH") %>% pull(Value)
    
    if(length(dh) < 2 | length(wl) < 2) next
    
    pv <- wilcox.test(dh, wl)$p.value
    if(pv >= 0.05) next
    
    sig_label <- ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", "*"))
    
    y_max <- max(sub$Value, na.rm = TRUE)
    y_min <- min(sub$Value, na.rm = TRUE)
    y_pos <- y_max + (y_max - y_min) * 0.15
    
    sig_df <- rbind(sig_df, data.frame(
      Index = index_name,
      Season = s,
      x = 1.5,
      y = y_pos,
      label = sig_label
    ))
  }
  
  return(sig_df)
}

seasons <- c("Spring", "Summer", "Autumn")
indices <- c("Shannon", "Pielou", "Margalef")

sig_all <- data.frame()
for(idx in indices) {
  sig_temp <- calculate_significance(all_long, idx, seasons)
  sig_all <- rbind(sig_all, sig_temp)
}

# =====================================================
# 创建分面图：行=指数，列=季节
# =====================================================

# 创建组合分面变量
all_long$facet_row <- factor(all_long$Index, 
                             levels = c("Shannon", "Pielou", "Margalef"),
                             labels = c("Shannon", "Pielou", "Margalef"))
all_long$facet_col <- factor(all_long$Season, 
                             levels = c("Spring", "Summer", "Autumn"))

# 计算每个分面的Y轴范围
y_ranges <- all_long %>%
  group_by(Index, Season) %>%
  summarise(
    y_max = max(Value, na.rm = TRUE),
    y_min = min(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_upper = y_max + (y_max - y_min) * 0.2
  )

# 主图
p <- ggplot(all_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 21, size = 2, alpha = 0.7, width = 0.15) +
  facet_grid(rows = vars(facet_row), cols = vars(facet_col), scales = "free_y") +
  scale_fill_manual(values = lake_colors, name = "Lake") +
  labs(x = NULL, y = "Value") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# 添加显著性标记
if(nrow(sig_all) > 0) {
  for(i in 1:nrow(sig_all)) {
    # 获取对应分面的Y轴范围
    y_range_info <- y_ranges %>% 
      filter(Index == sig_all$Index[i], Season == sig_all$Season[i])
    
    if(nrow(y_range_info) > 0) {
      y_offset <- (y_range_info$y_max - y_range_info$y_min) * 0.08
      
      p <- p +
        annotate("segment", 
                 x = 1, xend = 2, 
                 y = sig_all$y[i] - y_offset, 
                 yend = sig_all$y[i] - y_offset,
                 linewidth = 0.5) +
        annotate("text", 
                 x = sig_all$x[i], 
                 y = sig_all$y[i], 
                 label = sig_all$label[i],
                 size = 4.5, fontface = "bold")
    }
  }
}

# =====================================================
# 添加Y轴标签（为每一行单独设置）
# =====================================================

# 创建Y轴标签的data frame
y_labels <- data.frame(
  facet_row = c("Shannon", "Pielou", "Margalef"),
  y_label = c(
    expression("Shannon diversity (H'"[log2]*")"),
    expression("Pielou evenness (J'"[log2]*")"),
    expression("Margalef richness (D"[log2]*")")
  )
)

# 使用ggplot2的标签功能
p <- p + 
  ylab(NULL) +
  theme(strip.text.y = element_text(face = "bold", size = 11))

# =====================================================
# 保存为可编辑PDF和预览图
# =====================================================

# 保存PDF（文字可编辑，矢量格式）
pdf(file.path(output_dir, "S2A_Diversity_Editable.pdf"), 
    width = 12, height = 9, 
    useDingbats = FALSE)
print(p)
dev.off()

# 保存PNG用于预览
ggsave(file.path(output_dir, "S2A_Diversity_Preview.png"), 
       p, width = 12, height = 9, dpi = 300, bg = "white")

cat("\n✅ 完成！生成一个综合分面图：\n")
cat("   - 列分面：三个季节（Spring, Summer, Autumn）\n")
cat("   - 行分面：三个指数（Shannon, Pielou, Margalef）\n")
cat("   - 每个小图：左侧岱海(DH)，右侧永定河(WLSH)\n")
cat("   - 箱线图 + 抖动点\n")
cat("   - 显著性标记（*, **, ***）\n\n")
cat("📁 输出文件夹: S2A\n")
cat("   输出文件:\n")
cat("   - S2A_Diversity_Editable.pdf (一个图，3×3分面布局)\n")
cat("   - S2A_Diversity_Preview.png (预览图)\n")

# =====================================================
# Source: S6形态学.R
# =====================================================

# =====================================================
# 形态学数据：修复 + 多样性分析 + 出图
# 按照论文公式：
#    Shannon: H = -∑P_i × log₂(P_i)
#    Margalef: D = (S-1) / log₂(N)
#    Pielou: J = H / ln(S)
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magick)


# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/S6形态")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 配色
lake_colors <- c("DH"="#D6857E", "WLSH"="#5A7D9A")

# 时间点
timepoints <- c("DH04", "DH07", "DH10", "WLSH04", "WLSH07", "WLSH10")

# =====================================================
# 按照论文公式定义函数
# =====================================================

# Shannon指数: 使用 log2 (公式2-3)
shannon_log2 <- function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(NA)
  p <- x / sum(x)
  -sum(p * log2(p))
}

# Pielou指数: H / ln(S) (公式2-5)
pielou_formula <- function(x) {
  x <- x[x > 0]
  if(length(x) <= 1) return(NA)
  H <- shannon_log2(x)
  S <- length(x)
  H / log(S)  # ln(S)
}

# Margalef指数: (S-1) / log₂(N) (公式2-4)
margalef_log2 <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  N <- sum(x)
  if(N > 0) (S - 1) / log2(N) else NA
}

# =====================================================
# 第一步：修复所有形态学文件（合并重复物种）
# =====================================================

cat("========== 第一步：修复数据文件 ==========\n")

fix_and_merge <- function(file_path) {
  morpho <- read.table(file_path, 
                       header = TRUE, 
                       sep = "\t", 
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       fill = TRUE)
  
  tp_name <- basename(dirname(file_path))
  cat("\n处理:", tp_name, "\n")
  cat("  原始行数:", nrow(morpho), "\n")
  
  species_names <- morpho[[1]]
  
  # 移除空物种名的行
  empty_rows <- which(species_names == "" | is.na(species_names))
  if(length(empty_rows) > 0) {
    cat("  移除空物种名行:", length(empty_rows), "行\n")
    morpho <- morpho[-empty_rows, ]
    species_names <- morpho[[1]]
  }
  
  # 合并重复的物种名
  dup_names <- unique(species_names[duplicated(species_names)])
  
  if(length(dup_names) > 0) {
    cat("  发现重复物种名:", paste(dup_names, collapse=", "), "\n")
    
    for(dup in dup_names) {
      dup_rows <- which(species_names == dup)
      
      if(length(dup_rows) > 1) {
        sample_cols <- 2:ncol(morpho)
        
        for(i in 2:length(dup_rows)) {
          morpho[dup_rows[1], sample_cols] <- as.numeric(morpho[dup_rows[1], sample_cols]) + 
            as.numeric(morpho[dup_rows[i], sample_cols])
        }
        
        morpho <- morpho[-dup_rows[-1], ]
        cat("    合并:", dup, "- 共", length(dup_rows), "行\n")
        species_names <- morpho[[1]]
      }
    }
  }
  
  write.table(morpho, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  保存成功，新行数:", nrow(morpho), "\n")
  
  return(morpho)
}

# 修复所有文件
for(tp in timepoints) {
  file_path <- file.path(data_dir, tp, "morpho.txt")
  if(file.exists(file_path)) {
    fix_and_merge(file_path)
  } else {
    cat("⚠️ 文件不存在:", file_path, "\n")
  }
}

# =====================================================
# 第二步：读取修复后的数据并计算多样性
# =====================================================

cat("\n========== 第二步：计算多样性指数 ==========\n")

get_lake <- function(tp) {
  if(grepl("^DH", tp)) return("DH")
  if(grepl("^WLSH", tp)) return("WLSH")
  return(NA)
}

get_season <- function(tp) {
  if(grepl("04$", tp)) return("Spring")
  if(grepl("07$", tp)) return("Summer")
  if(grepl("10$", tp)) return("Autumn")
  return(NA)
}

all_diversity <- data.frame()

for(tp in timepoints) {
  file_path <- file.path(data_dir, tp, "morpho.txt")
  
  if(!file.exists(file_path)) {
    cat("⚠️ 缺失:", file_path, "\n")
    next
  }
  
  morpho <- read.table(file_path, 
                       header = TRUE, 
                       sep = "\t", 
                       check.names = FALSE,
                       row.names = 1,
                       stringsAsFactors = FALSE)
  
  morpho <- as.data.frame(lapply(morpho, function(x) as.numeric(as.character(x))))
  morpho[is.na(morpho)] <- 0
  
  lake <- get_lake(tp)
  season <- get_season(tp)
  
  cat("✓", tp, "- 湖泊:", lake, "- 季节:", season, "- 样本数:", ncol(morpho), "- 物种数:", nrow(morpho), "\n")
  
  for(sample_name in colnames(morpho)) {
    abu <- morpho[[sample_name]]
    abu <- abu[abu > 0]
    
    if(length(abu) > 0) {
      shannon <- shannon_log2(abu)
      richness <- length(abu)
      pielou <- pielou_formula(abu)
      margalef <- margalef_log2(abu)
    } else {
      shannon <- NA
      richness <- 0
      pielou <- NA
      margalef <- NA
    }
    
    all_diversity <- rbind(all_diversity, data.frame(
      Sample = sample_name,
      Lake = lake,
      Season = season,
      Shannon = shannon,
      Richness = richness,
      Pielou = pielou,
      Margalef = margalef,
      stringsAsFactors = FALSE
    ))
  }
}

# 检查数据
cat("\n数据统计:\n")
print(table(all_diversity$Lake, all_diversity$Season))

if(nrow(all_diversity) == 0) {
  stop("没有找到任何数据文件！")
}

# =====================================================
# 第三步：转为长格式
# =====================================================

all_long <- all_diversity %>%
  pivot_longer(cols = c(Shannon, Pielou, Margalef), 
               names_to = "Index", 
               values_to = "Value") %>%
  filter(!is.na(Value))

all_long$Season <- factor(all_long$Season, levels = c("Spring", "Summer", "Autumn"))
all_long$Index <- factor(all_long$Index, levels = c("Shannon", "Pielou", "Margalef"))
all_long$Lake <- factor(all_long$Lake, levels = c("DH", "WLSH"))

all_long$facet <- factor(
  paste(all_long$Season, all_long$Index, sep = " - "),
  levels = c(
    "Spring - Shannon", "Spring - Pielou", "Spring - Margalef",
    "Summer - Shannon", "Summer - Pielou", "Summer - Margalef",
    "Autumn - Shannon", "Autumn - Pielou", "Autumn - Margalef"
  )
)

# =====================================================
# 第四步：计算显著性
# =====================================================

cat("\n========== 第三步：显著性检验 ==========\n")

sig_data <- data.frame()

for(s in c("Spring", "Summer", "Autumn")) {
  for(i in c("Shannon", "Pielou", "Margalef")) {
    sub <- filter(all_long, Season == s, Index == i)
    
    if(nrow(sub) == 0) next
    
    dh_vals <- sub$Value[sub$Lake == "DH"]
    wl_vals <- sub$Value[sub$Lake == "WLSH"]
    
    if(length(dh_vals) < 2 | length(wl_vals) < 2) next
    
    pv <- wilcox.test(dh_vals, wl_vals)$p.value
    
    if(pv >= 0.05) next
    
    sig <- ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", "*"))
    
    y_max <- max(sub$Value, na.rm = TRUE)
    y_min <- min(sub$Value, na.rm = TRUE)
    y_pos <- y_max + (y_max - y_min) * 0.15
    
    sig_data <- rbind(sig_data, data.frame(
      facet = paste(s, i, sep = " - "),
      x = 1.5,
      y = y_pos,
      label = sig
    ))
    
    cat("  ", s, "-", i, ": p =", format(pv, scientific = TRUE), sig, "\n")
  }
}

if(nrow(sig_data) == 0) {
  cat("  无显著差异 (p < 0.05)\n")
}

# =====================================================
# 第五步：绘图 - 直接保存为可编辑PDF
# =====================================================

cat("\n========== 第四步：绘图并保存可编辑PDF ==========\n")

p <- ggplot(all_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 21, size = 1.4, alpha = 0.7, width = 0.1) +
  facet_wrap(~facet, nrow = 3, ncol = 3, scales = "free") +
  scale_fill_manual(values = lake_colors, name = "Lake") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10)
  )

if(nrow(sig_data) > 0) {
  p <- p +
    geom_text(data = sig_data, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 4.5, fontface = "bold") +
    geom_segment(data = sig_data, 
                 aes(x = 1, xend = 2, y = y - 0.05, yend = y - 0.05),
                 inherit.aes = FALSE, linewidth = 0.4)
}

# ========== 直接保存可编辑的PDF（不经过PNG裁剪） ==========
# 直接保存可编辑的PDF（不经过PNG裁剪）
# 使用默认字体（Windows下通常是 Arial 或 sans）
pdf(file.path(output_dir, "FigS6_Morphology_Editable.pdf"), 
    width = 10, height = 12, 
    family = "sans",  # 改成 "sans" 或直接删除这行
    useDingbats = FALSE)
print(p)
dev.off()

# 方法2：同时保存PNG版本（用于预览）
ggsave(file.path(output_dir, "FigS6_Morphology.png"), 
       p, width = 10, height = 12, dpi = 300, bg = "white")

cat("\n✅ 可编辑PDF已保存！\n")
cat("📁 文件位置:", file.path(output_dir, "FigS6_Morphology_Editable.pdf"), "\n")
cat("\n注意：这个PDF的文字在AI中可以直接编辑！\n")

# =====================================================
# 完成
# =====================================================

cat("\n========== 完成！ ==========\n")
cat("📁 输出目录:", output_dir, "\n")
cat("\n按照论文公式计算：\n")
cat("  公式(2-3) Shannon: H = -∑P_i × log₂(P_i)\n")
cat("  公式(2-4) Margalef: D = (S-1) / log₂(N)\n")
cat("  公式(2-5) Pielou: J = H / ln(S)\n\n")
cat("输出文件:\n")
cat("  - FigS6_Morphology.png (预览用)\n")
cat("  - FigS6_Morphology_Editable.pdf (AI可编辑文字)\n")
cat("\n顺序:\n")
cat("  第1行: Spring-Shannon, Spring-Pielou, Spring-Margalef\n")
cat("  第2行: Summer-Shannon, Summer-Pielou, Summer-Margalef\n")
cat("  第3行: Autumn-Shannon, Autumn-Pielou, Autumn-Margalef\n")

# =====================================================
# Source: S6多样性.R
# =====================================================

# =====================================================
# eDNA多样性分析 - 完全按照论文公式
# 公式(2-3): Shannon = -∑P_i × log₂(P_i)
# 公式(2-4): Margalef = (S-1) / log₂(N)
# 公式(2-5): Pielou = H / ln(S)
# 输出：可编辑PDF（AI中可直接修改文字）
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/S6")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 配色
lake_colors <- c("DH"="#D6857E", "WLSH"="#5A7D9A")

# =====================================================
# 按照论文公式定义函数
# =====================================================

# Shannon指数: 使用 log2 (公式2-3)
shannon_log2 <- function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(NA)
  p <- x / sum(x)
  -sum(p * log2(p))
}

# Pielou指数: H / ln(S) (公式2-5)
pielou_formula <- function(x) {
  x <- x[x > 0]
  if(length(x) <= 1) return(NA)
  H <- shannon_log2(x)
  S <- length(x)
  H / log(S)
}

# Margalef指数: (S-1) / log₂(N) (公式2-4)
margalef_log2 <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  N <- sum(x)
  if(N > 0) (S - 1) / log2(N) else NA
}

# =====================================================
# 计算多样性指数
# =====================================================

timepoints <- c("DH04","DH07","DH10","WLSH04","WLSH07","WLSH10")
season_map <- c("Spring","Summer","Autumn","Spring","Summer","Autumn")
names(season_map) <- timepoints

all_diversity <- data.frame()

for(tp in timepoints) {
  asv_file <- file.path(data_dir, tp, "ASV.txt")
  if(!file.exists(asv_file)) next
  
  asv <- read.table(asv_file, header=T, sep="\t", check.names=F)
  sample_cols <- setdiff(colnames(asv), c("ASV_ID","taxonomy"))
  if(length(sample_cols)==0) next
  
  asv_t <- as.data.frame(t(asv[,sample_cols]))
  colnames(asv_t) <- asv$ASV_ID
  
  lake <- ifelse(grepl("DH",tp),"DH","WLSH")
  season <- season_map[tp]
  
  cat("处理:", tp, "-", lake, "-", season, "- 样本数:", nrow(asv_t), "\n")
  
  for(s in rownames(asv_t)){
    abu <- as.numeric(asv_t[s,])
    
    sha <- shannon_log2(abu)
    ric <- specnumber(abu)
    pie <- pielou_formula(abu)
    mar <- margalef_log2(abu)
    
    all_diversity <- rbind(all_diversity, data.frame(
      Sample=s, Lake=lake, Season=season, 
      Shannon=sha, Pielou=pie, Margalef=mar,
      Richness=ric
    ))
  }
}

# 查看结果统计
cat("\n多样性指数统计（按论文公式）:\n")
cat("公式(2-3) Shannon: H = -∑P_i × log₂(P_i)\n")
cat("公式(2-4) Margalef: D = (S-1) / log₂(N)\n")
cat("公式(2-5) Pielou: J = H / ln(S)\n\n")
print(summary(all_diversity[,c("Shannon","Pielou","Margalef")]))

# =====================================================
# 转为长格式
# =====================================================

all_long <- all_diversity %>%
  pivot_longer(cols=c(Shannon,Pielou,Margalef), names_to="Index", values_to="Value") %>%
  filter(!is.na(Value))

# 设置因子顺序
all_long$Season <- factor(all_long$Season, levels = c("Spring", "Summer", "Autumn"))
all_long$Index <- factor(all_long$Index, levels = c("Shannon", "Pielou", "Margalef"))
all_long$Lake <- factor(all_long$Lake, levels = c("DH", "WLSH"))

# 创建分面变量
all_long$facet <- factor(
  paste(all_long$Season, all_long$Index, sep = " - "),
  levels = c(
    "Spring - Shannon", "Spring - Pielou", "Spring - Margalef",
    "Summer - Shannon", "Summer - Pielou", "Summer - Margalef",
    "Autumn - Shannon", "Autumn - Pielou", "Autumn - Margalef"
  )
)

# =====================================================
# 计算显著性
# =====================================================

sig_data <- data.frame()

for(s in c("Spring", "Summer", "Autumn")) {
  for(i in c("Shannon", "Pielou", "Margalef")) {
    sub <- filter(all_long, Season == s, Index == i)
    dh <- filter(sub, Lake == "DH")$Value
    wl <- filter(sub, Lake == "WLSH")$Value
    
    if(length(dh) < 2 | length(wl) < 2) next
    
    pv <- wilcox.test(dh, wl)$p.value
    if(pv >= 0.05) next
    
    sig <- ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", "*"))
    
    y_max <- max(sub$Value, na.rm = TRUE)
    y_min <- min(sub$Value, na.rm = TRUE)
    y_pos <- y_max + (y_max - y_min) * 0.15
    
    sig_data <- rbind(sig_data, data.frame(
      facet = paste(s, i, sep = " - "),
      x = 1.5,
      y = y_pos,
      label = sig
    ))
  }
}

# =====================================================
# 绘图
# =====================================================

p <- ggplot(all_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 21, size = 1.4, alpha = 0.7, width = 0.1) +
  facet_wrap(~facet, nrow = 3, ncol = 3, scales = "free") +
  scale_fill_manual(values = lake_colors, name = "Lake") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10)
  )

if(nrow(sig_data) > 0) {
  p <- p +
    geom_text(data = sig_data, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 4.5, fontface = "bold") +
    geom_segment(data = sig_data, 
                 aes(x = 1, xend = 2, y = y - 0.05, yend = y - 0.05),
                 inherit.aes = FALSE, linewidth = 0.4)
}

# =====================================================
# 直接保存为可编辑的PDF
# =====================================================

# 保存PDF（文字可编辑，矢量格式）
pdf(file.path(output_dir, "FigS6_Editable.pdf"), 
    width = 10, height = 12, 
    useDingbats = FALSE)
print(p)
dev.off()

# 同时保存PNG用于预览
ggsave(file.path(output_dir, "FigS6_Preview.png"), 
       p, width = 10, height = 12, dpi = 300, bg = "white")

cat("\n✅ 完成！完全按照论文公式：\n")
cat("   公式(2-3) Shannon: H = -∑P_i × log₂(P_i)\n")
cat("   公式(2-4) Margalef: D = (S-1) / log₂(N)\n")
cat("   公式(2-5) Pielou: J = H / ln(S)\n\n")
cat("📁 输出文件:\n")
cat("   - FigS6_Editable.pdf (AI中可编辑文字)\n")
cat("   - FigS6_Preview.png (预览用)\n\n")
cat("   布局顺序:\n")
cat("     第1行: Spring-Shannon, Spring-Pielou, Spring-Margalef\n")
cat("     第2行: Summer-Shannon, Summer-Pielou, Summer-Margalef\n")
cat("     第3行: Autumn-Shannon, Autumn-Pielou, Autumn-Margalef\n")
