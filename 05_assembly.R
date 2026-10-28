# 05 Community assembly


# =====================================================
# Source: 图6A.R
# =====================================================

# =====================================================
# 基于Stegen et al. (2013)方法的群落构建过程分析
# 包含βNTI和RCbray计算及ANOSIM季节性差异检验
# 置换次数：99次
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(ape)
library(picante)
library(ggpubr)

# ========== 路径设置 ==========
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/零模型")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ========== 配色 ==========
lake_colors <- c("DH" = "#C82423", "WLSH" = "#56B4E9")
season_colors <- c("Spring" = "#2E8B57", "Summer" = "#DAA520", "Autumn" = "#CD5C5C")

# ========== 文件夹映射 ==========
time_list <- c("DH04", "DH07", "DH10", "WLSH04", "WLSH07", "WLSH10")
season_map <- c(
  DH04 = "Spring", DH07 = "Summer", DH10 = "Autumn",
  WLSH04 = "Spring", WLSH07 = "Summer", WLSH10 = "Autumn"
)

# =====================================================
# 1. 读取数据函数
# =====================================================
read_edna_data <- function(tp) {
  asv_file <- file.path(data_dir, tp, "ASV.txt")
  if(!file.exists(asv_file)) return(NULL)
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE)
  sample_cols <- !colnames(asv) %in% c("ASV_ID", "taxonomy")
  asv_matrix <- as.matrix(asv[, sample_cols])
  rownames(asv_matrix) <- asv$ASV_ID
  asv_t <- as.data.frame(t(asv_matrix))
  asv_t[] <- lapply(asv_t, as.numeric)
  asv_t <- asv_t[rowSums(asv_t) > 0, , drop = FALSE]
  return(asv_t)
}

read_phylogenetic_tree <- function(tp) {
  tree_file <- file.path(data_dir, tp, "tree.nwk")
  if(file.exists(tree_file)) {
    return(read.tree(tree_file))
  }
  return(NULL)
}

# =====================================================
# 2. 计算βMNTD（参考Stegen et al. 2013）
# =====================================================
calc_beta_mntd <- function(comm, tree) {
  n_samples <- nrow(comm)
  if(n_samples < 2) return(NULL)
  
  common_taxa <- intersect(colnames(comm), tree$tip.label)
  if(length(common_taxa) < 3) {
    return(as.matrix(vegdist(comm, method = "bray")))
  }
  
  comm_sub <- comm[, common_taxa, drop = FALSE]
  tree_sub <- keep.tip(tree, common_taxa)
  phy_dist <- cophenetic(tree_sub)
  
  mntd_mat <- tryCatch({
    as.matrix(comdistnt(comm_sub, phy_dist, abundance.weighted = TRUE))
  }, error = function(e) {
    as.matrix(vegdist(comm_sub, method = "bray"))
  })
  return(mntd_mat)
}

# 单次置换（打乱物种丰度）
one_perm_beta_mntd <- function(comm, tree) {
  comm_null <- comm
  for(j in 1:ncol(comm_null)) {
    comm_null[, j] <- sample(comm_null[, j])
  }
  return(calc_beta_mntd(comm_null, tree))
}

# =====================================================
# 3. 计算βNTI（每个样本对单独计算零分布，99次置换）
# =====================================================
calc_bNTI_matrix <- function(comm, tree, n_perm = 99) {
  n_samples <- nrow(comm)
  if(n_samples < 3 || is.null(tree)) return(NULL)
  
  cat("    样本数:", n_samples, "\n")
  
  # 观测MNTD
  obs_mntd <- calc_beta_mntd(comm, tree)
  if(is.null(obs_mntd)) return(NULL)
  
  # 为每个样本对准备存储空间
  n_pairs <- n_samples * (n_samples - 1) / 2
  null_list <- vector("list", n_pairs)
  
  idx <- 1
  for(i in 1:(n_samples-1)) {
    for(j in (i+1):n_samples) {
      null_list[[idx]] <- numeric(n_perm)
      idx <- idx + 1
    }
  }
  
  # n_perm次置换
  cat("    开始", n_perm, "次置换...\n")
  for(perm in 1:n_perm) {
    if(perm %% 20 == 0) cat("      已完成", perm, "/", n_perm, "\n")
    
    perm_mntd <- one_perm_beta_mntd(comm, tree)
    if(!is.null(perm_mntd)) {
      idx <- 1
      for(i in 1:(n_samples-1)) {
        for(j in (i+1):n_samples) {
          null_list[[idx]][perm] <- perm_mntd[i, j]
          idx <- idx + 1
        }
      }
    }
  }
  
  # 计算每个样本对的βNTI
  bNTI_mat <- matrix(NA, n_samples, n_samples)
  rownames(bNTI_mat) <- rownames(comm)
  colnames(bNTI_mat) <- rownames(comm)
  
  idx <- 1
  for(i in 1:(n_samples-1)) {
    for(j in (i+1):n_samples) {
      null_vals <- null_list[[idx]]
      mean_null <- mean(null_vals, na.rm = TRUE)
      sd_null <- sd(null_vals, na.rm = TRUE)
      
      if(sd_null > 0) {
        bNTI_val <- (obs_mntd[i, j] - mean_null) / sd_null
      } else {
        bNTI_val <- 0
      }
      # 限制在[-5, 5]
      bNTI_val <- pmin(pmax(bNTI_val, -5), 5)
      
      bNTI_mat[i, j] <- bNTI_val
      bNTI_mat[j, i] <- bNTI_val
      idx <- idx + 1
    }
  }
  
  return(bNTI_mat)
}

# =====================================================
# 4. 计算RCbray（每个样本对单独计算零分布，99次置换）
# =====================================================
calc_RCbray_matrix <- function(comm, n_perm = 99) {
  n <- nrow(comm)
  if(n < 2) return(NULL)
  
  cat("    计算RCbray...\n")
  
  # 观测Bray-Curtis距离
  obs_bc <- as.matrix(vegdist(comm, method = "bray"))
  
  # 为每个样本对准备存储空间
  n_pairs <- n * (n - 1) / 2
  null_list <- vector("list", n_pairs)
  
  idx <- 1
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      null_list[[idx]] <- numeric(n_perm)
      idx <- idx + 1
    }
  }
  
  # n_perm次置换
  for(perm in 1:n_perm) {
    if(perm %% 20 == 0) cat("      已完成", perm, "/", n_perm, "\n")
    
    # 打乱丰度
    comm_null <- comm
    for(j in 1:ncol(comm_null)) {
      comm_null[, j] <- sample(comm_null[, j])
    }
    null_bc <- as.matrix(vegdist(comm_null, method = "bray"))
    
    idx <- 1
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        null_list[[idx]][perm] <- null_bc[i, j]
        idx <- idx + 1
      }
    }
  }
  
  # 计算每个样本对的RCbray
  RCbray_mat <- matrix(NA, n, n)
  rownames(RCbray_mat) <- rownames(comm)
  colnames(RCbray_mat) <- rownames(comm)
  
  idx <- 1
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      null_vals <- null_list[[idx]]
      mean_null <- mean(null_vals, na.rm = TRUE)
      sd_null <- sd(null_vals, na.rm = TRUE)
      
      if(sd_null > 0) {
        rc_val <- (obs_bc[i, j] - mean_null) / sd_null
      } else {
        rc_val <- 0
      }
      # 限幅到[-1, 1]
      rc_val <- pmin(pmax(rc_val, -1), 1)
      
      RCbray_mat[i, j] <- rc_val
      RCbray_mat[j, i] <- rc_val
      idx <- idx + 1
    }
  }
  
  return(RCbray_mat)
}

# =====================================================
# 5. 提取下三角样本对
# =====================================================
extract_lower_tri <- function(bNTI_mat, RCbray_mat, samples, season, lake) {
  n <- length(samples)
  result <- data.frame()
  
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      result <- rbind(result, data.frame(
        Season = season,
        Lake = lake,
        Sample1 = samples[i],
        Sample2 = samples[j],
        bNTI = bNTI_mat[i, j],
        RCbray = RCbray_mat[j, i],
        stringsAsFactors = FALSE
      ))
    }
  }
  return(result)
}

# =====================================================
# 6. 主循环（每个湖独立计算，99次置换）
# =====================================================

n_perm <- 99

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("eDNA零模型分析（基于Stegen et al. 2013方法）\n")
cat("置换次数:", n_perm, "\n")
cat(paste(rep("=", 60), collapse=""), "\n")

all_results <- data.frame()

for(tp in time_list) {
  season <- season_map[tp]
  lake <- ifelse(grepl("^DH", tp), "DH", "WLSH")
  
  cat("\n========================================\n")
  cat("处理:", tp, "(", lake, "-", season, ")\n")
  cat("========================================\n")
  
  # 读取数据
  comm <- read_edna_data(tp)
  if(is.null(comm)) {
    cat("  ❌ 跳过：ASV文件不存在\n")
    next
  }
  
  tree <- read_phylogenetic_tree(tp)
  if(is.null(tree)) {
    cat("  ❌ 跳过：树文件不存在\n")
    next
  }
  
  if(nrow(comm) < 3) {
    cat("  ❌ 跳过：样本数不足", nrow(comm), "\n")
    next
  }
  
  cat("  ✅ 样本数:", nrow(comm), "，ASV数:", ncol(comm), "\n")
  
  # 计算βNTI
  bNTI_mat <- calc_bNTI_matrix(comm, tree, n_perm)
  if(is.null(bNTI_mat)) {
    cat("  ❌ βNTI计算失败\n")
    next
  }
  
  # 计算RCbray
  RCbray_mat <- calc_RCbray_matrix(comm, n_perm)
  if(is.null(RCbray_mat)) {
    cat("  ❌ RCbray计算失败\n")
    next
  }
  
  # 提取结果
  samples <- rownames(comm)
  result_df <- extract_lower_tri(bNTI_mat, RCbray_mat, samples, season, lake)
  all_results <- rbind(all_results, result_df)
  
  cat("  ✅ 完成：", nrow(result_df), "个样本对\n")
}

# =====================================================
# 7. 检查结果
# =====================================================

if(nrow(all_results) == 0) {
  stop("❌ 没有生成任何数据！")
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("数据汇总\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("总样本对数:", nrow(all_results), "\n")
cat("bNTI范围:", round(range(all_results$bNTI, na.rm = TRUE), 3), "\n")
cat("RCbray范围:", round(range(all_results$RCbray, na.rm = TRUE), 3), "\n")

# =====================================================
# 8. ANOSIM分析（99次置换检验）评估群落季节性差异
# =====================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ANOSIM分析（99次置换检验）评估群落季节性差异\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# 准备Bray-Curtis距离矩阵用于ANOSIM
all_comm_list <- list()
all_metadata <- data.frame()

for(tp in time_list) {
  comm <- read_edna_data(tp)
  if(is.null(comm)) next
  
  season <- season_map[tp]
  lake <- ifelse(grepl("^DH", tp), "DH", "WLSH")
  
  rownames_old <- rownames(comm)
  rownames(comm) <- paste(tp, rownames_old, sep = "_")
  
  all_comm_list[[tp]] <- comm
  
  temp_meta <- data.frame(
    SampleID = rownames(comm),
    Season = season,
    Lake = lake,
    TimePoint = tp,
    stringsAsFactors = FALSE
  )
  all_metadata <- rbind(all_metadata, temp_meta)
}

# 合并所有ASV表
all_taxa <- unique(unlist(lapply(all_comm_list, colnames)))
full_comm <- matrix(0, nrow = nrow(all_metadata), ncol = length(all_taxa))
rownames(full_comm) <- all_metadata$SampleID
colnames(full_comm) <- all_taxa

for(tp in names(all_comm_list)) {
  comm_tp <- all_comm_list[[tp]]
  for(samp in rownames(comm_tp)) {
    otu_cols <- match(colnames(comm_tp), all_taxa)
    full_comm[samp, otu_cols] <- as.numeric(comm_tp[samp, ])
  }
}

cat("总样本数:", nrow(full_comm), "\n")
cat("总OTU数:", ncol(full_comm), "\n")

# 计算Bray-Curtis距离矩阵
bc_dist <- vegdist(full_comm, method = "bray")

# 按湖泊分别进行ANOSIM
anosim_results <- list()

for(lake in c("DH", "WLSH")) {
  cat("\n", lake, "湖:\n")
  lake_samples <- all_metadata$SampleID[all_metadata$Lake == lake]
  
  if(length(lake_samples) < 2) {
    cat("  样本数不足，跳过\n")
    next
  }
  
  bc_dist_matrix <- as.matrix(bc_dist)
  lake_dist_matrix <- bc_dist_matrix[lake_samples, lake_samples]
  lake_dist <- as.dist(lake_dist_matrix)
  
  lake_season <- all_metadata$Season[all_metadata$Lake == lake]
  
  set.seed(123)
  anosim_res <- anosim(lake_dist, lake_season, permutations = 99)
  anosim_results[[lake]] <- anosim_res
  
  cat("  R值:", round(anosim_res$statistic, 4), "\n")
  cat("  P值:", anosim_res$signif, "\n")
  
  if(anosim_res$signif < 0.01) {
    cat("  结论: 群落季节差异极显著 (p < 0.01)\n")
  } else if(anosim_res$signif < 0.05) {
    cat("  结论: 群落季节差异显著 (p < 0.05)\n")
  } else {
    cat("  结论: 群落季节差异不显著 (p > 0.05)\n")
  }
}

# =====================================================
# 9. 绘制散点图（三个季节分面）
# =====================================================

all_results$Season <- factor(all_results$Season, levels = c("Spring", "Summer", "Autumn"))

season_counts <- table(all_results$Season)

# 主图：βNTI vs RCbray
p1 <- ggplot(all_results, aes(x = bNTI, y = RCbray, color = Lake)) +
  geom_point(size = 1.8, alpha = 0.6) +
  scale_color_manual(values = lake_colors, name = "湖泊") +
  facet_wrap(~ Season, ncol = 3, 
             labeller = labeller(Season = function(x) {
               count <- ifelse(x %in% names(season_counts), season_counts[x], 0)
               paste(x, "(n =", count, "对, 99次置换)")
             })) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.3) +
  annotate("text", x = -3.5, y = 0.8, label = "同质选择", size = 3, color = "darkred") +
  annotate("text", x = 3.5, y = 0.8, label = "异质选择", size = 3, color = "darkred") +
  annotate("text", x = 0, y = 0.5, label = "随机性过程", size = 3, color = "gray30") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-1.2, 1.2)) +
  labs(x = expression(beta * NTI), 
       y = "RCbray",
       title = "群落构建过程 | βNTI vs RCbray (99次置换)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# 箱线图：βNTI的季节变化
p2 <- ggplot(all_results, aes(x = Season, y = bNTI, fill = Season)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") +
  scale_fill_manual(values = season_colors) +
  facet_wrap(~ Lake, ncol = 2) +
  labs(y = expression(beta * NTI), x = "季节",
       title = "βNTI季节变化 (99次置换)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# 箱线图：RCbray的季节变化
p3 <- ggplot(all_results, aes(x = Season, y = RCbray, fill = Season)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", color = "red") +
  scale_fill_manual(values = season_colors) +
  facet_wrap(~ Lake, ncol = 2) +
  labs(y = "RCbray", x = "季节",
       title = "RCbray季节变化 (99次置换)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# 组合图形
combined_plot <- (p1) / (p2 | p3) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

# 保存图片
output_file1 <- file.path(output_dir, "bNTI_RCbray_ThreeSeasons_99perm.png")
ggsave(output_file1, p1, width = 14, height = 5.5, dpi = 300)

output_file2 <- file.path(output_dir, "Boxplot_SeasonalVariation_99perm.png")
ggsave(output_file2, combined_plot, width = 14, height = 10, dpi = 300)

cat("\n✅ 图片保存至:\n")
cat("   -", output_file1, "\n")
cat("   -", output_file2, "\n")

# =====================================================
# 10. 群落构建过程统计
# =====================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("群落构建过程统计（基于Stegen et al. 2013方法）\n")
cat("注：|βNTI| > 2 为确定性过程主导；|βNTI| < 2 为随机性过程主导\n")
cat("置换次数：99次\n")
cat(paste(rep("=", 60), collapse=""), "\n")

process_summary <- data.frame()

for(season in c("Spring", "Summer", "Autumn")) {
  season_data <- all_results[all_results$Season == season, ]
  if(nrow(season_data) == 0) next
  
  for(lake in c("DH", "WLSH")) {
    lake_data <- season_data[season_data$Lake == lake, ]
    if(nrow(lake_data) == 0) next
    
    total <- nrow(lake_data)
    
    deterministic <- sum(abs(lake_data$bNTI) > 2, na.rm = TRUE)
    deterministic_pct <- deterministic / total * 100
    
    homogeneous <- sum(lake_data$bNTI < -2, na.rm = TRUE)
    homogeneous_pct <- homogeneous / total * 100
    
    heterogeneous <- sum(lake_data$bNTI > 2, na.rm = TRUE)
    heterogeneous_pct <- heterogeneous / total * 100
    
    stochastic <- sum(abs(lake_data$bNTI) <= 2, na.rm = TRUE)
    stochastic_pct <- stochastic / total * 100
    
    stoch_data <- lake_data[abs(lake_data$bNTI) <= 2, ]
    if(nrow(stoch_data) > 0) {
      homogenizing_dispersal <- sum(stoch_data$RCbray < -0.95, na.rm = TRUE)
      dispersal_limitation <- sum(stoch_data$RCbray > 0.95, na.rm = TRUE)
      undominated <- sum(abs(stoch_data$RCbray) <= 0.95, na.rm = TRUE)
      hd_pct <- homogenizing_dispersal / total * 100
      dl_pct <- dispersal_limitation / total * 100
      und_pct <- undominated / total * 100
    } else {
      homogenizing_dispersal <- 0
      dispersal_limitation <- 0
      undominated <- 0
      hd_pct <- 0
      dl_pct <- 0
      und_pct <- 0
    }
    
    temp_df <- data.frame(
      Season = season,
      Lake = lake,
      TotalPairs = total,
      Deterministic_pct = deterministic_pct,
      Homogeneous_pct = homogeneous_pct,
      Heterogeneous_pct = heterogeneous_pct,
      Stochastic_pct = stochastic_pct,
      HomogenizingDispersal_pct = hd_pct,
      DispersalLimitation_pct = dl_pct,
      Undominated_pct = und_pct
    )
    process_summary <- rbind(process_summary, temp_df)
    
    cat("\n", season, "-", lake, "湖 (", total, "个样本对):\n", sep = "")
    cat("  确定性过程: ", deterministic, " (", round(deterministic_pct, 1), "%)\n", sep = "")
    cat("    - 同质选择 (βNTI < -2): ", homogeneous, " (", round(homogeneous_pct, 1), "%)\n", sep = "")
    cat("    - 异质选择 (βNTI > 2): ", heterogeneous, " (", round(heterogeneous_pct, 1), "%)\n", sep = "")
    cat("  随机性过程: ", stochastic, " (", round(stochastic_pct, 1), "%)\n", sep = "")
    cat("    - 同质扩散 (RCbray < -0.95): ", homogenizing_dispersal, " (", round(hd_pct, 1), "%)\n", sep = "")
    cat("    - 扩散限制 (RCbray > 0.95): ", dispersal_limitation, " (", round(dl_pct, 1), "%)\n", sep = "")
    cat("    - 非主导过程 (|RCbray| < 0.95): ", undominated, " (", round(und_pct, 1), "%)\n", sep = "")
  }
}

# 保存统计结果
stats_file <- file.path(output_dir, "CommunityAssembly_Statistics_99perm.csv")
write.csv(process_summary, stats_file, row.names = FALSE)
cat("\n✅ 统计结果保存至:", stats_file, "\n")

# 保存原始数据
results_file <- file.path(output_dir, "bNTI_RCbray_Results_99perm.csv")
write.csv(all_results, results_file, row.names = FALSE)
cat("✅ 原始数据保存至:", results_file, "\n")

cat("\n🎉 全部完成！\n")

# ==============================
# 【最终版】600 DPI 高清去白底
# 处理：Fig7eDNA.tif + Fig7形态.tif
# 输出：超高清透明PNG（600DPI）
# ==============================

# 安装包（第一次运行）
install.packages("magick")
library(magick)

# --------------------------
# 1. 处理 Fig7eDNA.tif
# --------------------------
img1 <- image_read("E:/写东西、杂物/论文/Yan/测序文件整理/输出图片/Fig7eDNA.tif")

# 去白色背景
img1_trans <- image_transparent(img1, color = "white", fuzz = 10)

# 导出 600 DPI 高清透明图
image_write(
  img1_trans,
  path = "E:/写东西、杂物/论文/Yan/测序文件整理/输出图片/Fig7eDNA_600DPI_透明.png",
  format = "png",
  density = 600  # 600DPI 高清
)

# --------------------------
# 2. 处理 Fig7形态.tif
# --------------------------
img2 <- image_read("E:/写东西、杂物/论文/Yan/测序文件整理/输出图片/Fig7形态.tif")

# 去白色背景
img2_trans <- image_transparent(img2, color = "white", fuzz = 10)

# 导出 600 DPI 高清透明图
image_write(
  img2_trans,
  path = "E:/写东西、杂物/论文/Yan/测序文件整理/输出图片/Fig7形态_600DPI_透明.png",
  format = "png",
  density = 600  # 600DPI 高清
)

# --------------------------
# 完成提示
# --------------------------
cat("✅ 全部完成！\n")
cat("已生成 600 DPI 透明背景图片：\n")
cat("1. Fig7eDNA_600DPI_透明.png\n")
cat("2. Fig7形态_600DPI_透明.png\n")
