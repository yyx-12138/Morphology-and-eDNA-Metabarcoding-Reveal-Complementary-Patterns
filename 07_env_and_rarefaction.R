# 07 Environmental factors and rarefaction


# =====================================================
# Source: 图2.R
# =====================================================

# =====================================================
# 水质参数两湖对比箱线图（3×4分面布局）
# 不分季节，每个小图左DH右WLSH，自由纵坐标
# 颜色：Daihai=#FE6100, Wuliangsuhai=#077E98
# =====================================================

rm(list=ls())

library(ggplot2)
library(reshape2)
library(patchwork)

# 设置路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/水质参数对比")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==================== 1. 定义颜色 ====================
lake_colors <- c(
  "DH" = "#D6857E",
  "WLSH" = "#5A7D9A"
)

# ==================== 2. 读取数据（使用 env.txt）====================
read_env_data <- function(timepoints, lake_name) {
  all_data <- data.frame()
  
  for(tp in timepoints) {
    env_file <- file.path(data_dir, tp, "env.txt")
    
    if(!file.exists(env_file)) {
      cat("文件不存在:", env_file, "\n")
      next
    }
    
    env <- read.table(env_file, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, check.names = FALSE)
    
    colnames(env)[1] <- "SiteID"
    env$SiteID <- gsub("-", "", env$SiteID)
    env$Lake <- lake_name
    env$TimePoint <- tp
    
    all_data <- rbind(all_data, env)
  }
  
  return(all_data)
}

# 定义时间点
dh_timepoints <- c("DH04", "DH07", "DH10")
wlsh_timepoints <- c("WLSH04", "WLSH07", "WLSH10")

# 读取数据
cat("========== 读取环境数据 ==========\n")
dh_env <- read_env_data(dh_timepoints, "Daihai")
wlsh_env <- read_env_data(wlsh_timepoints, "Wuliangsuhai")

# 合并两个湖的数据
all_env <- rbind(dh_env, wlsh_env)
cat("Daihai数据:", nrow(dh_env), "行\n")
cat("Wuliangsuhai数据:", nrow(wlsh_env), "行\n")

# ==================== 关键修改：把湖名改成 DH / WLSH ====================
all_env$Lake <- ifelse(all_env$Lake == "Daihai", "DH", "WLSH")

# ==================== 3. 定义12个水质参数 ====================
env_factors <- c("Temp", "DO", "Ph", "Sal", "Chl-a", "CODcr", 
                 "TN", "TP", "NO3--N", "NH4+-N", "PO43--P", "F")

env_labels <- c(
  "Temp" = "Temp (°C)",
  "DO" = "DO (mg/L)",
  "Ph" = "pH",
  "Sal" = "Sal (‰)",
  "Chl-a" = "Chl-a (μg/L)",
  "CODcr" = "CODcr (mg/L)",
  "TN" = "TN (mg/L)",
  "TP" = "TP (mg/L)",
  "NO3--N" = "NO₃⁻-N (mg/L)",
  "NH4+-N" = "NH₄⁺-N (mg/L)",
  "PO43--P" = "PO₄³⁻-P (mg/L)",
  "F" = "F⁻ (mg/L)"
)

# ==================== 4. 转换长格式 ====================
data_long <- melt(all_env, 
                  id.vars = c("SiteID", "Lake", "TimePoint"),
                  measure.vars = env_factors,
                  variable.name = "Parameter",
                  value.name = "Value")

# 确保顺序 DH -> WLSH
data_long$Lake <- factor(data_long$Lake, levels = c("DH", "WLSH"))

# 参数顺序（3×4布局）
param_order <- c("Temp", "DO", "Ph", "Sal", 
                 "Chl-a", "CODcr", "TN", "TP",
                 "NO3--N", "NH4+-N", "PO43--P", "F")
data_long$Parameter <- factor(data_long$Parameter, levels = param_order)

# ==================== 5. 创建基础图形 ====================
cat("========== 生成水质参数对比图 ==========\n")

# 创建图形（自由纵坐标）
p <- ggplot(data_long, aes(x = Lake, y = Value, fill = Lake)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(aes(color = Lake), width = 0.2, size = 1.5, alpha = 0.6, shape = 16) +
  scale_fill_manual(values = lake_colors) +
  scale_color_manual(values = lake_colors) +
  facet_wrap(~Parameter, ncol = 4, scales = "free_y",
             labeller = labeller(Parameter = env_labels)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "grey90", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 0, size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# ==================== 6. 添加显著性标记 ====================
# 为每个参数单独计算显著性并添加
params <- unique(data_long$Parameter)

# 创建显著性数据框
sig_data <- data.frame(
  Parameter = character(),
  x = numeric(),
  y = numeric(),
  label = character(),
  stringsAsFactors = FALSE
)

for(param in params) {
  param_data <- data_long[data_long$Parameter == param, ]
  
  # 提取两个湖的数据
  dh_vals <- param_data[param_data$Lake == "DH", "Value"]
  wlsh_vals <- param_data[param_data$Lake == "WLSH", "Value"]
  
  if(length(dh_vals) > 0 && length(wlsh_vals) > 0) {
    # Wilcoxon检验
    test_result <- wilcox.test(dh_vals, wlsh_vals)
    p_val <- test_result$p.value
    
    # 确定显著性标记
    if(p_val < 0.001) sig <- "***"
    else if(p_val < 0.01) sig <- "**"
    else if(p_val < 0.05) sig <- "*"
    else sig <- NA
    
    if(!is.na(sig)) {
      # 计算该参数的y轴位置
      y_max <- max(param_data$Value, na.rm = TRUE)
      y_min <- min(param_data$Value, na.rm = TRUE)
      y_range <- y_max - y_min
      y_pos <- y_max + y_range * 0.15
      
      sig_data <- rbind(sig_data, data.frame(
        Parameter = param,
        x = 1.5,
        y = y_pos,
        label = sig,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 添加显著性标记（使用geom_text，不继承任何aes）
if(nrow(sig_data) > 0) {
  p <- p + geom_text(data = sig_data, 
                     aes(x = x, y = y, label = label),
                     size = 5, fontface = "bold",
                     inherit.aes = FALSE)
}

# ==================== 7. 保存图片 ====================
pdf_file <- file.path(output_dir, "WaterQuality_Comparison_S1.pdf")
ggsave(pdf_file, p, width = 12, height = 10, dpi = 300)
cat("PDF已保存:", pdf_file, "\n")

png_file <- file.path(output_dir, "WaterQuality_Comparison_S2.png")
ggsave(png_file, p, width = 12, height = 10, dpi = 300)
cat("PNG已保存:", png_file, "\n")

# ==================== 8. 保存整理后的数据 ====================
write.table(all_env, file.path(output_dir, "water_quality_data.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n数据统计:\n")
print(table(all_env$Lake))

cat("\n========== 完成 ==========\n")

rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggsignif)

# 设置路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/水质参数对比")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==================== 1. 读取数据 ====================
read_all_env <- function() {
  all_data <- data.frame()
  timepoints <- c("DH04", "DH07", "DH10", "WLSH04", "WLSH07", "WLSH10")
  season_map <- c("DH04"="Spring", "DH07"="Summer", "DH10"="Autumn",
                  "WLSH04"="Spring", "WLSH07"="Summer", "WLSH10"="Autumn")
  for(tp in timepoints) {
    env_file <- file.path(data_dir, tp, "env.txt")
    if(!file.exists(env_file)) next
    env <- read.table(env_file, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, check.names = FALSE)
    colnames(env)[1] <- "SiteID"
    env$SiteID <- gsub("-", "", env$SiteID)
    env$TimePoint <- tp
    env$Season <- season_map[tp]
    env$Lake <- ifelse(grepl("DH", tp), "Daihai", "Wuliangsuhai")
    all_data <- rbind(all_data, env)
  }
  return(all_data)
}
all_env <- read_all_env()

# ==================== 2. WQI/TLI 计算函数 ====================
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
  if(!is.na(chla_val) && chla_val > 0) { TLI_Chla <- 10 * (2.5 + 1.086 * log(chla_val)) } else { TLI_Chla <- NA }
  tn_val <- as.numeric(row["TN"])
  if(!is.na(tn_val) && tn_val > 0) { TLI_TN <- 10 * (5.453 + 1.694 * log(tn_val)) } else { TLI_TN <- NA }
  tp_val <- as.numeric(row["TP"])
  if(!is.na(tp_val) && tp_val > 0) { TLI_TP <- 10 * (9.436 + 1.624 * log(tp_val)) } else { TLI_TP <- NA }
  if(!is.na(TLI_Chla) && !is.na(TLI_TP) && !is.na(TLI_TN)) {
    return(0.4205 * TLI_Chla + 0.2967 * TLI_TP + 0.2828 * TLI_TN)
  } else { return(NA) }
}

# ==================== 3. 数据整理 ====================
all_env$WQI <- apply(all_env, 1, calculate_WQI)
all_env$TLI <- apply(all_env, 1, calculate_TLI)
all_env$WQI <- as.numeric(all_env$WQI)
all_env$TLI <- as.numeric(all_env$TLI)
all_env$Group <- paste(all_env$Season, all_env$Lake, sep = "_")
all_env$Group <- gsub("_Daihai", "-DH", all_env$Group)
all_env$Group <- gsub("_Wuliangsuhai", "-WLSH", all_env$Group)

group_order <- c("Spring-DH", "Summer-DH", "Autumn-DH", "Spring-WLSH", "Summer-WLSH", "Autumn-WLSH")
all_env$Group <- factor(all_env$Group, levels = group_order)

six_colors <- c(
  "Spring-DH"     = "#E8B5B0",
  "Summer-DH"     = "#D6857E",
  "Autumn-DH"     = "#B86B64",
  "Spring-WLSH"   = "#91A3B8",
  "Summer-WLSH"   = "#5A7D9A",
  "Autumn-WLSH"   = "#3F5A73"
)

group_labels <- c("Spring\nDH", "Summer\nDH", "Autumn\nDH", "Spring\nWLSH", "Summer\nWLSH", "Autumn\nWLSH")

# 图例标签：DH Spring 格式
legend_labels <- c("DH Spring", "DH Summer", "DH Autumn", "WLSH Spring", "WLSH Summer", "WLSH Autumn")

wqi_data <- all_env[!is.na(all_env$WQI), c("Group", "WQI")]
colnames(wqi_data)[2] <- "Value"
tli_data <- all_env[!is.na(all_env$TLI), c("Group", "TLI")]
colnames(tli_data)[2] <- "Value"

# ==================== 4. 显著性计算 ====================
all_comparisons <- list(
  c("Spring-DH", "Summer-DH"), c("Spring-DH", "Autumn-DH"), c("Summer-DH", "Autumn-DH"),
  c("Spring-WLSH", "Summer-WLSH"), c("Spring-WLSH", "Autumn-WLSH"), c("Summer-WLSH", "Autumn-WLSH"),
  c("Spring-DH", "Spring-WLSH"), c("Summer-DH", "Summer-WLSH"), c("Autumn-DH", "Autumn-WLSH"),
  c("Spring-DH", "Summer-WLSH"), c("Spring-DH", "Autumn-WLSH"),
  c("Summer-DH", "Spring-WLSH"), c("Summer-DH", "Autumn-WLSH"),
  c("Autumn-DH", "Spring-WLSH"), c("Autumn-DH", "Summer-WLSH")
)

get_sig_list <- function(data) {
  sig_list <- list()
  for(comp in all_comparisons) {
    g1 <- data[data$Group == comp[1], "Value"]
    g2 <- data[data$Group == comp[2], "Value"]
    if(length(g1) > 0 && length(g2) > 0) {
      pval <- wilcox.test(g1, g2)$p.value
      if(pval < 0.001) lab <- "***"
      else if(pval < 0.01) lab <- "**"
      else if(pval < 0.05) lab <- "*"
      else lab <- NA
      if(!is.na(lab)) {
        x1 <- which(group_order == comp[1])
        x2 <- which(group_order == comp[2])
        sig_list[[length(sig_list)+1]] <- list(x1 = x1, x2 = x2, label = lab)
      }
    }
  }
  return(sig_list)
}

assign_levels <- function(sig_list) {
  if(length(sig_list) == 0) return(list(sig_list = sig_list, n_levels = 0))
  sig_list <- sig_list[order(sapply(sig_list, function(x) x$x1))]
  levels <- list()
  for(sig in sig_list) {
    placed <- FALSE
    for(i in seq_along(levels)) {
      overlap <- any(sapply(levels[[i]], function(s) !(sig$x2 < s$x1 || sig$x1 > s$x2)))
      if(!overlap) { levels[[i]] <- c(levels[[i]], list(sig)); placed <- TRUE; break }
    }
    if(!placed) levels[[length(levels)+1]] <- list(sig)
  }
  result <- list()
  for(i in seq_along(levels)) {
    for(sig in levels[[i]]) {
      sig$level <- i
      result[[length(result)+1]] <- sig
    }
  }
  return(list(sig_list = result, n_levels = length(levels)))
}

wqi_sig_raw <- get_sig_list(wqi_data)
tli_sig_raw <- get_sig_list(tli_data)
wqi_res <- assign_levels(wqi_sig_raw)
tli_res <- assign_levels(tli_sig_raw)
wqi_sig_list <- wqi_res$sig_list
tli_sig_list <- tli_res$sig_list

cat("WQI显著对比数:", length(wqi_sig_list), "/15，层数:", wqi_res$n_levels, "\n")
cat("TLI显著对比数:", length(tli_sig_list), "/15，层数:", tli_res$n_levels, "\n")

# ==================== 5. 分面图（WQI在左，TLI在右）====================
combined_data <- bind_rows(
  wqi_data %>% mutate(Index = "Water Quality Index (WQI)"),
  tli_data %>% mutate(Index = "Trophic Level Index (TLI)")
)

# 设置因子顺序：WQI在前，TLI在后
combined_data$Index <- factor(combined_data$Index, levels = c("Water Quality Index (WQI)", "Trophic Level Index (TLI)"))

p <- ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, shape = 16) +
  facet_wrap(~Index, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = six_colors, labels = legend_labels) +
  scale_x_discrete(labels = group_labels) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

# WQI 显著性（左图）
if(length(wqi_sig_list) > 0) {
  step <- 2.2
  base_y <- 93
  for(i in 1:length(wqi_sig_list)) {
    s <- wqi_sig_list[[i]]
    ypos <- base_y + (s$level - 1) * step
    p <- p + ggsignif::geom_signif(
      comparisons = list(c(group_order[s$x1], group_order[s$x2])),
      y_position = ypos, annotations = s$label, tip_length = 0.01,
      size = 0.4, textsize = 3.5, vjust = 0.5,
      data = combined_data %>% filter(Index == "Water Quality Index (WQI)")
    )
  }
}

# TLI 显著性（右图）
if(length(tli_sig_list) > 0) {
  step <- 1.4
  base_y <- 69
  for(i in 1:length(tli_sig_list)) {
    s <- tli_sig_list[[i]]
    ypos <- base_y + (s$level - 1) * step
    p <- p + ggsignif::geom_signif(
      comparisons = list(c(group_order[s$x1], group_order[s$x2])),
      y_position = ypos, annotations = s$label, tip_length = 0.01,
      size = 0.4, textsize = 3.5, vjust = 0.5,
      data = combined_data %>% filter(Index == "Trophic Level Index (TLI)")
    )
  }
}

# ==================== 保存 ====================
ggsave(file.path(output_dir, "WQI_TLI_facet_final.png"), p, width = 14, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "WQI_TLI_facet_final.pdf"), p, width = 14, height = 5.5)

cat("\n✅ 完成！左边WQI，右边TLI，图例格式为 DH Spring 等\n")

# =====================================================
# Source: S1S2环境因子.R
# =====================================================

# =====================================================
# 环境因子季节动态箱线图（分面图风格）+ Kruskal-Wallis显著性标识
# 显著性标识以星号线形式标注在箱线图上方中央位置（整体下移）
# WLSH和DH分开画，分面图，横坐标一致，右侧统一图例
# 无标题，字体放大两号，无网格线
# 输出各季节平均值、标准差等统计量
# =====================================================

rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(rstatix)
library(tidyr)
library(ggsignif)

# 设置路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/S1S2环境因子")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==================== 1. 定义颜色 ====================
season_colors <- c(
  "Spring" = "#8DA0CB",
  "Summer" = "#FC8D62",
  "Autumn" = "#66C2A5"
)

# ==================== 2. 读取数据函数 ====================
read_env_data <- function(lake, timepoints, season_map) {
  all_data <- data.frame()
  
  for(tp in timepoints) {
    env_file <- file.path(data_dir, tp, "env.txt")
    if(!file.exists(env_file)) {
      cat("文件不存在:", env_file, "\n")
      next
    }
    
    env <- read.table(env_file, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, check.names = FALSE)
    
    colnames(env)[1] <- "SiteID"
    env$Season <- season_map[tp]
    env$Lake <- lake
    env$TimePoint <- tp
    
    all_data <- rbind(all_data, env)
  }
  
  return(all_data)
}

# ==================== 3. 定义时间点和季节映射 ====================
wlsh_timepoints <- c("WLSH04", "WLSH07", "WLSH10")
wlsh_season_map <- c("WLSH04" = "Spring", "WLSH07" = "Summer", "WLSH10" = "Autumn")

dh_timepoints <- c("DH04", "DH07", "DH10")
dh_season_map <- c("DH04" = "Spring", "DH07" = "Summer", "DH10" = "Autumn")

# ==================== 4. 读取数据 ====================
cat("========== 读取环境数据 ==========\n")
wlsh_env <- read_env_data("WLSH", wlsh_timepoints, wlsh_season_map)
dh_env <- read_env_data("DH", dh_timepoints, dh_season_map)
cat("WLSH环境数据行数:", nrow(wlsh_env), "\n")
cat("DH环境数据行数:", nrow(dh_env), "\n")

# ==================== 5. 定义环境因子 ====================
env_factors <- c("Temp", "DO", "Ph", "Sal", "Chl-a", "CODcr", 
                 "TN", "TP", "NO3--N", "NH4+-N", "PO43--P", "F")

# 标准简写（用于分面标签）
env_short_labels <- c(
  "Temp" = "Temp (°C)",
  "DO" = "DO (mg/L)",
  "Ph" = "pH",
  "Sal" = "Sal (‰)",
  "Chl-a" = "Chl-a (μg/L)",
  "CODcr" = "CODcr (mg/L)",
  "TN" = "TN (mg/L)",
  "TP" = "TP (mg/L)",
  "NO3--N" = "NO₃⁻-N (mg/L)",
  "NH4+-N" = "NH₄⁺-N (mg/L)",
  "PO43--P" = "PO₄³⁻-P (mg/L)",
  "F" = "F⁻ (mg/L)"
)

# 因子中文名称（用于输出表格）
env_chinese_names <- c(
  "Temp" = "水温",
  "DO" = "溶解氧",
  "Ph" = "pH",
  "Sal" = "盐度",
  "Chl-a" = "叶绿素a",
  "CODcr" = "化学需氧量",
  "TN" = "总氮",
  "TP" = "总磷",
  "NO3--N" = "硝酸盐氮",
  "NH4+-N" = "氨氮",
  "PO43--P" = "磷酸盐",
  "F" = "氟化物"
)

# ==================== 6. 计算统计量和显著性 ====================
calculate_stats_and_signif <- function(data, lake_name) {
  
  # 转换数据为长格式
  data_long <- melt(data, 
                    id.vars = c("SiteID", "Season", "Lake", "TimePoint"),
                    measure.vars = env_factors,
                    variable.name = "Factor",
                    value.name = "Value")
  
  # 确保Season顺序
  data_long$Season <- factor(data_long$Season, levels = c("Spring", "Summer", "Autumn"))
  
  # ========== 6.1 计算各季节平均值、标准差、标准误 ==========
  stats_summary <- data_long %>%
    group_by(Lake, Factor, Season) %>%
    summarise(
      N = n(),
      Mean = round(mean(Value, na.rm = TRUE), 3),
      SD = round(sd(Value, na.rm = TRUE), 3),
      SE = round(sd(Value, na.rm = TRUE) / sqrt(n()), 3),
      Min = round(min(Value, na.rm = TRUE), 3),
      Max = round(max(Value, na.rm = TRUE), 3),
      .groups = "drop"
    )
  
  # ========== 6.2 Kruskal-Wallis检验 + Dunn事后检验 ==========
  # 存储整体KW检验结果
  kw_results <- data.frame()
  # 存储两两比较结果（用于画图）
  pairwise_for_plot <- data.frame()
  
  for(factor in env_factors) {
    factor_data <- data_long %>% 
      filter(Factor == factor) %>%
      filter(!is.na(Value))
    
    if(nrow(factor_data) < 3 || length(unique(factor_data$Season)) < 2) next
    
    # Kruskal-Wallis检验
    kw_test <- kruskal.test(Value ~ Season, data = factor_data)
    
    # 存储整体KW结果
    kw_results <- rbind(kw_results, data.frame(
      Lake = lake_name,
      Factor = factor,
      Factor_CN = env_chinese_names[factor],
      Kruskal_p = round(kw_test$p.value, 5),
      Kruskal_signif = case_when(
        kw_test$p.value < 0.001 ~ "***",
        kw_test$p.value < 0.01 ~ "**",
        kw_test$p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ))
    
    # 如果整体显著，进行Dunn事后检验
    if(kw_test$p.value < 0.05) {
      dunn_res <- tryCatch({
        factor_data %>% dunn_test(Value ~ Season, p.adjust.method = "bonferroni")
      }, error = function(e) NULL)
      
      if(!is.null(dunn_res) && nrow(dunn_res) > 0) {
        for(i in 1:nrow(dunn_res)) {
          if(dunn_res$p.adj[i] < 0.05) {
            season_levels <- c("Spring", "Summer", "Autumn")
            x1 <- which(season_levels == dunn_res$group1[i])
            x2 <- which(season_levels == dunn_res$group2[i])
            
            signif_label <- case_when(
              dunn_res$p.adj[i] < 0.001 ~ "***",
              dunn_res$p.adj[i] < 0.01 ~ "**",
              dunn_res$p.adj[i] < 0.05 ~ "*"
            )
            
            pairwise_for_plot <- rbind(pairwise_for_plot, data.frame(
              Factor = factor,
              group1 = dunn_res$group1[i],
              group2 = dunn_res$group2[i],
              xmin = x1,
              xmax = x2,
              label = signif_label,
              p_adj = dunn_res$p.adj[i],
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  return(list(stats = stats_summary, kw_results = kw_results, pairwise = pairwise_for_plot))
}

# ==================== 7. 计算两个湖泊的结果 ====================
cat("\n========== 计算统计量和显著性 ==========\n")
wlsh_results <- calculate_stats_and_signif(wlsh_env, "Wuliangsuhai")
dh_results <- calculate_stats_and_signif(dh_env, "Daihai")

# ==================== 8. 输出统计表格到控制台 ====================
cat("\n========== WLSH 各季节环境因子统计表 ==========\n")
print(wlsh_results$stats)
cat("\n========== WLSH Kruskal-Wallis检验结果 ==========\n")
print(wlsh_results$kw_results)
cat("\n========== WLSH 两两比较显著性 ==========\n")
if(nrow(wlsh_results$pairwise) > 0) print(wlsh_results$pairwise) else cat("无显著两两差异\n")

cat("\n========== DH 各季节环境因子统计表 ==========\n")
print(dh_results$stats)
cat("\n========== DH Kruskal-Wallis检验结果 ==========\n")
print(dh_results$kw_results)
cat("\n========== DH 两两比较显著性 ==========\n")
if(nrow(dh_results$pairwise) > 0) print(dh_results$pairwise) else cat("无显著两两差异\n")

# ==================== 9. 保存统计表格到CSV ====================
write.csv(wlsh_results$stats, 
          file.path(output_dir, "WLSH_env_seasonal_stats.csv"), 
          row.names = FALSE)
write.csv(wlsh_results$kw_results, 
          file.path(output_dir, "WLSH_env_kruskal_results.csv"), 
          row.names = FALSE)
write.csv(wlsh_results$pairwise, 
          file.path(output_dir, "WLSH_env_pairwise_signif.csv"), 
          row.names = FALSE)
write.csv(dh_results$stats, 
          file.path(output_dir, "DH_env_seasonal_stats.csv"), 
          row.names = FALSE)
write.csv(dh_results$kw_results, 
          file.path(output_dir, "DH_env_kruskal_results.csv"), 
          row.names = FALSE)
write.csv(dh_results$pairwise, 
          file.path(output_dir, "DH_env_pairwise_signif.csv"), 
          row.names = FALSE)

cat("\n✓ 统计表格已保存到:", output_dir, "\n")

# ==================== 10. 绘制箱线图函数（显著性标识整体下移） ====================
plot_env_boxplot <- function(data, lake_name, pairwise_df, output_dir) {
  
  # 转换数据为长格式
  data_long <- melt(data, 
                    id.vars = c("SiteID", "Season", "Lake", "TimePoint"),
                    measure.vars = env_factors,
                    variable.name = "Factor",
                    value.name = "Value")
  
  # 确保Season的顺序
  data_long$Season <- factor(data_long$Season, levels = c("Spring", "Summer", "Autumn"))
  
  # 因子顺序（3×4布局，按行排列）
  factor_order <- c("Temp", "DO", "Ph", "Sal", 
                    "Chl-a", "CODcr", "TN", "TP",
                    "NO3--N", "NH4+-N", "PO43--P", "F")
  data_long$Factor <- factor(data_long$Factor, levels = factor_order)
  
  # 创建基础图形
  p <- ggplot(data_long, aes(x = Season, y = Value, fill = Season)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, shape = 16) +
    facet_wrap(~Factor, ncol = 4, scales = "free_y",
               labeller = labeller(Factor = env_short_labels)) +
    scale_fill_manual(values = season_colors, name = "Season") +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey70"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 0, size = 10),
      axis.text.y = element_text(size = 9),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  # 添加显著性标识（整体下移）
  if(nrow(pairwise_df) > 0) {
    # 为每个因子单独添加
    for(factor_name in unique(pairwise_df$Factor)) {
      factor_pairwise <- pairwise_df[pairwise_df$Factor == factor_name, ]
      
      # 获取当前因子的y范围
      factor_data <- data_long[data_long$Factor == factor_name, ]
      y_max <- max(factor_data$Value, na.rm = TRUE)
      y_min <- min(factor_data$Value, na.rm = TRUE)
      y_range <- y_max - y_min
      
      # 构建comparisons列表
      signif_map <- list()
      for(i in 1:nrow(factor_pairwise)) {
        comp <- c(as.character(factor_pairwise$group1[i]), 
                  as.character(factor_pairwise$group2[i]))
        signif_map[[i]] <- comp
      }
      
      # 提取标签
      signif_labels <- factor_pairwise$label
      
      # 计算y位置 - 整体下移，起始位置降低
      # 从 y_max + y_range * 0.02 开始（原来是 0.08）
      y_positions <- y_max + y_range * seq(0.02, 0.02 + 0.05*(length(signif_labels)-1), by = 0.05)
      
      # 添加显著性标识
      p <- p + geom_signif(
        data = subset(data_long, Factor == factor_name),
        comparisons = signif_map,
        map_signif_level = FALSE,
        annotations = signif_labels,
        y_position = y_positions,
        textsize = 3.5,        # 星号字体大小
        tip_length = 0.015,     # 横线长度
        vjust = 0.5,
        step_increase = 0.05
      )
    }
  }
  
  # 保存PDF
  pdf_file <- file.path(output_dir, paste0(lake_name, "_S1_env_seasonal_boxplot.pdf"))
  ggsave(pdf_file, p, width = 12, height = 10, dpi = 300)
  cat("PDF已保存:", pdf_file, "\n")
  
  # 保存PNG
  png_file <- file.path(output_dir, paste0(lake_name, "_S2_env_seasonal_boxplot.png"))
  ggsave(png_file, p, width = 12, height = 10, dpi = 300)
  cat("PNG已保存:", png_file, "\n")
  
  return(p)
}

# ==================== 11. 生成WLSH图 ====================
cat("\n========== 生成WLSH箱线图 ==========\n")
wlsh_plot <- plot_env_boxplot(wlsh_env, "Wuliangsuhai", wlsh_results$pairwise, output_dir)

# ==================== 12. 生成DH图 ====================
cat("\n========== 生成DH箱线图 ==========\n")
dh_plot <- plot_env_boxplot(dh_env, "Daihai", dh_results$pairwise, output_dir)

# ==================== 13. 输出综合对比表（两湖季节平均值对比） ====================
cat("\n========== 两湖环境因子季节平均值对比 ==========\n")

# 宽格式输出，方便对比
compare_table <- rbind(
  wlsh_results$stats %>% mutate(Lake = "Wuliangsuhai"),
  dh_results$stats %>% mutate(Lake = "Daihai")
) %>%
  select(Lake, Factor, Season, Mean, SD, N) %>%
  arrange(Factor, Season, Lake)

print(compare_table)

# 保存对比表
write.csv(compare_table, 
          file.path(output_dir, "TwoLakes_env_seasonal_comparison.csv"), 
          row.names = FALSE)

# ==================== 14. 输出样本数统计 ====================
cat("\n========== 样本数统计 ==========\n")
cat("WLSH各季节样本数:\n")
print(table(wlsh_env$Season))
cat("\nDH各季节样本数:\n")
print(table(dh_env$Season))

# ==================== 15. 保存整理后的原始数据 ====================
write.table(wlsh_env, file.path(output_dir, "WLSH_env_data.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dh_env, file.path(output_dir, "DH_env_data.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n✓ 整理后的数据已保存到:", output_dir, "\n")
cat("\n========== 完成 ==========\n")

# =====================================================
# Source: 稀释曲线.R
# =====================================================

# =====================================================
# eDNA 稀释曲线（Rarefaction Curve）
# 数据：数据/时间点/ASV.txt
# 输出：Fig_Rarefaction_eDNA.pdf/png
# 布局：3 个季节分面（Spring / Summer / Autumn），DH 红、WLSH 蓝
# =====================================================

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/稀释曲线")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 配色（与现有图一致）
lake_colors <- c("DH" = "#C82423", "WLSH" = "#5A7D9A")

timepoints <- c("DH04","DH07","DH10","WLSH04","WLSH07","WLSH10")
season_map <- c("Spring","Summer","Autumn","Spring","Summer","Autumn")
names(season_map) <- timepoints

# =====================================================
# 对单个样本矩阵（行=ASV，列=样本）计算稀释曲线
# 返回长格式：Sample, Lake, Season, Depth, Richness
# =====================================================
calc_rarefaction <- function(comm, lake, season, step = 2000) {
  comm <- as.matrix(comm)
  comm <- comm[rowSums(comm) > 0, , drop = FALSE]
  if(ncol(comm) == 0) return(NULL)
  
  depths <- colSums(comm)
  depths <- depths[depths > 0]
  if(length(depths) == 0) return(NULL)
  
  res_list <- list()
  for(s in names(depths)) {
    n_max <- depths[s]
    # 采样深度序列：小样本逐点，大样本按 step 取点（至少包含 1 和 n_max）
    if(n_max <= 100) {
      depth_seq <- 1:n_max
    } else {
      depth_seq <- unique(c(seq(1, n_max, by = step), n_max))
    }
    # vegan::rarefy 计算每个深度下的期望 ASV 数
    rich <- as.numeric(rarefy(t(comm[, s, drop = FALSE]), sample = depth_seq))
    res_list[[s]] <- data.frame(
      Sample = s,
      Lake = lake,
      Season = season,
      Depth = depth_seq,
      Richness = rich,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res_list)
}

# =====================================================
# 计算所有样本的稀释曲线
# =====================================================
cat("========== eDNA 稀释曲线 ==========\n")

edna_rare <- data.frame()

for(tp in timepoints) {
  asv_file <- file.path(data_dir, tp, "ASV.txt")
  if(!file.exists(asv_file)) { cat("缺失:", asv_file, "\n"); next }
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t",
                    check.names = FALSE, stringsAsFactors = FALSE)
  sample_cols <- setdiff(colnames(asv), c("ASV_ID","taxonomy"))
  if(length(sample_cols) == 0) next
  
  asv_matrix <- as.matrix(asv[, sample_cols, drop = FALSE])
  rownames(asv_matrix) <- asv$ASV_ID
  asv_matrix <- asv_matrix[rowSums(asv_matrix) > 0, , drop = FALSE]
  
  lake <- ifelse(grepl("^DH", tp), "DH", "WLSH")
  season <- season_map[tp]
  
  rare_tp <- calc_rarefaction(asv_matrix, lake, season, step = 2000)
  if(!is.null(rare_tp)) {
    edna_rare <- rbind(edna_rare, rare_tp)
    cat("✓", tp, "- 样本:", ncol(asv_matrix), "- ASV:", nrow(asv_matrix), "\n")
  }
}

if(nrow(edna_rare) == 0) stop("没有可用的 eDNA 数据！")

edna_rare$Season <- factor(edna_rare$Season, levels = c("Spring","Summer","Autumn"))
edna_rare$Lake <- factor(edna_rare$Lake, levels = c("DH","WLSH"))

# =====================================================
# 绘图
# =====================================================
p_edna <- ggplot(edna_rare, aes(x = Depth, y = Richness, color = Lake, group = Sample)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~Season, nrow = 1, scales = "free_x") +
  scale_color_manual(values = lake_colors, name = "Lake") +
  labs(x = "Number of reads", y = "Observed ASV richness") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold")
  )

# 保存
pdf(file.path(output_dir, "Fig_Rarefaction_eDNA.pdf"),
    width = 12, height = 4.5, useDingbats = FALSE)
print(p_edna)
dev.off()

ggsave(file.path(output_dir, "Fig_Rarefaction_eDNA.png"),
       p_edna, width = 12, height = 4.5, dpi = 300, bg = "white")

cat("\n✅ 完成！输出文件：\n")
cat("  ", file.path(output_dir, "Fig_Rarefaction_eDNA.pdf"), "\n")
cat("  ", file.path(output_dir, "Fig_Rarefaction_eDNA.png"), "\n")

# =====================================================
# 不分面的稀释曲线（两个版本）
# 前置：已运行前面的 edna_rare 计算部分
# =====================================================

# ---------- 版本一：所有样本叠加，一图多线，颜色区分湖泊 ----------
p_all <- ggplot(edna_rare, aes(x = Depth, y = Richness, color = Lake, group = Sample)) +
  geom_line(linewidth = 0.4, alpha = 0.5) +
  scale_color_manual(values = lake_colors, name = "Lake") +
  labs(x = "Number of reads", y = "Observed ASV richness") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "Fig_Rarefaction_eDNA_all.pdf"),
       p_all, width = 8, height = 5, useDingbats = FALSE)
ggsave(file.path(output_dir, "Fig_Rarefaction_eDNA_all.png"),
       p_all, width = 8, height = 5, dpi = 300, bg = "white")


# ---------- 版本二：按 湖泊×季节 聚合为均值曲线 + 误差带 ----------
rare_summary <- edna_rare %>%
  group_by(Lake, Season, Depth) %>%
  summarise(
    mean_rich = mean(Richness, na.rm = TRUE),
    se_rich   = sd(Richness, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 用线型区分季节（颜色区分湖泊）
p_mean <- ggplot(rare_summary, aes(x = Depth, y = mean_rich,
                                   color = Lake, fill = Lake,
                                   linetype = Season, group = interaction(Lake, Season))) +
  geom_ribbon(aes(ymin = mean_rich - se_rich, ymax = mean_rich + se_rich),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = lake_colors, name = "Lake") +
  scale_fill_manual(values = lake_colors, name = "Lake") +
  scale_linetype_manual(values = c(Spring = "solid", Summer = "dashed", Autumn = "dotted"),
                        name = "Season") +
  labs(x = "Number of reads", y = "Observed ASV richness") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "Fig_Rarefaction_eDNA_mean.pdf"),
       p_mean, width = 8, height = 5, useDingbats = FALSE)
ggsave(file.path(output_dir, "Fig_Rarefaction_eDNA_mean.png"),
       p_mean, width = 8, height = 5, dpi = 300, bg = "white")

cat("\n✅ 不分面稀释曲线已保存：\n")
cat("  版本一（所有样本）: Fig_Rarefaction_eDNA_all.pdf/png\n")
cat("  版本二（均值+误差带）: Fig_Rarefaction_eDNA_mean.pdf/png\n")
