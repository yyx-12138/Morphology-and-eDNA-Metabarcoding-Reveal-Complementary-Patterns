# 03 Density and abundance


# =====================================================
# Source: S5.R
# =====================================================

# =====================================================
# 计算原生动物密度显著性
# =====================================================

rm(list=ls())

library(dplyr)
library(tidyr)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")

# =====================================================
# 第一步：读取所有形态学数据并计算总个体数
# =====================================================

sample_config <- list(
  list(folder = "DH04", prefix = "C4DH", season = "Spring", lake = "DH", 
       samples = c("01", "02", "03", "04", "05", "06", "07", "09")),
  list(folder = "DH07", prefix = "X7DH", season = "Summer", lake = "DH",
       samples = c("01", "02", "03", "04", "05", "06", "08", "09")),
  list(folder = "DH10", prefix = "Q10DH", season = "Autumn", lake = "DH",
       samples = c("01", "02", "03", "04", "05", "06", "07", "08", "09")),
  list(folder = "WLSH04", prefix = "C4WLSH", season = "Spring", lake = "WLSH",
       samples = as.character(3:21)),
  list(folder = "WLSH07", prefix = "X7WLSH", season = "Summer", lake = "WLSH",
       samples = as.character(3:21)),
  list(folder = "WLSH10", prefix = "Q10WLSH", season = "Autumn", lake = "WLSH",
       samples = as.character(1:21))
)

all_counts <- data.frame()

for(cfg in sample_config) {
  file_path <- file.path(data_dir, cfg$folder, "morpho.txt")
  
  if(!file.exists(file_path)) {
    cat("⚠️ 文件不存在:", file_path, "\n")
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
  
  for(sample_num in cfg$samples) {
    sample_name <- paste0(cfg$prefix, sample_num)
    
    if(sample_name %in% colnames(morpho)) {
      total_inds <- sum(morpho[[sample_name]], na.rm = TRUE)
      
      all_counts <- rbind(all_counts, data.frame(
        Sample = sample_name,
        Lake = cfg$lake,
        Season = cfg$season,
        Total_Individuals = total_inds,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat("✓", cfg$folder, "- 样本数:", length(cfg$samples), "\n")
}

# =====================================================
# 第二步：密度参数（V1, V2, V3）
# =====================================================

params <- data.frame(
  Sample = c(
    "C4DH01", "C4DH02", "C4DH03", "C4DH04", "C4DH05", "C4DH06", "C4DH07", "C4DH09",
    "X7DH01", "X7DH02", "X7DH03", "X7DH04", "X7DH05", "X7DH06", "X7DH08", "X7DH09",
    "Q10DH01", "Q10DH02", "Q10DH03", "Q10DH04", "Q10DH05", "Q10DH06", "Q10DH07", "Q10DH08", "Q10DH09",
    "C4WLSH3", "C4WLSH4", "C4WLSH5", "C4WLSH6", "C4WLSH7", "C4WLSH8", "C4WLSH9",
    "C4WLSH10", "C4WLSH11", "C4WLSH12", "C4WLSH13", "C4WLSH14", "C4WLSH15", "C4WLSH16",
    "C4WLSH17", "C4WLSH18", "C4WLSH19", "C4WLSH20", "C4WLSH21",
    "X7WLSH3", "X7WLSH4", "X7WLSH5", "X7WLSH6", "X7WLSH7", "X7WLSH8", "X7WLSH9",
    "X7WLSH10", "X7WLSH11", "X7WLSH12", "X7WLSH13", "X7WLSH14", "X7WLSH15", "X7WLSH16",
    "X7WLSH17", "X7WLSH18", "X7WLSH19", "X7WLSH20", "X7WLSH21",
    "Q10WLSH1", "Q10WLSH2", "Q10WLSH3", "Q10WLSH4", "Q10WLSH5", "Q10WLSH6", "Q10WLSH7",
    "Q10WLSH8", "Q10WLSH9", "Q10WLSH10", "Q10WLSH11", "Q10WLSH12", "Q10WLSH13", "Q10WLSH14",
    "Q10WLSH15", "Q10WLSH16", "Q10WLSH17", "Q10WLSH18", "Q10WLSH19", "Q10WLSH20", "Q10WLSH21"
  ),
  V1 = c(
    2, 1, 1, 1, 1, 1, 1, 1,
    1, 0.5, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 1, 1, 1, 1, 1, 5, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5,
    1, 0.1, 0.1, 0.1, 0.1, 1, 0.5, 0.1, 1, 0.1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
    1, 1, 0.1, 1, 0.5, 1, 1, 1, 1, 0.5, 0.1, 1, 1, 0.5, 1, 0.5, 1, 1, 1, 1, 1
  ),
  V2 = 100,
  V3 = 30
)

# =====================================================
# 第三步：合并并计算密度
# =====================================================

density_data <- all_counts %>%
  left_join(params, by = "Sample") %>%
  mutate(
    Density = Total_Individuals * (V2 / V1) / V3
  ) %>%
  filter(!is.na(Density))

density_data$Season <- factor(density_data$Season, levels = c("Spring", "Summer", "Autumn"))

# =====================================================
# 第四步：统计检验
# =====================================================

# 分别对两个湖泊进行统计检验
for(lake_name in c("DH", "WLSH")) {
  cat("\n\n")
  cat("==========================================\n")
  cat("湖泊:", lake_name, "\n")
  cat("==========================================\n")
  
  lake_data <- density_data[density_data$Lake == lake_name, ]
  
  # 描述性统计
  cat("\n【描述性统计】\n")
  desc_stats <- lake_data %>%
    group_by(Season) %>%
    summarise(
      n = n(),
      Mean = round(mean(Density, na.rm = TRUE), 2),
      SD = round(sd(Density, na.rm = TRUE), 2),
      Median = round(median(Density, na.rm = TRUE), 2),
      Min = round(min(Density, na.rm = TRUE), 2),
      Max = round(max(Density, na.rm = TRUE), 2)
    )
  print(desc_stats)
  
  # 正态性检验（Shapiro-Wilk）
  cat("\n【正态性检验（Shapiro-Wilk）】\n")
  shapiro_res <- lake_data %>%
    group_by(Season) %>%
    summarise(
      W = round(shapiro.test(Density)$statistic, 4),
      p_value = shapiro.test(Density)$p.value,
      normal = ifelse(p_value > 0.05, "是", "否")
    )
  print(shapiro_res)
  
  # 方差齐性检验（Levene检验）
  library(car)
  levene_res <- leveneTest(Density ~ Season, data = lake_data)
  cat("\n【方差齐性检验（Levene）】\n")
  cat("F =", round(levene_res$`F value`[1], 4), "\n")
  cat("p =", levene_res$`Pr(>F)`[1], "\n")
  cat("方差齐性:", ifelse(levene_res$`Pr(>F)`[1] > 0.05, "是", "否"), "\n")
  
  # 选择参数或非参数检验
  # 检查是否满足正态性和方差齐性
  all_normal <- all(shapiro_res$p_value > 0.05)
  var_equal <- levene_res$`Pr(>F)`[1] > 0.05
  
  cat("\n【检验方法选择】\n")
  cat("所有组符合正态分布:", all_normal, "\n")
  cat("方差齐性:", var_equal, "\n")
  
  if(all_normal & var_equal) {
    cat("使用 ANOVA + Tukey HSD\n")
    
    # ANOVA
    anova_res <- aov(Density ~ Season, data = lake_data)
    cat("\n【ANOVA结果】\n")
    print(summary(anova_res))
    
    # Tukey HSD事后检验
    tukey_res <- TukeyHSD(anova_res)
    cat("\n【Tukey HSD事后检验】\n")
    tukey_df <- as.data.frame(tukey_res$Season)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$signif <- ifelse(tukey_df$`p adj` < 0.001, "***",
                              ifelse(tukey_df$`p adj` < 0.01, "**",
                                     ifelse(tukey_df$`p adj` < 0.05, "*", "n.s.")))
    print(tukey_df[, c("comparison", "p adj", "signif")])
    
  } else {
    cat("使用 Kruskal-Wallis + Dunn检验\n")
    
    # Kruskal-Wallis检验
    kruskal_res <- kruskal.test(Density ~ Season, data = lake_data)
    cat("\n【Kruskal-Wallis结果】\n")
    cat("Chi-squared =", kruskal_res$statistic, "\n")
    cat("df =", kruskal_res$parameter, "\n")
    cat("p-value =", kruskal_res$p.value, "\n")
    
    # Dunn事后检验
    if(kruskal_res$p.value < 0.05) {
      library(FSA)
      dunn_res <- dunnTest(Density ~ Season, data = lake_data, method = "bonferroni")
      cat("\n【Dunn事后检验（Bonferroni校正）】\n")
      dunn_df <- dunn_res$res
      dunn_df$signif <- ifelse(dunn_df$P.adj < 0.001, "***",
                               ifelse(dunn_df$P.adj < 0.01, "**",
                                      ifelse(dunn_df$P.adj < 0.05, "*", "n.s.")))
      print(dunn_df[, c("Comparison", "P.adj", "signif")])
    } else {
      cat("\n组间无显著差异（p > 0.05）\n")
    }
  }
  
  # 输出每个季节的密度值（方便查看）
  cat("\n【原始密度数据】\n")
  lake_data %>%
    select(Season, Density) %>%
    arrange(Season) %>%
    print()
}

cat("\n\n========== 计算完成 ==========\n")

# =====================================================
# Source: eDNA丰富度比较.R
# =====================================================

# =====================================================
# 最终完整版：eDNA 总丰度 + Margalef丰富度 全套计算
# 适配结构：数据/DH04/ASV.txt  数据/WLSH04/ASV.txt ...
# =====================================================
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)

# ==================== 固定路径（和你项目完全一致）====================
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/最终版_eDNA_丰度_丰富度")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==================== 季节配色（统一论文风格）====================
season_colors <- c(
  "Spring" = "#8DA0CB",
  "Summer" = "#FC8D62",
  "Autumn" = "#66C2A5"
)

# ==================== 所有采样文件夹配置 ====================
time_config <- list(
  list(dir = "DH04",   lake = "DH",   season = "Spring"),
  list(dir = "DH07",   lake = "DH",   season = "Summer"),
  list(dir = "DH10",   lake = "DH",   season = "Autumn"),
  list(dir = "WLSH04", lake = "WLSH", season = "Spring"),
  list(dir = "WLSH07", lake = "WLSH", season = "Summer"),
  list(dir = "WLSH10", lake = "WLSH", season = "Autumn")
)

# ==================== 循环读取所有 ASV 并计算 ====================
final_result <- data.frame()

for(cfg in time_config){
  asv_file <- file.path(data_dir, cfg$dir, "ASV.txt")
  
  if(!file.exists(asv_file)){
    cat("⚠️ 跳过不存在文件：", asv_file, "\n")
    next
  }
  
  cat("✅ 正在计算：", asv_file, "\n")
  asv <- read.delim(asv_file, header=T, sep="\t", check.names=F, stringsAsFactors=F)
  
  # 自动匹配样本列
  sample_cols <- grep("^C4|^X7|^Q10|^DH|^WLSH", colnames(asv), value=T)
  
  # 转为长数据
  asv_long <- asv %>%
    select(all_of(sample_cols), ASV_ID) %>%
    pivot_longer(-ASV_ID, names_to="Sample", values_to="Reads")
  
  # 按样本计算：总丰度 + 物种数 + Margalef
  sample_calc <- asv_long %>%
    group_by(Sample) %>%
    summarise(
      Total_Reads = sum(Reads, na.rm=T),          # 总丰度（用于比较）
      S = sum(Reads > 0, na.rm=T),                # ASV 物种数
      .groups = "drop"
    ) %>%
    mutate(
      N = Total_Reads,                            # 公式用 N=总reads
      Margalef = (S - 1) / log(N),                # 严格按你给的公式
      Lake = cfg$lake,
      Season = cfg$season
    )
  
  final_result <- rbind(final_result, sample_calc)
}

# ==================== 统一季节顺序 ====================
final_result$Season <- factor(final_result$Season, 
                              levels = c("Spring","Summer","Autumn"))

# ==================== 导出完整表格 ====================
write.csv(final_result,
          file.path(output_dir, "eDNA_总丰度_Margalef_完整结果.csv"),
          row.names = FALSE)

cat("\n🎉 计算全部完成！结果表已保存：\n")
print(file.path(output_dir, "eDNA_总丰度_Margalef_完整结果.csv"))
print(head(final_result, 10))

# =====================================================
# 图1：总丰度比较箱线图
# =====================================================
p1 <- ggplot(final_result, aes(x=Season, y=Total_Reads, fill=Season)) +
  geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.8) +
  geom_jitter(width=0.18, size=1.6, alpha=0.7) +
  facet_wrap(~Lake, ncol=2, scales="free_y") +
  scale_fill_manual(values=season_colors) +
  labs(x=NULL, y="Total eDNA Reads (Relative Abundance)") +
  theme_bw(base_size=14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face="bold", size=13),
    axis.text.x = element_text(face="bold", size=12),
    axis.text.y = element_text(face="bold", size=11)
  )

ggsave(file.path(output_dir, "图1_eDNA总丰度比较.pdf"), p1, 
       width=10, height=5.5, dpi=300, bg="white")

# =====================================================
# 图2：Margalef 丰富度比较箱线图
# =====================================================
p2 <- ggplot(final_result, aes(x=Season, y=Margalef, fill=Season)) +
  geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.8) +
  geom_jitter(width=0.18, size=1.6, alpha=0.7) +
  facet_wrap(~Lake, ncol=2, scales="free_y") +
  scale_fill_manual(values=season_colors) +
  labs(x=NULL, y="Margalef Richness Index") +
  theme_bw(base_size=14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face="bold", size=13),
    axis.text.x = element_text(face="bold", size=12),
    axis.text.y = element_text(face="bold", size=11)
  )

ggsave(file.path(output_dir, "图2_Margalef丰富度比较.pdf"), p2, 
       width=10, height=5.5, dpi=300, bg="white")

cat("\n✅ 两张论文图已保存完成！")
