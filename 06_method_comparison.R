# 06 Method comparison


# =====================================================
# Source: 深度计算物种检出比例等.R
# =====================================================

# =====================================================
# 形态学 vs eDNA 方法学比较统计报告（修正版）
# 修正：形态学数据过滤掉全0行（未实际检出的物种）
# =====================================================

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)

# =====================================================
# 路径设置
# =====================================================
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/Fig3_最终版/方法学比较统计_修正后")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =====================================================
# 时间点列表
# =====================================================
timepoints <- list(
  list(folder = "WLSH04", lake = "WLSH", season = "Spring"),
  list(folder = "WLSH07", lake = "WLSH", season = "Summer"),
  list(folder = "WLSH10", lake = "WLSH", season = "Autumn"),
  list(folder = "DH04", lake = "DH", season = "Spring"),
  list(folder = "DH07", lake = "DH", season = "Summer"),
  list(folder = "DH10", lake = "DH", season = "Autumn")
)

# =====================================================
# 排除列表（eDNA专用）
# =====================================================
exclude_phyla <- c(
  "Dinoflagellata",
  "Peronosporomycetes", "Hyphochytriomycetes",
  "Cryptomycota", "Chytridiomycota", "Ascomycota",
  "Basidiomycota", "Aphelidea", "Mucoromycota", "Zoopagomycota",
  "Porphyridiophyceae", "Kathablepharidae", "Chlorophyta", 
  "Ochrophyta", "Bacillariophyta", "Rhodophyta", "Cryptophyta", 
  "Haptophyta", "Phragmoplastophyta", "Pavlovophyceae", "Prymnesiophyceae",
  "Labyrinthulomycetes", "MAST-2", "MAST-12", "MAST-6", "LKM15",
  "Arthropoda", "Rotifera", "Nematoda", "Annelida",
  "Streptophyta", "Magnoliophyta", "Pinophyta"
)

exclude_genera <- c(
  "Unclassified", "uncultured", "Pythium", "Phytophthora",
  "Porphyridium", "MAST-2", "A31", "Ceratium", "MAST-12C", "NA",
  "Symbiodinium", "Pfiesteria", "Woloszynskia", "Pelagodinium",
  "Durinskia", "Amphidinium", "Scrippsiella", "Biecheleria",
  "Katablepharis", "Ebria", "Glaucocystis"
)

exclude_families <- c(
  "MAST-2", "Thraustochytriaceae", "MAST-12A", "Incertae_Sedis", "MAST-12C",
  "Cryptomonadales", "Pavlovophyceae", "Chlorophyceae", "Trebouxiophyceae",
  "Mediophyceae", "Mamiellales", "Chlorodendrales", "Pedinellales",
  "Chromulinales", "Bacillariophyceae", "Fragilariales", "Chrysophyceae",
  "Ochromonadales", "Sphaeropleales", "Ulvophyceae", "Ochrophyta",
  "Pedinomonadales", "Mamiellophyceae", "Eustigmatales", "Prymnesiales",
  "Zygnematophyceae", "Thoracosphaeraceae", "Isochrysidales", "Chattonellales",
  "Gonyaulacales", "Woloszynskia", "Porphyridiophyceae",
  "Peronosporomycetes", "LKM11", "Malasseziaceae", "Pleosporaceae",
  "Didymellaceae", "Cladosporiaceae", "Aphelidea", "Trichomeriaceae",
  "Rhinosporideacae", "Rhizophydiaceae", "Schizophyllaceae", "Kappamycetaceae",
  "Ichthyophonidae", "Filobasidiaceae", "Chytridiaceae", "Aspergillaceae",
  "Gromochytriaceae", "Debaryomycetaceae", "Sporormiaceae", "Phaffomycetaceae",
  "Bulleribasidiaceae", "Pleurotaceae", "Hydropus", "Pythiaceae",
  "Magnoliophyta", "Pinophyta", "uncultured", "NA", "Calanoida"
)

# =====================================================
# 辅助函数
# =====================================================
extract_last_num <- function(s) {
  s_as_char <- as.character(s)
  nums <- regmatches(s_as_char, gregexpr("[0-9]+$", s_as_char))
  num_str <- unlist(nums)
  return(as.numeric(num_str))
}

# =====================================================
# 第一部分：读取形态学数据（过滤全0行）
# =====================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║    形态学 vs eDNA 方法学比较统计报告（修正版）                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("【数据读取中...】\n")

# 存储形态学数据
morpho_all_data <- data.frame()

for(tp in timepoints) {
  file_path <- file.path(data_dir, tp$folder, "morpho.txt")
  
  if(!file.exists(file_path)) {
    cat("⚠️ 文件不存在:", file_path, "\n")
    next
  }
  
  morpho <- read.table(file_path, 
                       header = TRUE, 
                       sep = "\t", 
                       check.names = FALSE,
                       stringsAsFactors = FALSE)
  
  genus_names <- morpho[[1]]
  genus_names <- gsub(" sp\\.?$", "", genus_names)
  genus_names <- trimws(genus_names)
  
  morpho_matrix <- as.matrix(morpho[, -1])
  rownames(morpho_matrix) <- genus_names
  mode(morpho_matrix) <- "numeric"
  morpho_matrix[is.na(morpho_matrix)] <- 0
  
  # 关键修正：过滤掉所有样本计数均为0的物种
  row_sums <- rowSums(morpho_matrix)
  morpho_matrix <- morpho_matrix[row_sums > 0, ]
  
  # 如果没有检出任何物种，跳过
  if(nrow(morpho_matrix) == 0) {
    cat("⚠️", tp$folder, "无检出物种\n")
    next
  }
  
  # 聚合相同属（如有同属不同种的情况）
  morpho_agg <- aggregate(morpho_matrix, by = list(Genus = rownames(morpho_matrix)), FUN = sum)
  
  morpho_long <- morpho_agg %>%
    pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Count") %>%
    mutate(
      Lake = tp$lake,
      Season = tp$season,
      Sample_full = paste(tp$folder, Sample, sep = "_")
    )
  
  morpho_all_data <- rbind(morpho_all_data, morpho_long)
  cat("✓ 形态学:", tp$folder, "- 实际检出", nrow(morpho_agg), "属\n")
}

# =====================================================
# 第二部分：读取 eDNA 数据（排除非原生动物）
# =====================================================

asv_data_raw <- data.frame()

for(tp in timepoints) {
  asv_file <- file.path(data_dir, tp$folder, "ASV.txt")
  tax_file <- file.path(data_dir, tp$folder, "tax.txt")
  
  if(!file.exists(asv_file) | !file.exists(tax_file)) {
    cat("⚠️ 文件不存在:", asv_file, "或", tax_file, "\n")
    next
  }
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  tax <- read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = c("NA", "", "uncultured"))
  
  tax$Phylum <- gsub("^p__", "", tax$Phylum)
  tax$Class <- gsub("^c__", "", tax$Class)
  tax$Order <- gsub("^o__", "", tax$Order)
  tax$Family <- gsub("^f__", "", tax$Family)
  tax$Genus <- gsub("^g__", "", tax$Genus)
  
  asv$ASV_ID <- as.character(asv$ASV_ID)
  tax$ASV_ID <- as.character(tax$ASV_ID)
  
  asv_counts <- if("taxonomy" %in% colnames(asv)) asv[, colnames(asv) != "taxonomy"] else asv
  
  asv_long <- asv_counts %>%
    pivot_longer(cols = -ASV_ID, names_to = "Sample", values_to = "Count") %>%
    filter(Count > 0) %>%
    left_join(tax, by = "ASV_ID") %>%
    filter(!(Phylum %in% exclude_phyla)) %>%
    filter(!is.na(Phylum)) %>%
    filter(!(Genus %in% exclude_genera)) %>%
    filter(!(Family %in% exclude_families)) %>%
    mutate(
      Lake = tp$lake,
      Season = tp$season,
      Sample_full = paste(tp$folder, Sample, sep = "_")
    )
  
  asv_data_raw <- rbind(asv_data_raw, asv_long)
  cat("✓ eDNA:", tp$folder, "- 保留ASV数:", length(unique(asv_long$ASV_ID)), "\n")
}

# =====================================================
# 按湖泊分别处理数据
# =====================================================

cat("\n【数据处理中...】\n")

# 形态学数据按湖泊拆分
morpho_wlsh <- morpho_all_data %>% filter(Lake == "WLSH")
morpho_dh <- morpho_all_data %>% filter(Lake == "DH")

# 形态学所有检出属（按湖泊合并去重）
morpho_all_wlsh <- unique(morpho_wlsh$Genus)
morpho_all_dh <- unique(morpho_dh$Genus)

cat("\n形态学实际检出属数（合并去重后）:\n")
cat("  WLSH:", length(morpho_all_wlsh), "属\n")
cat("  DH:", length(morpho_all_dh), "属\n")
cat("  WLSH检出属:", paste(morpho_all_wlsh, collapse = ", "), "\n")
cat("  DH检出属:", paste(morpho_all_dh, collapse = ", "), "\n")

# eDNA 属级别数据
edna_by_genus <- asv_data_raw %>%
  filter(!is.na(Genus), Genus != "", Genus != "Unclassified", Genus != "uncultured") %>%
  group_by(Lake, Genus) %>%
  summarise(TotalCount = sum(Count), .groups = "drop")

edna_all_wlsh <- edna_by_genus %>% filter(Lake == "WLSH") %>% pull(Genus) %>% unique()
edna_all_dh <- edna_by_genus %>% filter(Lake == "DH") %>% pull(Genus) %>% unique()

# eDNA ASV级别统计
asv_wlsh <- asv_data_raw %>% filter(Lake == "WLSH")
asv_dh <- asv_data_raw %>% filter(Lake == "DH")

asv_total_wlsh <- length(unique(asv_wlsh$ASV_ID))
asv_total_dh <- length(unique(asv_dh$ASV_ID))

# =====================================================
# 一、基础数据
# =====================================================

cat("\n\n【一、数据基础】\n")
cat("──────────────────────────────────────────────────────────────────────\n")

cat(sprintf("  WLSH - 形态学总检出属数: %d\n", length(morpho_all_wlsh)))
cat(sprintf("  WLSH - eDNA总检出属数: %d\n", length(edna_all_wlsh)))
cat(sprintf("  WLSH - eDNA总ASV数: %d\n", asv_total_wlsh))
cat("\n")
cat(sprintf("  DH   - 形态学总检出属数: %d\n", length(morpho_all_dh)))
cat(sprintf("  DH   - eDNA总检出属数: %d\n", length(edna_all_dh)))
cat(sprintf("  DH   - eDNA总ASV数: %d\n", asv_total_dh))

# =====================================================
# 二、共有与特有类群分析
# =====================================================

cat("\n\n【二、共有与特有类群分析】\n")
cat("──────────────────────────────────────────────────────────────────────\n")

shared_wlsh <- intersect(morpho_all_wlsh, edna_all_wlsh)
shared_dh <- intersect(morpho_all_dh, edna_all_dh)

overlap_rate_wlsh <- length(shared_wlsh) / length(morpho_all_wlsh) * 100
overlap_rate_dh <- length(shared_dh) / length(morpho_all_dh) * 100

cat("\n【2.1 共有类群】\n")
cat(sprintf("  WLSH 共有属: %d 个 (重叠率 %.1f%%)\n", length(shared_wlsh), overlap_rate_wlsh))
if(length(shared_wlsh) > 0) {
  cat(sprintf("          %s\n", paste(shared_wlsh, collapse = ", ")))
}
cat(sprintf("  DH   共有属: %d 个 (重叠率 %.1f%%)\n", length(shared_dh), overlap_rate_dh))
if(length(shared_dh) > 0) {
  cat(sprintf("          %s\n", paste(shared_dh, collapse = ", ")))
}

# 形态学特有
morpho_only_wlsh <- setdiff(morpho_all_wlsh, edna_all_wlsh)
morpho_only_dh <- setdiff(morpho_all_dh, edna_all_dh)

cat("\n【2.2 形态学特有类群（仅形态学检出）】\n")
cat(sprintf("  WLSH 形态学特有属: %d 个\n", length(morpho_only_wlsh)))
if(length(morpho_only_wlsh) > 0) {
  cat(sprintf("          %s\n", paste(morpho_only_wlsh, collapse = ", ")))
}
cat(sprintf("\n  DH   形态学特有属: %d 个\n", length(morpho_only_dh)))
if(length(morpho_only_dh) > 0) {
  cat(sprintf("          %s\n", paste(morpho_only_dh, collapse = ", ")))
}

# eDNA特有
edna_only_wlsh <- setdiff(edna_all_wlsh, morpho_all_wlsh)
edna_only_dh <- setdiff(edna_all_dh, morpho_all_dh)

edna_only_asv_wlsh <- asv_wlsh %>%
  filter(!is.na(Genus), Genus %in% edna_only_wlsh) %>%
  distinct(ASV_ID, Genus)

edna_only_asv_dh <- asv_dh %>%
  filter(!is.na(Genus), Genus %in% edna_only_dh) %>%
  distinct(ASV_ID, Genus)

cat("\n【2.3 eDNA特有类群（仅eDNA检出）】\n")
cat(sprintf("  WLSH eDNA特有属: %d 个\n", length(edna_only_wlsh)))
if(length(edna_only_wlsh) > 0 && length(edna_only_wlsh) <= 20) {
  cat(sprintf("          %s\n", paste(edna_only_wlsh, collapse = ", ")))
} else if(length(edna_only_wlsh) > 20) {
  cat(sprintf("          (前20个): %s\n", paste(head(edna_only_wlsh, 20), collapse = ", ")))
}
cat(sprintf("        eDNA特有ASV: %d 条\n", nrow(edna_only_asv_wlsh)))

cat(sprintf("\n  DH   eDNA特有属: %d 个\n", length(edna_only_dh)))
if(length(edna_only_dh) > 0) {
  cat(sprintf("          %s\n", paste(edna_only_dh, collapse = ", ")))
}
cat(sprintf("        eDNA特有ASV: %d 条\n", nrow(edna_only_asv_dh)))

# =====================================================
# 三、汇总统计表
# =====================================================

cat("\n\n【三、汇总统计表】\n")
cat("──────────────────────────────────────────────────────────────────────\n\n")

cat("┌────────────────────────────────┬──────────┬──────────┐\n")
cat("│        统计指标               │   WLSH   │    DH    │\n")
cat("├────────────────────────────────┼──────────┼──────────┤\n")
cat(sprintf("│ 形态学总检出属数               │    %3d    │    %3d    │\n", length(morpho_all_wlsh), length(morpho_all_dh)))
cat(sprintf("│ eDNA总检出属数                 │    %3d    │    %3d    │\n", length(edna_all_wlsh), length(edna_all_dh)))
cat(sprintf("│ eDNA总ASV数                    │    %3d    │    %3d    │\n", asv_total_wlsh, asv_total_dh))
cat("├────────────────────────────────┼──────────┼──────────┤\n")
cat(sprintf("│ 共有属数量                     │    %3d    │    %3d    │\n", length(shared_wlsh), length(shared_dh)))
cat(sprintf("│ 重叠率（以形态学为基准）       │   %5.1f%%  │   %5.1f%%  │\n", overlap_rate_wlsh, overlap_rate_dh))
cat("├────────────────────────────────┼──────────┼──────────┤\n")
cat(sprintf("│ 形态学特有属                   │    %3d    │    %3d    │\n", length(morpho_only_wlsh), length(morpho_only_dh)))
cat(sprintf("│ eDNA特有属                     │    %3d    │    %3d    │\n", length(edna_only_wlsh), length(edna_only_dh)))
cat(sprintf("│ eDNA特有ASV                    │    %3d    │    %3d    │\n", nrow(edna_only_asv_wlsh), nrow(edna_only_asv_dh)))
cat("└────────────────────────────────┴──────────┴──────────┘\n")

# =====================================================
# 四、保存结果
# =====================================================

cat("\n\n【四、保存详细结果到CSV】\n")
cat("──────────────────────────────────────────────────────────────────────\n")

write.csv(data.frame(
  Lake = c(rep("WLSH", length(shared_wlsh)), rep("DH", length(shared_dh))),
  Genus = c(shared_wlsh, shared_dh)
), file.path(output_dir, "1_共有属列表.csv"), row.names = FALSE)

write.csv(data.frame(
  Lake = c(rep("WLSH", length(morpho_only_wlsh)), rep("DH", length(morpho_only_dh))),
  Genus = c(morpho_only_wlsh, morpho_only_dh)
), file.path(output_dir, "2_形态学特有属.csv"), row.names = FALSE)

edna_specific <- bind_rows(
  edna_only_asv_wlsh %>% mutate(Lake = "WLSH"),
  edna_only_asv_dh %>% mutate(Lake = "DH")
)
write.csv(edna_specific, file.path(output_dir, "3_eDNA特有属及ASV.csv"), row.names = FALSE)

summary_table <- data.frame(
  指标 = c("形态学总检出属数", "eDNA总检出属数", "eDNA总ASV数",
         "共有属数量", "重叠率(%)",
         "形态学特有属", "eDNA特有属", "eDNA特有ASV"),
  WLSH = c(length(morpho_all_wlsh), length(edna_all_wlsh), asv_total_wlsh,
           length(shared_wlsh), round(overlap_rate_wlsh, 1),
           length(morpho_only_wlsh), length(edna_only_wlsh), nrow(edna_only_asv_wlsh)),
  DH = c(length(morpho_all_dh), length(edna_all_dh), asv_total_dh,
         length(shared_dh), round(overlap_rate_dh, 1),
         length(morpho_only_dh), length(edna_only_dh), nrow(edna_only_asv_dh))
)
write.csv(summary_table, file.path(output_dir, "0_汇总统计表.csv"), row.names = FALSE)

cat("✓ 结果已保存至:", output_dir, "\n")

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                          报告生成完成                                ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")

# =====================================================
# 审稿意见 2.10 + 4.2 回应图
# 图1：属级别 Venn 图（两圆等大，重叠率标在交集数字下方）
# 图2：属级别 UpSet 图（季节 × 方法）
# UpSet 导出：svglite 导出 SVG（AI 可正常显示）+ PNG 预览
# 配色：Morphology = #D6857E，eDNA = #5A7D9A
# =====================================================

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(VennDiagram)
library(UpSetR)
library(grid)
library(futile.logger)
library(svglite)

# 路径
base_dir <- "E:/写东西、杂物/论文/Yan/测序文件整理"
data_dir <- file.path(base_dir, "数据")
output_dir <- file.path(base_dir, "输出图片/Reviewer_2.10_4.2")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

flog.threshold(ERROR, name = "VennDiagramLogger")

method_colors <- c("Morphology" = "#D6857E", "eDNA" = "#5A7D9A")
method_edges  <- c("Morphology" = "#A9524C", "eDNA" = "#3E5C75")

timepoints <- list(
  list(folder = "WLSH04", lake = "WLSH", season = "Spring"),
  list(folder = "WLSH07", lake = "WLSH", season = "Summer"),
  list(folder = "WLSH10", lake = "WLSH", season = "Autumn"),
  list(folder = "DH04", lake = "DH", season = "Spring"),
  list(folder = "DH07", lake = "DH", season = "Summer"),
  list(folder = "DH10", lake = "DH", season = "Autumn")
)

exclude_phyla <- c(
  "Dinoflagellata", "Peronosporomycetes", "Hyphochytriomycetes",
  "Cryptomycota", "Chytridiomycota", "Ascomycota", "Basidiomycota",
  "Aphelidea", "Mucoromycota", "Zoopagomycota", "Porphyridiophyceae",
  "Kathablepharidae", "Chlorophyta", "Ochrophyta", "Bacillariophyta",
  "Rhodophyta", "Cryptophyta", "Haptophyta", "Phragmoplastophyta",
  "Pavlovophyceae", "Prymnesiophyceae", "Labyrinthulomycetes",
  "MAST-2", "MAST-12", "MAST-6", "LKM15",
  "Arthropoda", "Rotifera", "Nematoda", "Annelida",
  "Streptophyta", "Magnoliophyta", "Pinophyta"
)
exclude_genera <- c(
  "Unclassified", "uncultured", "Pythium", "Phytophthora",
  "Porphyridium", "MAST-2", "A31", "Ceratium", "MAST-12C", "NA",
  "Symbiodinium", "Pfiesteria", "Woloszynskia", "Pelagodinium",
  "Durinskia", "Amphidinium", "Scrippsiella", "Biecheleria",
  "Katablepharis", "Ebria", "Glaucocystis"
)
exclude_families <- c(
  "MAST-2", "Thraustochytriaceae", "MAST-12A", "Incertae_Sedis", "MAST-12C",
  "Cryptomonadales", "Pavlovophyceae", "Chlorophyceae", "Trebouxiophyceae",
  "Mediophyceae", "Mamiellales", "Chlorodendrales", "Pedinellales",
  "Chromulinales", "Bacillariophyceae", "Fragilariales", "Chrysophyceae",
  "Ochromonadales", "Sphaeropleales", "Ulvophyceae", "Ochrophyta",
  "Pedinomonadales", "Mamiellophyceae", "Eustigmatales", "Prymnesiales",
  "Zygnematophyceae", "Thoracosphaeraceae", "Isochrysidales", "Chattonellales",
  "Gonyaulacales", "Woloszynskia", "Porphyridiophyceae",
  "Peronosporomycetes", "LKM11", "Malasseziaceae", "Pleosporaceae",
  "Didymellaceae", "Cladosporiaceae", "Aphelidea", "Trichomeriaceae",
  "Rhinosporideacae", "Rhizophydiaceae", "Schizophyllaceae", "Kappamycetaceae",
  "Ichthyophonidae", "Filobasidiaceae", "Chytridiaceae", "Aspergillaceae",
  "Gromochytriaceae", "Debaryomycetaceae", "Sporormiaceae", "Phaffomycetaceae",
  "Bulleribasidiaceae", "Pleurotaceae", "Hydropus", "Pythiaceae",
  "Magnoliophyta", "Pinophyta", "uncultured", "NA", "Calanoida"
)

# ========== 读取形态学数据 ==========
cat("========== 读取形态学数据 ==========\n")
morpho_all_data <- data.frame()

for(tp in timepoints) {
  file_path <- file.path(data_dir, tp$folder, "morpho.txt")
  if(!file.exists(file_path)) next
  
  morpho <- read.table(file_path, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE)
  genus_names <- gsub(" sp\\.?$", "", trimws(morpho[[1]]))
  morpho_matrix <- as.matrix(morpho[, -1])
  rownames(morpho_matrix) <- genus_names
  mode(morpho_matrix) <- "numeric"
  morpho_matrix[is.na(morpho_matrix)] <- 0
  morpho_matrix <- morpho_matrix[rowSums(morpho_matrix) > 0, , drop = FALSE]
  if(nrow(morpho_matrix) == 0) next
  
  morpho_agg <- aggregate(morpho_matrix, by = list(Genus = rownames(morpho_matrix)), FUN = sum)
  
  morpho_long <- morpho_agg %>%
    pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Count") %>%
    mutate(Lake = tp$lake, Season = tp$season,
           Sample_full = paste(tp$folder, Sample, sep = "_"))
  morpho_all_data <- rbind(morpho_all_data, morpho_long)
  cat("✓ 形态学:", tp$folder, "- 检出属:", nrow(morpho_agg), "\n")
}

# ========== 读取 eDNA 数据 ==========
cat("\n========== 读取 eDNA 数据 ==========\n")
asv_data_raw <- data.frame()

for(tp in timepoints) {
  asv_file <- file.path(data_dir, tp$folder, "ASV.txt")
  tax_file <- file.path(data_dir, tp$folder, "tax.txt")
  if(!file.exists(asv_file) || !file.exists(tax_file)) next
  
  asv <- read.table(asv_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  tax <- read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    na.strings = c("NA", "", "uncultured"))
  
  tax$Phylum <- gsub("^p__", "", tax$Phylum)
  tax$Class  <- gsub("^c__", "", tax$Class)
  tax$Order  <- gsub("^o__", "", tax$Order)
  tax$Family <- gsub("^f__", "", tax$Family)
  tax$Genus  <- gsub("^g__", "", tax$Genus)
  
  asv$ASV_ID <- as.character(asv$ASV_ID)
  tax$ASV_ID <- as.character(tax$ASV_ID)
  
  asv_counts <- if("taxonomy" %in% colnames(asv)) asv[, colnames(asv) != "taxonomy"] else asv
  
  asv_long <- asv_counts %>%
    pivot_longer(cols = -ASV_ID, names_to = "Sample", values_to = "Count") %>%
    filter(Count > 0) %>%
    left_join(tax, by = "ASV_ID") %>%
    filter(!(Phylum %in% exclude_phyla)) %>%
    filter(!is.na(Phylum)) %>%
    filter(!(Genus %in% exclude_genera)) %>%
    filter(!(Family %in% exclude_families)) %>%
    mutate(Lake = tp$lake, Season = tp$season,
           Sample_full = paste(tp$folder, Sample, sep = "_"))
  
  asv_data_raw <- rbind(asv_data_raw, asv_long)
  cat("✓ eDNA:", tp$folder, "- 保留ASV:", length(unique(asv_long$ASV_ID)), "\n")
}

# ========== 按湖泊汇总属列表 ==========
morpho_genus_by_lake <- morpho_all_data %>%
  filter(Count > 0) %>%
  group_by(Lake, Genus) %>%
  summarise(Total = sum(Count), .groups = "drop")

edna_genus_by_lake <- asv_data_raw %>%
  filter(!is.na(Genus), Genus != "", Genus != "Unclassified", Genus != "uncultured") %>%
  group_by(Lake, Genus) %>%
  summarise(Total = sum(Count), .groups = "drop")

# ========== 计算重叠率 ==========
overlap_stats <- data.frame()

for(lk in c("WLSH", "DH")) {
  m <- morpho_genus_by_lake %>% filter(Lake == lk) %>% pull(Genus) %>% unique()
  e <- edna_genus_by_lake %>% filter(Lake == lk) %>% pull(Genus) %>% unique()
  shared <- intersect(m, e)
  overlap_rate <- ifelse(length(m) > 0, length(shared) / length(m) * 100, NA)
  overlap_stats <- rbind(overlap_stats, data.frame(
    Lake = lk,
    Morphology = length(m),
    eDNA = length(e),
    Shared = length(shared),
    Overlap_rate = round(overlap_rate, 1)
  ))
}
write.csv(overlap_stats, file.path(output_dir, "Overlap_Statistics_Genus.csv"), row.names = FALSE)
cat("\n重叠率统计:\n")
print(overlap_stats)

# =====================================================
# 图1：属级别 Venn 图
# =====================================================
cat("\n========== 绘制 Venn 图 ==========\n")

plot_venn_lake <- function(lake_name) {
  m <- morpho_genus_by_lake %>% filter(Lake == lake_name) %>% pull(Genus) %>% unique()
  e <- edna_genus_by_lake %>% filter(Lake == lake_name) %>% pull(Genus) %>% unique()
  
  st <- overlap_stats %>% filter(Lake == lake_name)
  fill_cols <- unname(method_colors)
  edge_cols <- unname(method_edges)
  
  venn_plot <- venn.diagram(
    x = list("Morphology" = m, "eDNA" = e),
    category.names = c("Morphology", "eDNA"),
    filename = NULL,
    output = TRUE,
    scaled = FALSE,
    lwd = 3,
    lty = "solid",
    col = edge_cols,
    fill = fill_cols,
    alpha = 0.45,
    cex = 1.8, fontface = "bold",
    cat.cex = 1.5, cat.fontface = "bold", cat.col = edge_cols,
    main = lake_name,
    main.cex = 1.4, main.fontface = "bold"
  )
  
  png(file.path(output_dir, paste0("Venn_", lake_name, "_Genus.png")),
      width = 8, height = 8, units = "in", res = 300, bg = "white")
  grid.draw(venn_plot)
  grid.text(paste0("Overlap = ", st$Overlap_rate, "%"),
            x = 0.5, y = 0.47,
            gp = gpar(fontsize = 12, fontface = "bold", col = "grey20"))
  dev.off()
  
  pdf(file.path(output_dir, paste0("Venn_", lake_name, "_Genus.pdf")),
      width = 8, height = 8, useDingbats = FALSE)
  grid.draw(venn_plot)
  grid.text(paste0("Overlap = ", st$Overlap_rate, "%"),
            x = 0.5, y = 0.47,
            gp = gpar(fontsize = 12, fontface = "bold", col = "grey20"))
  dev.off()
  
  cat("✓ Venn 已保存:", lake_name,
      "| Shared =", st$Shared,
      "| Overlap =", st$Overlap_rate, "%\n")
}

plot_venn_lake("WLSH")
plot_venn_lake("DH")

# =====================================================
# 图2：属级别 UpSet 图（季节 × 方法）
# 导出：svglite 导出 SVG（AI 可直接打开、逐对象编辑）+ PNG 预览
# =====================================================
cat("\n========== 绘制 UpSet 图 ==========\n")

build_upset_matrix <- function(lake_name) {
  m_df <- morpho_all_data %>%
    filter(Lake == lake_name, Count > 0) %>%
    mutate(Set = paste(Season, "Morphology", sep = "_")) %>%
    select(Genus, Set) %>% distinct() %>% mutate(Value = 1)
  
  e_df <- asv_data_raw %>%
    filter(Lake == lake_name, !is.na(Genus), Genus != "",
           Genus != "Unclassified", Genus != "uncultured", Count > 0) %>%
    mutate(Set = paste(Season, "eDNA", sep = "_")) %>%
    select(Genus, Set) %>% distinct() %>% mutate(Value = 1)
  
  all_df <- bind_rows(m_df, e_df)
  if(nrow(all_df) == 0) return(NULL)
  
  mat <- all_df %>%
    pivot_wider(names_from = Set, values_from = Value, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- mat$Genus
  mat$Genus <- NULL
  return(mat)
}

save_upset <- function(lake_name) {
  mat <- build_upset_matrix(lake_name)
  if(is.null(mat) || nrow(mat) == 0) {
    cat("⚠️", lake_name, "UpSet 矩阵为空，跳过\n")
    return(NULL)
  }
  
  set_order <- c("Spring_Morphology", "Summer_Morphology", "Autumn_Morphology",
                 "Spring_eDNA", "Summer_eDNA", "Autumn_eDNA")
  set_order <- intersect(set_order, colnames(mat))
  mat <- mat[, set_order, drop = FALSE]
  
  bar_colors <- ifelse(grepl("eDNA", set_order),
                       method_colors["eDNA"],
                       method_colors["Morphology"])
  
  # ★ SVG 导出（AI 可正常显示、逐对象可编辑）
  svglite(file.path(output_dir, paste0("UpSet_", lake_name, "_Genus.svg")),
          width = 9, height = 5.5)
  up <- upset(mat,
              nsets = ncol(mat),
              nintersects = 30,
              sets.bar.color = bar_colors,
              order.by = "freq",
              main.bar.color = "grey30",
              matrix.color = "grey30",
              mainbar.y.label = "Shared genera",
              sets.x.label = "Genera per set",
              text.scale = c(1.2, 1.1, 1.0, 1.0, 1.1, 1.0))
  print(up)
  dev.off()
  
  # PNG 预览（白底）
  png(file.path(output_dir, paste0("UpSet_", lake_name, "_Genus.png")),
      width = 9, height = 5.5, units = "in", res = 300, bg = "white", type = "cairo")
  print(up)
  dev.off()
  
  cat("✓ UpSet 已保存 (SVG + PNG):", lake_name, "\n")
}

save_upset("WLSH")
save_upset("DH")

cat("\n✅ 全部完成！输出目录:", output_dir, "\n")
cat("  Venn: Venn_WLSH_Genus / Venn_DH_Genus (.pdf/.png)\n")
cat("  UpSet: UpSet_WLSH_Genus / UpSet_DH_Genus (.svg/.png)\n")
cat("  提示：UpSet 的 SVG 直接在 Illustrator 中打开即可编辑；如需 PDF，在 AI 里另存为 PDF。\n")
