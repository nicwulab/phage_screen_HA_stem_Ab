library(ggplot2)
library(stringr)
library(reshape2)
library(patchwork)
library(plotly)

ChainPair_tb <- read.table("Result/IGH_IGL.tsv", header = TRUE, sep = "\t")
ChainPair_tb <- cbind(ChainPair_tb, data.frame(str_split_fixed(ChainPair_tb$X, ":", 2)))
ChainPair_tb <- ChainPair_tb[, -1]
colnames(ChainPair_tb) <- c("Count", "P", 'Pair')
ChainPair_tb[-grep("/", ChainPair_tb$Pair),] -> ChainPair_tb
ChainPair_tb <- reshape(ChainPair_tb, idvar = "Pair", timevar = "P", direction = "wide")

colnames(ChainPair_tb) <- str_remove(colnames(ChainPair_tb), "Count.")
ChainPair_tb <- ChainPair_tb[-3]

# assign NA to 0
ChainPair_tb[is.na(ChainPair_tb)] <- 0
# sum the last two columns and sorted by them
ChainPair_tb <- ChainPair_tb[order(rowSums(ChainPair_tb[c('P3_cm', 'P3_wt')]), decreasing = TRUE),]
# Finally, split the Pair column into IGH and IGL
ChainPair_tb <- cbind(ChainPair_tb, data.frame(str_split_fixed(ChainPair_tb$Pair, ":", 2)))
colnames(ChainPair_tb)[7:8] <- c("IGH", "IGL")

ChainPair_tb <- ChainPair_tb[-c(3:4)]
ChainPair_tb <- ChainPair_tb[-1]
# calculate the ratio
ratio_tb <- ChainPair_tb[1:3]/colSums(ChainPair_tb[1:3])
paste(colnames(ratio_tb), "_Ratio", sep='') -> colnames(ratio_tb)
cbind(ChainPair_tb, ratio_tb) -> ChainPair_tb

ChainPair_tb$P3_cmR <- ChainPair_tb$P3_cm_Ratio / ChainPair_tb$P0_Ab_Ratio
ChainPair_tb$P3_wtR <- ChainPair_tb$P3_wt_Ratio / ChainPair_tb$P0_Ab_Ratio

# filter the top 20
head(unique(ChainPair_tb$IGH),20) -> topIGH
head(unique(ChainPair_tb$IGL),20) -> topIGL

ChainPairFilter_tb  <- ChainPair_tb[ChainPair_tb$IGH %in% topIGH & ChainPair_tb$IGL %in% topIGL,]

write.csv(ChainPair_tb, 'Result/IGH_IGL_plot.csv', row.names = FALSE)

p1 <- ggplot(ChainPairFilter_tb, aes(x = IGL, y = IGH)) +
  geom_point(aes(size = P3_cm)) +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "IGH-IGL Pairing", x = "IGL", y = "IGH") 

p2 <- ggplot(ChainPairFilter_tb, aes(x = IGL, y = IGH)) +
  geom_point(aes(size = P3_wt)) +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "IGH-IGL Pairing", x = "IGL", y = "IGH") 

p1/p2
ggsave("Picture/IGH_IGL.svg", width = 10.9, height = 8.44)


p1 <- ggplot(ChainPairFilter_tb, aes(x = IGL, y = IGH)) +
  geom_point(aes(size = P0_Ab)) +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "P0: IGH-IGL Pairing", x = "IGL", y = "IGH")

p2 <- ggplot(ChainPairFilter_tb, aes(x = IGL, y = IGH)) +
  geom_point(aes(size = P3_wt)) +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "P3: IGH-IGL Pairing", x = "IGL", y = "IGH")


ggsave("Picture/P0_IGH_IGL.svg", p1/p2, width = 10.9, height = 8.44)

