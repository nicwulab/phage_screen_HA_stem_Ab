library(ggplot2)

TB <- read.table("Result/Paired_aa_counts.tsv", header = TRUE, sep = "\t")

tb <- head(TB, 10)

tb$Name <- paste('Top', row.names(tb), sep = "_")
tb <- rbind(tb, data.frame(Count = sum(TB[-c(1:10),]), Name = "Other")) 

tb$Name <- factor(tb$Name, levels = tb$Name)


# pi plot
p <- ggplot(tb, aes(x = "", y = Count, fill = Name)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Top 10 Paired Amino Acid Counts") +
  scale_fill_brewer(palette = "Set3")


ggsave("Picture/Paired_aa_counts.svg", plot = p, width = 5.72, height = 4.22)
