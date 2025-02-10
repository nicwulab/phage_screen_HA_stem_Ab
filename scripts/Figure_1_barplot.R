library(ggplot2)
library(stringr)
library(reshape2)
library(patchwork)

theme_text <- theme(axis.text = element_text(family = "Arial", face = "bold", size = 7),
        strip.text = element_text(family = "Arial", face = "bold", size = 7),
        axis.title = element_text(family = "Arial", face = "bold", size = 7),
        legend.text = element_text(family = "Arial", face = "bold", size = 7),
        legend.title = element_text(family = "Arial", face = "bold", size = 7),
        plot.title = element_text(family = "Arial", face = "bold", size = 7, hjust = .5))

TB <- read.csv('Result/V_class_family_count.tsv', sep ='', header = F)
TB <- cbind(TB, data.frame(str_split_fixed(TB$V1, "_", 2)))

TB$Count <- TB$V3
TB$Count[is.na(TB$V3)] <- as.numeric(TB$V2[is.na(TB$V3)])
TB$V2[is.na(TB$V3)] <- "Other"

#TB = TB[!TB$X2=='id',]
#TB <- na.omit(TB)
TB = TB[TB$X1!='sequence',]
Counts <- read.table('Result/Counts_total.csv', header = 0)
Counts <- cbind(Counts, data.frame(str_split_fixed(Counts$V2, "_", 2)))

TB['ratio'] = 0
for(sample in unique(TB$X1)){
    TB$ratio[TB$X1==sample]  = TB$Count[TB$X1==sample]/Counts$V3[Counts$X1==sample]
}

HC <- TB[TB$V2 %in%  unique(TB$V2)[grep("HV", unique(TB$V2))],]
LC <- TB[TB$V2 %in%  unique(TB$V2)[grep("LV", unique(TB$V2))],]
KC <- TB[TB$V2 %in%  unique(TB$V2)[grep("KV", unique(TB$V2))],]

Enrich_cal <- function(HC){
  New_TB <- data.frame() 
  for(id in HC$V2[HC$X1=="P0"]){
    tmp = HC[HC$V2 == id,]
    ct = tmp[tmp$X1 == 'P0',]
    tmp = tmp[tmp$X1 != 'P0',]
    tmp$Count <- tmp$Count/ ct$Count
    New_TB <- rbind(New_TB, tmp)
  }
  New_TB$V2 <- paste("IG", New_TB$V2, sep = "")
  New_TB$V2 <- factor(New_TB$V2, 
    levels = sort(unique(New_TB$V2)))
  return(New_TB)
}

ggplot(Enrich_cal(HC), aes(V2, log(Count), fill = V1)) + geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 270, hjust =0, vjust =.5))
ggsave("plot/HC_bar_enrichment.png", w = 8.77, h = 3.78)
ggsave("plot/HC_bar_enrichment.svg", w = 8.77, h = 3.78)
ggplot(Enrich_cal(LC), aes(V2, log(Count), fill = V1, alpha = ratio)) + geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 270, hjust =0, vjust =.5))
ggsave("plot/LC_bar_enrichment.png", w = 8.77, h = 3.78)
ggsave("plot/LC_bar_enrichment.svg", w = 8.77, h = 3.78)
ggplot(Enrich_cal(KC), aes(V2, log(Count), fill = V1, alpha = ratio)) + geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 270, hjust =0, vjust =.5))
ggsave("plot/KC_bar_enrichment.png", w = 8.77, h = 3.78)
ggsave("plot/KC_bar_enrichment.svg", w = 8.77, h = 3.78)

ggplot(Enrich_cal(HC), aes(V2, log(Count), color = V1, size = ratio)) + 
  geom_point() + 
  theme_bw() + theme(axis.text.x = element_text(angle = 270, hjust =0, vjust =.5))
ggsave("plot/HC_buble_enrichment.png", w = 8.77, h = 3.78)
ggsave("plot/HC_buble_enrichment.svg", w = 8.77, h = 3.78)
ggplot(Enrich_cal(LC), aes(V2, log(Count), color = V1, size = ratio)) + 
  geom_point() + 
  theme_bw() + theme(axis.text.x = element_text(angle = 270, hjust =0, vjust =.5))
ggsave("plot/LC_buble_enrichment.png", w = 8.77, h = 3.78)
ggsave("plot/LC_buble_enrichment.svg", w = 8.77, h = 3.78)
ggplot(Enrich_cal(KC), aes(V2, log(Count), color = V1, size = ratio)) + 
  geom_point() + 
  theme_bw() + theme(axis.text.x = element_text(angle = 270, hjust =0, vjust =.5))
ggsave("plot/KC_buble_enrichment.png", w = 8.77, h = 3.78)
ggsave("plot/KC_buble_enrichment.svg", w = 8.77, h = 3.78)

# now, because there are too many entry in it, we want to do the filtering and show the most important data


TopCut <- function(HC, Thr = 25){
  tb.r <- reshape(HC[c('V1', 'V2', 'ratio')], timevar = 'V1', idvar = c('V2'), direction = 'wide')
  tb.r[is.na(tb.r)] <- 0
  colnames(tb.r) <- str_remove(colnames(tb.r), "ratio.")
  return(head(tb.r$V2[order(tb.r$P3_wt, decreasing = T)], Thr))
}
P1 <- ggplot(Enrich_cal(HC[HC$V2 %in% TopCut(HC),]), aes(V2, log(Count+1), color = V1, size = ratio)) + 
  geom_point() + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =.5),
    axis.title.x = element_blank()) 
P2 <- ggplot(Enrich_cal(LC[LC$V2 %in% TopCut(LC),]), aes(V2, log(Count+1), color = V1, size = ratio)) + 
  geom_point() + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =.5),
    axis.title.x = element_blank())
P3 <- ggplot(Enrich_cal(KC[KC$V2 %in% TopCut(KC),]), aes(V2, log(Count+1), color = V1, size = ratio)) + 
  geom_point() + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =.5))
P1/P2/P3

ggsave("plot/HC_buble_enrichment_top.png", w= 8.2, h = 7.76)
ggsave("plot/HC_buble_enrichment_top.svg", w= 8.2, h = 7.76)

