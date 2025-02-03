library(ggplot2)
library(stringr)
library(reshape2)

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
  New_TB$V2 <- factor(New_TB$V2, 
    levels = unique(New_TB$V2[order(log(New_TB$Count), decreasing = T)]))
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

ggplot(HC[HC$Count>1,], aes(X1, ratio, fill = X2)) +
    geom_bar(stat = 'identity', position = 'dodge')+ facet_wrap(~V2, scales = 'free_y', ncol = 15) +
    theme_bw() 
ggsave('plot/HV_Bar_split_enrichment.png', w = 20, h = 11)
ggsave('plot/HV_Bar_split_enrichment.svg', w = 20, h = 11)


ggplot(LC[LC$Count>1,], aes(X1, ratio, fill = X2)) +
    geom_bar(stat = 'identity', position = 'dodge')+ facet_wrap(~V2, scales = 'free_y', ncol = 15) +
    theme_bw() 
ggsave('plot/LV_Bar_split_enrichment.png', w = 20, h = 11)
ggsave('plot/LV_Bar_split_enrichment.svg', w = 20, h = 11)

ggplot(KC[KC$Count>1,], aes(X1, ratio, fill = X2)) +
    geom_bar(stat = 'identity', position = 'dodge')+ facet_wrap(~V2, scales = 'free_y', ncol = 15) +
    theme_bw() 
ggsave('plot/KV_Bar_split_enrichment.png', w = 20, h = 11)
ggsave('plot/KV_Bar_split_enrichment.svg', w = 20, h = 11)
