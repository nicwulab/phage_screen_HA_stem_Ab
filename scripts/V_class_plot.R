library(ggplot2)
library(stringr) 

TB <- read.csv('Result/V_class_family_count.tsv', sep ='', header = F)
TB <- cbind(TB, data.frame(str_split_fixed(TB$V1, "_", 2)))
#TB = TB[!TB$X2=='id',]
TB <- na.omit(TB)

Counts <- read.table('Result/Counts_total.csv', header = 0)
Counts <- cbind(Counts, data.frame(str_split_fixed(Counts$V2, "_", 2)))

TB['ratio'] = 0
for(sample in unique(TB$X1)){
    TB$ratio[TB$X1==sample]  = TB$V3[TB$X1==sample]/Counts$V3[Counts$X1==sample]
}

ggplot(TB[TB$V3>1,], aes(X1, ratio, fill = X2)) +
    geom_bar(stat = 'identity', position = 'dodge')+ facet_wrap(~V2, scales = 'free_y', ncol = 15) +
    theme_bw() 

ggsave('Picture/Family_all_bar.png', w = 20, h = 11)
