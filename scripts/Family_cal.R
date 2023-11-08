library(stringr)
library(reshape2)
library(ggplot2)

args<-commandArgs(TRUE)

Family <- args[1]
Top_N <- args[2]
OUT <- args[3]

TB = read.table(paste('Result/', Family,  '/IG', Family, '_ReadsID_counts.csv', sep =''), header = 0)
Count <- read.table('Result/Counts_total.csv', header = 0)
for(i in Count$V2){
    TB$V1[TB$V3==i] = TB$V1[TB$V3==i]/Count$V3[Count$V2==i]
}

TB$V4 <- as.data.frame(str_split_fixed(TB$V3, "_", 2))[[1]]

TB_m <- reshape(TB[-3], idvar =  "V2",timevar  = 'V4', direction =  'wid' )
colnames(TB_m) <- c('ID', 'P3', "P2", "P1")

TB_m[is.na(TB_m)] <- 0
TB_m$Sum <- apply(TB_m[-1], 1, sum)
TB_m <- TB_m[order(TB_m$Sum, decreasing = T),]

TB$V2 <- factor(TB$V2, levels = TB_m$ID)
colnames(TB) <- c('ratio', 'ID', 'Sample', "Stage")

TMP <- TB[TB$ID %in% head(TB_m$ID, Top_N),]
TMP$ID <- factor(TMP$ID, levels = unique(TMP$ID))

ggplot(TMP, aes(Stage, ratio)) + geom_bar(stat = 'identity') + theme_bw()+ 
    facet_wrap(~ID, scales = 'free_y', ncol = 13)+ ggtitle(paste( Family, 'Reads enrichment (Top', Top_N, 'according to the total ratio)'))+  
    theme( plot.title = element_text(hjust = .5))
ggsave(paste("Picture/IG", Family, "_Alltop", Top_N, ".png", sep=''), w = 20, h = 11)
write.table(head(TB_m$ID, Top_N), paste(OUT, "/IG", Family, "_Alltop", Top_N, ".list", sep = ''),  row.names = F, col.names = F)