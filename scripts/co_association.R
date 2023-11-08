library(reshape2)
library(ggplot2)
library(stringr)

TB <- read.csv("Result/ExpRe.csv") 
TB = TB[TB$exp> 1.2,]
TB_l <- reshape(TB[-1], timevar = 'Type', idvar = 'ID', direction = 'wide')
TB_l <- na.omit(TB_l)
TB_l$dif <- abs(TB_l$exp.wt - TB_l$exp.cm)
TB_l$dif_rf <- abs(1-(TB_l$exp.wt/ TB_l$exp.cm))

Count_ID <- read.table("Result/Counts_duplciate.csv", header = 0)
Count <- read.table('Result/Counts_total.csv', header = 0)
Count_ID$ratio = 0
for(i in Count$V2){
    Count_ID$ratio[Count_ID$V2==i] = Count_ID$V3[Count_ID$V2==i]/Count$V3[Count$V2==i]
}
colnames(Count_ID) <- c('ID', 'Sample', 'Counts', 'ratio')
Count_ID <- cbind(Count_ID, as.data.frame(str_split_fixed(Count_ID$Sample, '_', 2)))
colnames(Count_ID)[5:6] <- c('Stage', 'Type')

Count_plot <- Count_ID[Count_ID$ID %in% head(TB_l$ID[order(TB_l$dif_rf)], 30),]
Count_plot$ID <- factor(Count_plot$ID, levels=TB_l$ID[order(TB_l$dif_rf)])

ggplot(Count_plot, aes(Stage, ratio, fill = Type)) +  
    geom_bar(stat = 'identity', position = 'dodge') + facet_wrap(~ID, scales = 'free_y')+ 
    theme_linedraw() + ggtitle('Exponential growing association') + theme(plot.title = element_text(hjust = .5))


