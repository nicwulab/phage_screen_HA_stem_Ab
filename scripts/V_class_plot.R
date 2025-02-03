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

ggplot(TB[TB$Count>1,], aes(X1, ratio, fill = X2)) +
    geom_bar(stat = 'identity', position = 'dodge')+ facet_wrap(~V2, scales = 'free_y', ncol = 15) +
    theme_bw() 

ggsave('Picture/Family_all_bar.png', w = 20, h = 11)


##

TB_m <- reshape(TB[c(1,2,7)], timevar = 'V1', idvar = 'V2', direction = 'wide')
TB_m[is.na(TB_m)] <- 0

row.names(TB_m) <- TB_m[[1]]
TB_m <- TB_m[-1]
TB_m$P0 = TB_m[[1]]
TB_m$P2 = apply(TB_m[2:3], 1, mean)
TB_m$P3 = apply(TB_m[4:5], 1, mean)


TB_m$Co <- NA
for(i in c(1:nrow(TB_m))){
    tmp = data.frame(t(TB_m[i,c("P0", "P2", "P3")]))
    colnames(tmp) <- 'Y'
    tmp$X = c(0,2,3)
    TB_m$Co[i] <- lm(Y~X, tmp)$coefficients[[2]]
}

























