library(ggplot2)

TB <- read.csv('Result/Ratio_matrix.csv') 

#TB2 <- head(TB[order(TB$ratio, decreasing = T),], 100) 

TB2 <- TB[TB$ID<50,]

ggplot(TB2, aes(x = Round, y = ratio, fill = Type)) + geom_bar(stat = 'identity', position = 'dodge') + facet_wrap(~ID, scales = 'free_y')+ 
    theme_linedraw()

ggsave('Picture/First_100_bar.png', w = 20, h = 11)






