# AA_counts_all_density

library(ggplot2)
TB <- read.csv('Result/Du_AA_counts.csv')
TB[is.na(TB)] <- 0 

ggplot(data = TB) + geom_density(aes(x = Counts)) + theme_bw() + ggtitle('AA Counts for all 3 stages') + theme(plot.title = element_text(hjust = .5))

ggsave('Picture/Stat/AA_counts_all_density.png', w = 3.86, h =2.09)

TB$P23 = TB$P2_cm + TB$P2_wt + TB$P3_cm + TB$P3_wt
TB$P23[TB$P23>=10]=10
TB$P23[TB$P23<10]=0

ggplot(data = TB, aes(x = 'All', fill = as.factor(P23))) + geom_bar()+geom_text( 
        aes(label=paste( 100*signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)], digits=4), "%")), stat="count", 
        position=position_stack(vjust=0.5)) +  guides(fill=guide_legend(title="Threshold")) + 
        labs(y="Proportion", x = '') + theme_bw() + ggtitle('AA Counts\nP2 and P3 (wt & cm)') + theme(plot.title = element_text(hjust = .5)) 
ggsave('Picture/Stat/AA_counts_all_bar_ther10.png', w = 2.71, h = 5.77)

TB$P23 = TB$P2_cm + TB$P2_wt + TB$P3_cm + TB$P3_wt
TB$P23[TB$P23>=5]=5
TB$P23[TB$P23<5]=0

ggplot(data = TB, aes(x = 'All', fill = as.factor(P23))) + geom_bar()+geom_text( 
        aes(label=paste( 100*signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)], digits=4), "%")), stat="count", 
        position=position_stack(vjust=0.5)) +  guides(fill=guide_legend(title="Threshold")) + 
        labs(y="Proportion", x = '') + theme_bw() + ggtitle('AA Counts\nP2 and P3 (wt & cm)') + theme(plot.title = element_text(hjust = .5)) 
ggsave('Picture/Stat/AA_counts_all_bar_ther5.png', w = 2.71, h = 5.77)

   
