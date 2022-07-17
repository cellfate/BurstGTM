# import R package
library("ggplot2")
library("ggpubr")

rm(list = ls())
# read result file
genename <- read.table("data/genename.txt", sep = '\t', header = FALSE)
result.GM <- read.table("results/results_MEF.csv", header = FALSE, sep = ',')
result.GM <- cbind(genename, result.GM)
colnames(result.GM) <- c('genename', 'bf_GM', 'bs_GM')

data.all <- read.table("data/MEF_QC_all.csv", sep = ',', header = FALSE)
data.all[data.all == -1] <- NA
result.GM$mean_data <- apply(data.all,1,function(x) {mean(x, na.rm = T)})

result.ONOFF <- read.table("results/results_onoff_MEF.csv", header = FALSE, sep = ',')
result.ONOFF <- cbind(genename, result.ONOFF)
colnames(result.ONOFF) <- c('genename', 'kon', 'koff', 'mu', 'fval', 'bf_ONOFF', 'bs_ONOFF', 'mean_ONOFF', 'noise')

result.ONOFF$bf_ONOFF = 1/(1/result.ONOFF$kon + 1/result.ONOFF$koff)
result.ONOFF$bf_ONOFF2 = result.ONOFF$kon
result <- cbind(result.GM,result.ONOFF[,c('bf_ONOFF', 'bf_ONOFF2', 'bs_ONOFF','mean_ONOFF')])

p_A <- ggplot(result, aes(x = log10(bf_ONOFF), y = log10(bs_ONOFF), colour = log10(mean_data))) + 
  geom_point(size = 0.7)+
  scale_colour_gradient(low = "#2196f3", high = "#f44336") +
  labs(y = expression(log[10](bs)), x = expression(log[10](bf))) +
  xlim(-1.3,0.4) +
  ylim(0,2) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size = 0.25, lineend = 10),
        panel.grid = element_blank())
p_A

p_B <- ggplot(result, aes(x = log10(bf_GM),y = log10(bf_ONOFF))) +
  geom_point(size = 0.05, color = '#388ECD')+
  geom_abline(slope = 1, color = '#E73118', linetype = "dashed") +
  stat_cor(method = "pearson", size = 2) + 
  labs(x = "GM", y = "Telegraph model", title = expression(log[10](bf))) + 
  theme_bw() + 
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
p_B

p_C <- ggplot(result, aes(x = log10(bs_GM),y = log10(bs_ONOFF))) +
  geom_point(size = 0.05, color = '#388ECD')+
  geom_abline(slope = 1, color = '#E73118',linetype = "dashed") +
  stat_cor(method = "pearson", size = 2) + 
  labs(x = "GM", y = "Telegraph model", title = expression(log[10](bf))) + 
  theme_bw() + 
  theme(title = element_text(colour = 'black', size = 6),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6, color="black"),
        legend.key.size = unit(15, "pt"),
        axis.title = element_text(colour = 'black', size = 8),
        axis.text = element_text(colour = 'black', size = 6),
        axis.ticks = element_line(size=0.25, lineend = 10),
        panel.grid = element_blank()) 
p_C

pic1 <- ggarrange(p_A, p_B, p_C, ncol = 3, nrow = 1, widths=c(4), heights=c(4), align = "hv")
pic1

ggsave('Figure/Figure5.pdf', width = 5.1, height = 1.8, useDingbats = FALSE)


# rm(list = ls())
# # read result file
# genename <- read.table("data/genename.txt", sep = '\t', header = FALSE)
# result.GM <- read.table("results/results_MEF_bf2.csv", header = FALSE, sep = ',')
# result.GM <- cbind(genename, result.GM)
# colnames(result.GM) <- c('genename', 'bf_GM', 'bs_GM')
# 
# data.all <- read.table("data/MEF_QC_all.csv", sep = ',', header = FALSE)
# data.all[data.all == -1] <- NA
# result.GM$mean_data <- apply(data.all,1,function(x) {mean(x, na.rm = T)})
# 
# result.ONOFF <- read.table("results/results_onoff_MEF.csv", header = FALSE, sep = ',')
# result.ONOFF <- cbind(genename, result.ONOFF)
# colnames(result.ONOFF) <- c('genename', 'kon', 'koff', 'mu', 'fval', 'bf_ONOFF', 'bs_ONOFF', 'mean_ONOFF', 'noise')
# 
# result.ONOFF$bf_ONOFF = 1/(1/result.ONOFF$kon + 1/result.ONOFF$koff)
# result.ONOFF$bf_ONOFF2 = result.ONOFF$kon
# result <- cbind(result.GM,result.ONOFF[,c('bf_ONOFF', 'bf_ONOFF2', 'bs_ONOFF','mean_ONOFF')])
# 
# p_SA <- ggplot(result, aes(x = log10(bf_GM),y = log10(bf_ONOFF2))) +
#   geom_point(size = 0.05, color = '#388ECD')+
#   geom_abline(slope = 1, color = '#E73118', linetype = "dashed") +
#   stat_cor(method = "pearson", size = 2) + 
#   labs(x = "GM", y = "Telegraph model", title = expression(log[10](bf))) + 
#   theme_bw() + 
#   xlim(-2,0) + 
#   theme(title = element_text(colour = 'black', size = 6),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 6, color="black"),
#         legend.key.size = unit(15, "pt"),
#         axis.title = element_text(colour = 'black', size = 8),
#         axis.text = element_text(colour = 'black', size = 6),
#         axis.ticks = element_line(size=0.25, lineend = 10),
#         panel.grid = element_blank()) 
# p_SA
# 
# ggsave('../Figure/FigureS3.pdf', width = 1.7, height = 1.8, useDingbats = FALSE)