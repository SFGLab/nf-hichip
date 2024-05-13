setwd("~/Desktop/Figure.3_source/Fig_3A/.")
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)

#Data
file = read.csv(file = "stripes.csv")
colnames(file)[1] = "Sample"

file2 = pivot_longer(file,cols = c("Stripenn", "gStripes_MAPS","gStripes_HiCCUPS","gStripes_ChIA_PIPE"),
                     names_to = "Algorithm",
                     values_to = "Value") 

file2$Sample <- factor(file2$Sample , levels = c("HCC1599 SMC1 HiChIP","Rec1 SMC1 HiChIP", "GM12878 SMC1 HiChIP","HG00731 SMC1 dcHiChIP"))
#file2$Algorithm <- factor(file2$Algorithm , levels = c("Stripenn","gStripes_MAPS","gStripes_HiCCUPS","gStripes_ChIA_PIPE") , ordered = TRUE)
trt = file2$Algorithm
file2 %>% mutate(Algorithm = fct_rev(Algorithm))

barchart = ggplot()
barchart = barchart + geom_col(data = file2, aes(x = Value, y = Sample, fill = trt), position = position_dodge(.6), width = 0.50) 
barchart = barchart + scale_x_continuous(expand = c(0, 0, .05, 0)) +labs(x="", y="", fill ="") 
barchart = barchart + theme_bw() + scale_fill_brewer(palette="Spectral") 
barchart = barchart + theme(legend.position="top", axis.text.y = element_text(size = 10) , text=element_text(family="Arial")) 
barchart


ggsave("stripes.png", plot = barchart,dpi = 300)
ggsave("stripes.svg", plot = barchart)

barchart = ggplot()
barchart = barchart + geom_col(data = file2, aes(x = Value, y = Sample, fill = trt), position = position_dodge(.6), width = 0.50) 
barchart = barchart + scale_x_continuous(expand = c(0, 0, .05, 0)) +labs(x="", y="", fill ="") 
barchart = barchart + theme_bw() + scale_fill_brewer(palette="Spectral",
                                                     labels = c("                                        ",
                                                                "                                        ",
                                                                "                                        ",
                                                                "                                        ")) 
barchart = barchart + theme(legend.position="top",axis.text.x=element_blank(),axis.text.y=element_blank())
barchart
ggsave("stripes_withoutL&L.png", plot = barchart,dpi = 300)
ggsave("stripes_withoutL&L.svg", plot = barchart)