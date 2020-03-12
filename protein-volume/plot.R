install.packages("tidyverse")
library(ggplot2)
final = read.table("/Users/zhouwang/Documents/NOW/ECS 129/protein-volume/protein-volume/V_N.txt")
ggplot(final, aes(x = V4, V1)) + geom_ribbon(aes(ymin = V2, ymax=V3), alpha = 0.2) +  geom_line() + labs(title = "number vs value", x = "V", y = "S")


