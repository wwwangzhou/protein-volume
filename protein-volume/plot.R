install.packages("tidyverse")
library(ggplot2)
final = read.table("~/Desktop/129.txt")
ggplot(final, aes(x = V4, V1)) + geom_ribbon(aes(ymin = V2, ymax=V3), alpha = 0.2) +  geom_line() + labs(title = "number vs volume", x = "N", y = "V")


