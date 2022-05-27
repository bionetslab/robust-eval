library(data.table)
library(ggplot2)
library(tidyr)

runtime_table <- fread("../runtime_evaluation/runtime_results.csv")
runtime_table <- as.data.table(runtime_table %>% separate("paramters", into = c("Algorithm", "#Seeds", "#Trees"), sep = "_"))
colnames(runtime_table) <- c("Algorithm", "#Seeds", "#Trees", "Runtime")
runtime_table <- runtime_table[, `#Seeds` := as.numeric(`#Seeds`)]
runtime_table <- runtime_table[, `#Trees` := factor(`#Trees`, levels = c("0", "5", "10", "15", "20", "25", "30"))]

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(runtime_table, aes(y = Runtime, x = `#Seeds`, color = `#Trees`))+
  geom_point()+
  facet_wrap(~ Algorithm, scales = "free", ncol=1)+
  scale_colour_manual(values=colorBlindBlack8[c(1,6,3,2,8,4,5)], labels=c("n=0", "n=5", "n=10", "n=15", "n=20", "n=25", "n=30"))+
  geom_line()+
  labs(y = "Runtime [s]")+
  theme_bw()+
  theme(legend.position = "top", text=element_text(size = 14))

ggsave("../img/runtime.png", height = 12, width = 6)
