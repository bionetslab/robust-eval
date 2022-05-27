library(data.table)
library(ggplot2)
robustness_tests <- lapply(list.files("../robustness_comparison/ROBUSTOut_Hyperparameters/", full.names = T), fread)
names(robustness_tests) <- tstrsplit(tstrsplit(list.files("../robustness_comparison/ROBUSTOut_Hyperparameters/"), "ROBUST_", keep=2)[[1]], ".out", keep=1)[[1]]

robustness_tests <- rbindlist(robustness_tests, idcol="init_red")
library(tidyr)
robustness_tests <- robustness_tests %>% 
  separate(col="init_red", into=c("Initial Fraction", "Reduction Factor"), sep="_") %>%
  separate(col="seed_set", into=c("seed_set", "Threshold", "Nr of Trees"), sep="_")
robustness_tests[, Threshold := tstrsplit(Threshold, "key", keep=2)]
robustness_tests[, Threshold := as.factor(Threshold)]
robustness_tests[, `Initial Fraction` := `Initial Fraction`]
robustness_tests[, `Reduction Factor` := `Reduction Factor`]
robustness_tests[, `Nr of Trees` := factor(`Nr of Trees`, levels = c(5, 10, 15, 20, 25, 30))]

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(robustness_tests, aes(x=Threshold, y = `mean jaccard`, color = `Nr of Trees`))+
  geom_boxplot()+
  labs(x = expression(paste("Threshold ", tau)), 
       y=expression(paste("Distributions of robustness coefficients ", r[S], " over 100 seed sets S")))+
  facet_grid(`Reduction Factor` ~ `Initial Fraction`, labeller=label_bquote(
    rows=beta*"="*.(`Reduction Factor`),
    cols=alpha*"="*.(`Initial Fraction`))
    )+
  scale_colour_manual(name = "#Trees", values=colorBlindBlack8[c(6,3,2,8,4,5)], labels=c("n=5", "n=10", "n=15", "n=20", "n=25", "n=30"))+
  theme_bw()+
  theme(text=element_text(size = 13), legend.position = "top")
ggsave("../img/init_vs_red.png", width = 9, height = 12)

ggplot(robustness_tests[`Initial Fraction` == "0.25" & `Reduction Factor` == "0.9"])+
  geom_boxplot(aes(x = Threshold, y = `mean jaccard`, color = `Nr of Trees`))+
  labs(x = expression(paste("Threshold ", tau)),
       y=expression(atop(paste("Distributions of robustness coefficients ", r[S]), paste("over 100 seed sets S"))))+
  facet_grid(`Reduction Factor` ~ `Initial Fraction`, labeller=label_bquote(
    rows=beta*"="*.(`Reduction Factor`),
    cols=alpha*"="*.(`Initial Fraction`))
  )+
  scale_colour_manual(name = "#Trees", values=colorBlindBlack8[c(6,3,2,8,4,5)], labels=c("n=5", "n=10", "n=15", "n=20", "n=25", "n=30"))+
  theme_bw()+
  ylim(0.5,1)+
  theme(text=element_text(size = 20))
ggsave("../img/trees_025_09.png", width = 10, height = 6.5)
