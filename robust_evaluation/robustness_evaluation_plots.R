library(data.table)
library(ggplot2)

all_results_must <- fread("../robustness_comparison/MUSTOut/MUST.out", keepLeadingZeros = T)
colnames(all_results_must) <- c("seed_set", "MuST")

all_results_rmust <- fread("../robustness_comparison/RMuSTOut/RMUST.out")
all_results_rmust <- tidyr::separate(all_results_rmust, "seed_set", c("seed_set", "threshold"), sep = "_")
all_results_rmust[, threshold := tstrsplit(threshold, "thr", keep=2)]
colnames(all_results_rmust) <- c("seed_set", "threshold", "R-MuST")

all_results_diamond <- fread("../robustness_comparison/DIAMONDOut/DIAMOND.out", keepLeadingZeros = T)
colnames(all_results_diamond) <- c("seed_set", "DIAMOnD")

all_results_domino <- fread("../robustness_comparison/DOMINOOut/DOMINO.out", keepLeadingZeros = T)
colnames(all_results_domino) <- c("seed_set", "DOMINO")

all_results_robust <- fread("../robustness_comparison/ROBUSTOut/ROBUST_0.25_0.9.out")
all_results_robust <- tidyr::separate(all_results_robust, "seed_set", c("seed_set", "threshold"), sep = "_")
all_results_robust[, threshold := tstrsplit(threshold, "thr", keep=2)]
colnames(all_results_robust) <- c("seed_set", "threshold", "ROBUST")

all_results_wide <- merge(all_results_must, all_results_rmust[, c(1,3)], by = "seed_set")
all_results_wide <- merge(all_results_wide, all_results_diamond, by = "seed_set")
all_results_wide <- merge(all_results_wide, all_results_domino, by = "seed_set")
all_results_wide <- merge(all_results_wide, all_results_robust[, c(1,3)], by = "seed_set")
all_results <- melt(all_results_wide, id.vars = "seed_set", variable.name = "algorithm", value.name = "mean jaccard")
all_results[, algorithm := factor(algorithm, levels = c("DOMINO", "MuST", "R-MuST", "DIAMOnD", "ROBUST"))]

#exclude seed files used for the hyperparameter tuning
robust_hyperparameters <- fread("../robustness_comparison/ROBUSTOut_Hyperparameters/ROBUST_0.25_0.1.out")
robust_hyperparameters <- tidyr::separate(robust_hyperparameters, "seed_set",  c("seed_set", "threshold", "#trees"), sep = "_")
hyperparameter_seeds <- unique(robust_hyperparameters$seed_set)
all_results <- all_results[!seed_set %in% hyperparameter_seeds]
all_results_wide <- all_results_wide[!seed_set %in% hyperparameter_seeds]
print(length(unique(all_results$seed_set)))

colorBlind   <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#2 plots for paper
ggplot(all_results, aes(x = algorithm, y = `mean jaccard`, fill = algorithm))+
  geom_boxplot()+
  scale_fill_manual(name = "Algorithm", values = c("ROBUST" = colorBlind[6], "MuST" = colorBlind[9],
                                                   "DIAMOnD" = colorBlind[4], "DOMINO" = colorBlind[5], 
                                                   "R-MuST" = colorBlind[3]))+
  #scale_color_manual(name = "Algorithm", values = c("ROBUST" = colorBlind[6], "MuST" = colorBlind[9],"DIAMOnD" = colorBlind[4], "DOMINO" = colorBlind[5], "R-MuST" = colorBlind[3]))+
  theme_bw()+
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_text(hjust=2))+
  labs(x = "Algorithm", 
       y = expression(atop(paste("Distributions of robustness coefficients ", r[S]), paste("over 829 seed sets S"))))

ggsave("../img/all_robustness_results.png", width = 8, height = 6.5)

p.values = data.table(var1 = factor(), var2 = factor(), p.value = numeric())
for(i in seq(2, (length(all_results_wide)))){
  for(j in seq(2, (length(all_results_wide)))){
    if(colnames(all_results_wide)[i] != colnames(all_results_wide)[j]){
      print(paste(i, ", ", j, "->", colnames(all_results_wide)[i], ", ", colnames(all_results_wide)[j]))
      tmp <- data.table(
        var1 = colnames(all_results_wide)[i],
        var2 = colnames(all_results_wide)[j],
        p.value = wilcox.test(x = all_results_wide[[i]], y = all_results_wide[[j]], alternative = "greater")$p.value
      )
      p.values <- rbind(p.values, tmp)
    }
  }
}

p.values[, p.adj := p.adjust(p.value, method = "BH")]
p.values <- p.values[order(p.adj)]

library(xtable)
print(xtable(p.values, digits = -3))

