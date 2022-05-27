library(data.table)
library(ggplot2)
library(ggpubr)
full_results <- fread("full_results.csv")
full_results[, V1 := NULL]

results_robust<- fread("../amim_test_suite/results/all_results.csv")
results_robust[, V1 := NULL]
results_robust <- results_robust[algorithm_name == "ROBUST"]
results_must <- fread("../amim_test_suite/results/all_results.csv")
results_must[, V1 := NULL]
results_must <- results_must[algorithm_name == "MUST"]
results_must[,algorithm_name := "MuST"]

results_all <- rbind(results_robust, full_results)
results_all <- rbind(results_all, results_must)
results_all[algorithm_name == "DIAMOND", algorithm_name := "DIAMOnD"]
results_all <- results_all[algorithm_name %in% c("ROBUST", "MuST", "DOMINO", "DIAMOnD")]
results_all[, network_generator_name := as.factor(network_generator_name)]
results_all[, algorithm_name := as.factor(algorithm_name)]

colorBlind   <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

results_original <- results_all[network_generator_name == "ORIGINAL"]
results_original <- results_original[, -c("mean_mutual_information", "survival")]
results_original <- melt(results_original, measure.vars = c("disgenet_overlap", "neg_log_gsea_p_value"), 
                         variable.name = "score", value.name = "value")
results_original <- results_original[, score := gsub("disgenet_overlap", "Overlap with DisGeNET disease genes", score)]
results_original <- results_original[, score := gsub("neg_log_gsea_p_value", "KEGG gene set enrichment", score)]
results_original <- results_original[, algorithm_name := factor(algorithm_name, levels = c("ROBUST", "MuST", "DIAMOnD", "DOMINO"))]

per_disease_disgenet <- ggplot(results_original[score == "Overlap with DisGeNET disease genes"], aes(x = condition_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  theme(text = element_text(size=14))+
  labs(x = "Disease", y = "Overlap Coefficient")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])

per_disease_kegg <- ggplot(results_original[score == "KEGG gene set enrichment"], aes(x = condition_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  theme(text = element_text(size=14))+
  labs(x = "Disease", y = "-log 10(P-value)")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])

overall_disgenet <- ggplot(results_original[score == "Overlap with DisGeNET disease genes"], aes(x = algorithm_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  labs(x = "Algorithm", y = "Overlap Coefficient")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])+
  theme(text = element_text(size=14), legend.position = "none")

overall_kegg <- ggplot(results_original[score == "KEGG gene set enrichment"], aes(x = algorithm_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  labs(x = "Algorithm", y = "-log 10(P-value)")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])+
  theme(text = element_text(size=14), legend.position = "none")

#plot for the paper
ggarrange(per_disease_disgenet, per_disease_kegg, 
          overall_disgenet, overall_kegg,
          labels=c('A', 'B', 'C', 'D'),
          font.label=list(size=20),
          ncol=1,
          legend = "top",
          common.legend = T)
ggsave("../img/functional_relevance_all.png", height = 12, width = 6)


results_merged <- full_results[!algorithm_name %in% c("HOTNET", "NETCORE")]
results_merged <- rbind(results_robust, results_merged)
results_merged[, algorithm_name := as.factor(algorithm_name)]

compute_pvalues <- function(results, variable){
  return_list <- lapply(unique(results[, network_generator_name]), function(x){
    dt <- data.table(original = character(), generator = character(), algorithm = character(), p_value = numeric())
    for(y in unique(results[, algorithm_name])){
      print(paste(x, ",", y))
      p_val <- wilcox.test(results[network_generator_name == "ORIGINAL" & algorithm_name == y, get(variable)], 
                           results[network_generator_name == x & algorithm_name == y, get(variable)], alternative="greater")$p.value
      dt <- rbind(dt, data.table(original = "ORIGINAL", generator = x, algorithm = y, p_value = p_val))
    } 
    
    return(dt)
  })
  return_dt <- rbindlist(return_list)
  return_dt <- return_dt[generator != "ORIGINAL"]
  return_dt[, p_adj := p.adjust(p_value)]
  return_dt <- return_dt[, algorithm := factor(algorithm, levels = c("ROBUST", "DIAMOnD", "DOMINO", "GF", "COSINE", "KPM", "GXNA", "CLUSTEX2", "GIGA"))]
  return_dt <- return_dt[, generator := factor(generator, levels = c("REWIRED", "EXPECTED_DEGREE", 
                                                                           "SHUFFLED", "SCALE_FREE", "UNIFORM"))]
  return(return_dt)
}

compute_validity_scores <- function(results, variable){
  if(variable == "neg_log_gsea_p_value"){
    validity_scores <- results[network_generator_name == "ORIGINAL" & neg_log_gsea_p_value >= -log10(0.05), .N , by = algorithm_name]
  }else{
    validity_scores <- results[network_generator_name == "ORIGINAL" & get(variable) >= 0.2, .N , by = algorithm_name]
  }
  validity_scores <- merge(validity_scores, results[network_generator_name == "ORIGINAL", .N , by = algorithm_name], by = "algorithm_name", all=T)
  validity_scores$N.x <- nafill(validity_scores$N.x, type="const", fill=0)
  validity_scores$validity_score <- validity_scores$N.x / validity_scores$N.y
  return(validity_scores)
}

results_gsea <- compute_pvalues(results_merged, "neg_log_gsea_p_value")
validity_scores <- compute_validity_scores(results_merged, "neg_log_gsea_p_value")
results_gsea <- merge(results_gsea, validity_scores[, c(1,4)], by.x = "algorithm", by.y = "algorithm_name", all.x=T)

p_values_kegg <- ggplot(results_gsea[algorithm %in% c("DOMINO", "DIAMOnD", "ROBUST")], aes(x = generator, y = -log10(p_value), color = algorithm, size=validity_score))+
  geom_point()+
  scale_color_manual(name = "Algorithm", values = colorBlind[c(5,7,9)])+
  scale_size_continuous(breaks = seq(0.0, 1.0, 0.2),name = "Validity", limits = c(0.0, 1.0))+
  geom_hline(yintercept = -log10(0.05), color="red")+
  theme_bw()+
  theme(text = element_text(size=15))+
  labs(x = "Generator", y = "-log10(P-value)", title="KEGG gene set enrichment")


results_disgenet <- compute_pvalues(results_merged, "disgenet_overlap") 
validity_scores <- compute_validity_scores(results_merged, "disgenet_overlap")
results_disgenet <- merge(results_disgenet, validity_scores[, c(1,4)], by.x = "algorithm", by.y = "algorithm_name", all.x=T)
p_values_disgenet <- ggplot(results_disgenet[algorithm %in% c("DOMINO", "DIAMOnD", "ROBUST")], aes(x = generator, y = -log10(p_value), color = algorithm, size=validity_score))+
  geom_point()+
  scale_color_manual(name = "Algorithm", values = colorBlind[c(5,7,9)])+
  scale_size_continuous(breaks = seq(0.0, 1.0, 0.2),name = "Validity", limits = c(0.0, 0.4))+
  geom_hline(yintercept = -log10(0.05), color="red")+
  theme_bw()+
  theme(text = element_text(size=15))+
  labs(x = "Generator", y = "-log10(P-value)", title="Overlap with DisGeNET disease genes")

results_mi <- compute_pvalues(results_merged, "mean_mutual_information")
validity_scores <- compute_validity_scores(results_merged, "mean_mutual_information")
results_mi <- merge(results_mi, validity_scores[, c(1,4)], by.x = "algorithm", by.y = "algorithm_name", all.x=T)
p_values_mi <- ggplot(results_mi, aes(x = generator, y = -log10(p_value), color = algorithm, size=validity_score))+
  geom_point()+
  scale_color_manual(name = "Algorithm", values = colorBlind)+
  scale_size_continuous(breaks = seq(0,1,0.1))+
  geom_hline(yintercept = -log10(0.05), color="red")+
  theme_bw()+
  theme(text = element_text(size=15))+
  labs(x = "Generator", y = "-log10(P-value)", title="Mean mutual information with the phenotype")

results_survival <- compute_pvalues(results_merged, "survival")
validity_scores <- compute_validity_scores(results_merged, "survival")
results_survival <- merge(results_survival, validity_scores[, c(1,4)], by.x = "algorithm", by.y = "algorithm_name", all.x=T)
p_values_survival <- ggplot(results_survival, aes(x = generator, y = -log10(p_value), color = algorithm, size=validity_score))+
  geom_point()+
  scale_color_manual(name = "Algorithm", values = colorBlind)+
  scale_size_continuous(breaks = seq(0,1,0.1))+
  geom_hline(yintercept = -log10(0.05), color="red")+
  theme_bw()+
  theme(text = element_text(size=15))+
  labs(x = "Generator", y = "-log10(P-value)", title="Mean mutual information with the survival")

ggarrange(p_values_kegg, p_values_disgenet,
          labels=c("A", "B"),
          font.label=list(size=18),
          legend = "right",
          common.legend = F)

ggsave("../img/comparison_AMIM_paper.png", width=18, height = 4)
