import pandas as pd
import numpy as np
from os import listdir
import mygene
import gseapy
import time
import re
import sklearn.feature_selection as skf
import os


flatten = lambda l: [item for sublist in l for item in sublist]
def one_sample_one_tailed(sample_data, popmean, alpha=0.05, alternative='less'):
    t, p = stats.ttest_1samp(sample_data, popmean)

    if alternative == 'greater':
        return(p/2)
    elif alternative == 'less':
        return(1-p/2)

os.chdir("testsuite/")
from scipy import stats

def get_pathways(condition):
    '''
    Maps condition name to a pathway set
    '''
    if condition == "ALS":
        return ['hsa05014']
    elif condition == "LC":
        return ['hsa05223']
    elif condition== "UC":
        return ['hsa04060', 'hsa04630', 'hsa05321']
    elif condition== "HD":
        return ['hsa05016']
    elif condition== "CD":
        return ['hsa04621', 'hsa04060', 'hsa04630', 'hsa05321', 'hsa04140']
# 1. Reading the data and aggregating in a single dataframe

# path to the folder with individual algorithms results
mypath = "../results/"
all_results = listdir(mypath)

for i in range(len(all_results)):
    if all_results[i].endswith(".csv"):
        if i == 0:
            results = pd.read_csv(mypath + all_results[i])
        else:
            temp = pd.read_csv(mypath + all_results[i])
            #changed the naming during the project (1)
            #temp = temp.replace({"REWIRED": "EXPECTED_DEGREE"})
            #changed the naming during the project (2)
            #temp = temp.replace({"RDPN": "REWIRED"})
            results = pd.concat([results,temp], axis=0, ignore_index=True)
results = results.replace({"GSE112680":"ALS","GSE30219": "LC","GSE75214": "UC", "GSE75214_cd": "CD", "GSE3790": "HD"})

#excluded the method
results = results[results.algorithm_name != "PINNACLEZ"]


# 2. Converting gene IDs to symbols
all_genes = list(results.result_genes)
genes = []
for i in range(len(all_genes)):
    try:
        genes.append(all_genes[i].split(","))
    except:
        pass

genes = list(set(flatten(genes)))
mg = mygene.MyGeneInfo()
out = mg.querymany(genes, scopes="entrezgene", fields='symbol', species='human', verbose=False)
mapping = dict()
for line in out:
    try:
        mapping[line["query"]] = line["symbol"]
    except KeyError:
        print("{0} was not mapped to any gene name".format(line["query"]))
        mapping[line["query"]] = line["query"]

# 3. Recomputing p-values for the results where computation failed previously (p-values == -1)
kegg_pval = []

for i in results.index :
    if results.neg_log_gsea_p_value[i] == -1:
        gs = str(results.result_genes[i])
        gs = gs.split(",")
        gs = [mapping[x] for x in gs]

        res = False
        while res == False:
            try:
                enr = gseapy.enrichr(gene_list=gs, description='pathway', gene_sets='KEGG_2016', outdir="out")
                res = True
            except:
                print("sleeping")
                time.sleep(1)
        full_results = enr.results
        terms = list(full_results.Term)
        terms = [x.split(' ')[-1] for x in terms]
        p_values = []
        pathways = get_pathways(results.condition_name[i])
        for j in range(len(terms)):
            if terms[j] in pathways:
                p_values.append(-np.log10(full_results['Adjusted P-value'][j]))

        kegg_pval.append(np.mean(p_values))

        print(i)

#updating the table with recomputed p-values
next = 0
for i in results.index:

    if results.neg_log_gsea_p_value[i] == -1:
        num = kegg_pval[next]
        if np.isnan(num):
            num = 0
        results.neg_log_gsea_p_value[i] = 0

# changing generators names to categorical for visualization
results.network_generator_name=pd.Categorical(results.network_generator_name,categories=['ORIGINAL', 'REWIRED', 'EXPECTED_DEGREE',
                                                                                         "SHUFFLED",
                                                                                         "SCALE_FREE",
                                                                                         "UNIFORM"])
results=results.sort_values('network_generator_name')
# 4. Loading disgenet  network
disgenet = pd.read_csv("../data/networks/disgenet.csv")
dis_ids = {"ALS":"C0002736", "LC": "C1737250", "UC": "C0009324","HD":"C0020179","CD":"C0021390"}
disgenet = disgenet[disgenet["disease_id"].isin(list(dis_ids.values()))]


# computing overlaps for disgenet
overlaps = []
for i in results.index:
    try:
        gs = results.result_genes[i]
        condition = results["condition_name"][i]
        dis = disgenet[disgenet.disease_id == dis_ids[condition]]
        dis_genes = list(dis.gene)
        dis_genes = [str(x) for x in dis_genes]
        gs = gs.split(",")
        oc = len(set(gs).intersection(set(dis_genes)))/min(len(set(gs)), len(set(dis_genes)))
        overlaps.append(oc)
    except:
        overlaps.append(0)

results["disgenet_overlap"] = overlaps

# 5. Loading and preprocessing survival data
# ALS survival
survival_als = pd.read_csv("../data/expression/GSE112680/clinical.csv")
survival_als["characteristics_ch1.1"] = survival_als["characteristics_ch1.1"].apply(lambda x: x.split(" ")[1])
survival_als = survival_als[(survival_als["characteristics_ch1.1"] == "ALS")]["characteristics_ch1.4"]
expr_als = pd.read_csv("../data/expression/GSE112680/expr.csv.zip")
expr_als = expr_als.set_index("Unnamed: 0")
expr_als = expr_als.loc[survival_als.index]
expr_als = expr_als.reindex(survival_als.index)

# LC survival
survival_lc = pd.read_csv("../data/expression/GSE30219/clinical.csv")
d = survival_lc["title"].apply(lambda x: x[:-1])
ct = [re.sub("[\d-]","",x) for x in d]
survival_lc["title"] = ct
survival_lc = survival_lc[(survival_lc["title"] != "Small Cell Lung Carcinoma") & (survival_lc["title"] != "Lung cancer_Other") & (survival_lc["title"] != "Non Tumoral lung")]
survival_lc = survival_lc[["characteristics_ch1.9"]]
survival_lc = survival_lc.dropna()
survival_lc["characteristics_ch1.9"] = [x.split(" ")[-1] for x in survival_lc["characteristics_ch1.9"]]
survival_lc = survival_lc[survival_lc["characteristics_ch1.9"]!= "na"]

expr_ls = pd.read_csv("../data/expression/GSE30219/expr.csv.zip")
expr_ls = expr_ls.set_index("Unnamed: 0")
expr_ls = expr_ls.loc[survival_lc.index]
expr_ls = expr_ls.reindex(survival_lc.index)


s_als = list(survival_als)
s_als = np.array([int(float(x.split(" ")[1])*12) for x in s_als])
s_lc = list(survival_lc["characteristics_ch1.9"])
s_lc = np.array([int(x) for x in s_lc])

# subsetting data that does have survival information
surv_subset = results[(results.condition_name == "ALS") |(results.condition_name == "LC")]


# 6. Computing MMI wrt survival
survival = []
count = 0
for tup in results.itertuples():
    gs = tup[12]
    try:
        gs = gs.split(",")
        condition = tup[4]
        if condition == "ALS":
            mutual_information = skf.mutual_info_classif(expr_als[gs], s_als)
            survival.append(np.mean(mutual_information))
        elif condition == "LC":
            mutual_information = skf.mutual_info_classif(expr_ls[gs], s_lc)
            survival.append(np.mean(mutual_information))
        else:
            survival.append(0)
        if count%10 == 0:
            print(count)
        count += 1
    except:
        survival.append(0)


results["survival"] = survival

# in case you want to join your new data with old:

# results_old = pd.read_csv("../full_results.csv")
# results_full = pd.concat([results_old,results], axis=0, ignore_index=True)
#results_full.to_csv("../full_results.csv")

results.to_csv("../results/all_results.csv")

