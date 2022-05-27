import sys
from pcst_approach.utils.ppi import PpiInstance, read_terminals, UnitEdgeWeight, read_ppi, read_ppi_biasaware, BiasAwareEdgeWeight_Additive, BiasAwareEdgeWeight_Exponential
from pcst_approach.utils import ExpMinMaxDiverseSteinerTreeComputer
import networkx as nx
import numpy as np
import pandas as pd



def call_robust(path_to_graph, edge_cost_mode, node_namespace_mode, normalize_mode, lambda_value, path_to_seeds, outfile, init, red, numberOfSteinerTrees, threshold):
    import networkx as nx
    
    lambda_=lambda_value
    
    if lambda_>1 or lambda_<0:
        lambda_=0.00
    
        
    if normalize_mode=='BAIT_USAGE':
        if node_namespace_mode=='ENTREZ_GENE_ID':
            pathToStudyBiasData='./data/edgeweightdata/entrezID_baitUsage.txt'
        elif node_namespace_mode=='GENE_SYMBOL':
            pathToStudyBiasData='./data/edgeweightdata/geneSymbol_baitUsage.txt'
        elif node_namespace_mode=='UNIPROT_PROTEIN_ID':
            pathToStudyBiasData='./data/edgeweightdata/uniprotID_baitUsage.txt'
        flag=int(1)
            
    elif normalize_mode=='STUDY_ATTENTION':
        if node_namespace_mode=='ENTREZ_GENE_ID':
            pathToStudyBiasData='./data/edgeweightdata/entrezID_studyAttention.txt'
        elif node_namespace_mode=='GENE_SYMBOL':
            pathToStudyBiasData='./data/edgeweightdata/geneSymbol_studyAttentionn.txt'
        elif node_namespace_mode=='UNIPROT_PROTEIN_ID':
            pathToStudyBiasData='./data/edgeweightdata/uniprotID_studyAttention.txt'
        flag=int(2)
        
    elif normalize_mode=='CUSTOM':
        if node_namespace_mode=='ENTREZ_GENE_ID':
            pathToStudyBiasData='./data/edgeweightdata/custom.txt'
        elif node_namespace_mode=='GENE_SYMBOL':
            pathToStudyBiasData='./data/edgeweightdata/custom.txt'
        elif node_namespace_mode=='UNIPROT_PROTEIN_ID':
            pathToStudyBiasData='./data/edgeweightdata/custom.txt'
        flag=int(3)
        
    if edge_cost_mode=='UNIFORM':
        graph = read_ppi(path_to_graph)
    else:
        graph = read_ppi_biasaware(path_to_graph, pathToStudyBiasData, flag)
        
    terminals = read_terminals(path_to_seeds)
    #kick out terminals not in graph
    terminals = list(set(terminals).intersection(set(graph.nodes)))
        
    if edge_cost_mode=='UNIFORM':
        edge_weights = UnitEdgeWeight()
    else:
        if edge_cost_mode=='ADDITIVE':
            edge_weights = BiasAwareEdgeWeight_Additive(graph, lambda_)
        elif edge_cost_mode=='EXPONENTIAL':
            edge_weights = BiasAwareEdgeWeight_Exponential(graph, lambda_)
            
        
    ppi_instance = PpiInstance(graph, terminals, edge_weights)
    
    nx.write_graphml_lxml(graph, "TEMP.graphml")
    
    # 2. Solving the instance
    engine = ExpMinMaxDiverseSteinerTreeComputer(initial_fraction=init,
                                                 reduction_factor=red)
    # The most important parameter seems to be initial_fraction.
    steiner_trees = engine(ppi_instance, n=numberOfSteinerTrees)
    t = steiner_trees.get_occurrences(include_terminals=True)
    t = t[t["%occurrences"] >= threshold]
    subgraph = steiner_trees.get_subgraph(threshold=threshold)
    comp_idx = 0
    for comp in sorted(nx.connected_components(subgraph), key=len, reverse=True):
        for node in comp:
            subgraph.nodes[node]['connected_components_id'] = comp_idx
        comp_idx += 1
    print("Writing results...")
    if outfile.endswith(".csv"):
        t.to_csv(path_to_outfile)
    elif outfile.endswith(".graphml"):
        nx.write_graphml(subgraph, path_to_outfile)
    else:
        nx.write_edgelist(subgraph, path_to_outfile, data=False)


if __name__ == '__main__':
    # -----------------------------------------------------
    # Checking for input from the command line:
    # -----------------------------------------------------
    #
    # [1] file providing the network in the form of an edgelist
    #     (tab-separated table, columns 1 & 2 will be used)
    #
    # [2] file with the seed genes (if table contains more than one
    #     column they must be tab-separated; the first column will be
    #     used only)
    #
    # [3] path to output file
    #
    # [4] initial fraction
    # [5] reduction factor
    # [6] number of steiner trees to be computed
    # [7] threshold
    input_list = sys.argv
    print("Parsing input...")
    path_to_graph = str(input_list[1])
    edge_cost_mode=str(input_list[2])
    node_namespace_mode=str(input_list[3])
    normalize_mode=str(input_list[4])
    lambda_value=float(input_list[5])
    path_to_seeds = str(input_list[6])
    path_to_outfile = str(input_list[7])
    initial_fraction = float(input_list[8])
    reduction_factor = float(input_list[9])
    number_of_steiner_trees = int(input_list[10])
    threshold = float(input_list[11])

    print(f"Computing Steiner Trees with the following parameters: \n "
          f"graph: {path_to_graph}\n"
          f"edgeWeightMode: {path_to_graph}\n"
          f"seeds: {path_to_seeds}\n"
          f"outfile: {path_to_outfile}\n"
          f"initial fraction: {initial_fraction}\n"
          f"reduction factor: {reduction_factor}\n"
          f"number of steiner trees: {number_of_steiner_trees}\n"
          f"threshold: {threshold}")
    number_of_steiner_trees -=1
    call_robust(path_to_graph, edge_cost_mode, node_namespace_mode, normalize_mode, lambda_value, path_to_seeds, path_to_outfile, initial_fraction, reduction_factor, number_of_steiner_trees, threshold)
    