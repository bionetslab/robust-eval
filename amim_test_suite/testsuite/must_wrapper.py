import sys

from testsuite.algorithm_wrapper import AlgorithmWrapper
import networkx as nx
import subprocess
import pandas as pd


#just paste the run_algorithm method into custom_wrapper.py
class MuSTWrapper(AlgorithmWrapper):

    def run_algorithm(self, ggi_network, expression_data, phenotypes, seed_genes, p_values, indicator_matrix, prefix):
        """Runs the algorithm.

        Parameters
        ----------
        ggi_network : nx.Graph
            Possibly randomized GGI network.
        expression_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        phenotypes : np.array, shape (n_samples,)
            Phenotype data (indices are sample IDs).
        seed_genes : list of str
            Seed genes (entries are gene IDs).
        p_values : dict of str: float
            P-values for all genes (keys are gene IDs).
        indicator_matrix : pd.DataFrame
            Indicator matrix obtained from expression data (indices are sample IDs, column names are gene IDs).
        prefix : str
            Prefix to be used for temporary files and directories.

        Returns
        -------
        result_genes : list of str
            Set of genes computed by the algorithm.
        mean_degree : float
            Mean degree of the result genes.
        """

        # 1. Write GGI network in format required by your method
        path_to_network = f'../temp/{prefix}_MuST_ggi.txt'
        with open(path_to_network, 'w') as edge_list_file:
            gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
            edge_list_file.write('source_protein\ttarget_protein\n')
            for u, v in ggi_network.edges():
                edge_list_file.write(f'{gene_ids[u]}\t{gene_ids[v]}\n')

        # 2. Write expression data in format required by your method
        # path_to_expression = f'../temp/{prefix}_MuST_expr.txt'
        # expression_data["phenotype"] = phenotypes
        # expression_data.to_csv(path_to_expression, index=False)

        # 3. Write seed genes in format required by MuST.
        path_seeds = f'../temp/{prefix}_MuST_seed_genes.txt'
        AlgorithmWrapper.save_array(seed_genes, path_seeds, '\n', None)

        # 4. Insert the command to run your method, direct the output to path_to_output
        path_to_output = f'../temp/{prefix}_MuST_output.csv'

        # Run Biosteiner.
        covex_must = 'cd ../algorithms/must/; python must.py'
        command = f'{covex_must} ../{path_to_network} ../{path_seeds} 0.0 30 0 ../{path_to_output}'
        subprocess.call(command, shell=True, stdout=sys.stdout)

        # 5. Process results such that they are formatted as a list of strings (entez IDs)
        output_graph = pd.read_csv(path_to_output, sep="\t")
        result_genes = set(output_graph["source_node"]).union(set(output_graph["target_node"]))
        result_genes = list(map(str, result_genes))

        # 6. Delete temporary data.
        subprocess.call(f'rm ../temp/{prefix}_MuST_*', shell=True)

        return result_genes, AlgorithmWrapper.mean_degree(ggi_network, result_genes)