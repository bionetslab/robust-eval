import random
import string
import networkx as nx


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def shuffle_nodes(G):
    random_names = [id_generator(size=8) for _ in range(len(G.nodes()))]
    # create a random mapping old label -> new label
    node_mapping = dict(zip(G.nodes(), random_names))
    reverse_map = dict(reversed(item) for item in node_mapping.items())
    # build a new graph
    G_new = nx.relabel_nodes(G, node_mapping)
    return G_new, node_mapping, reverse_map


def run_robust(terminals, threshold, init, red, nr_of_trees):
    from amim_test_suite.algorithms.robust.pcst_approach.utils.ppi import PpiInstance, UnitEdgeWeight, read_ppi_shuffled
    from amim_test_suite.algorithms.robust.pcst_approach.utils import ExpMinMaxDiverseSteinerTreeComputer, SolutionSet
    graph = read_ppi_shuffled("robustness_comparison/data/protein-protein-interaction.txt", shuffle=True)
    graph, node_mapping, reverse_map = shuffle_nodes(graph)
    terminals = [node_mapping[key] for key in terminals]
    edge_weights = UnitEdgeWeight()
    ppi_instance = PpiInstance(graph, terminals, edge_weights)
    engine = ExpMinMaxDiverseSteinerTreeComputer(initial_fraction=init,
                                                 reduction_factor=red)
    if type(nr_of_trees) == list:
        return_dict = {}
        max_n = max(nr_of_trees)
        steiner_trees = engine(ppi_instance, n=max_n - 1)
        for nt in nr_of_trees:
            t = steiner_trees.get_occurrences(include_terminals=True, first_n=nt)
            t.loc[t["#occurrences"] > nt, "#occurrences"] = nt
            t["%occurrences"] = t["#occurrences"] / max(t["#occurrences"])
            if type(threshold) == list:
                for thr in threshold:
                    return_dict[f'{thr}_{nt}'] = t[t["%occurrences"] >= thr].index.values
            else:
                return_dict[nt] = t[t["%occurrences"] >= threshold].index.values
        return return_dict
    else:
        steiner_trees = engine(ppi_instance, n=nr_of_trees)
        t = steiner_trees.get_occurrences(include_terminals=True)
        t = t.rename(index=reverse_map)
        if type(threshold) == list:
            return_dict = {}
            for thr in threshold:
                return_dict[thr] = t[t["%occurrences"] >= thr].index.values
            return return_dict
        return set(t[t["%occurrences"] >= threshold].index.values)


def run_robust_original(terminals, threshold, init, red, nr_of_trees):
    from amim_test_suite.algorithms.robust.pcst_approach.utils.ppi import PpiInstance, UnitEdgeWeight, read_ppi_shuffled
    from amim_test_suite.algorithms.robust.pcst_approach.utils import ExpMinMaxDiverseSteinerTreeComputer, SolutionSet
    graph = read_ppi_shuffled("robustness_comparison/data/protein-protein-interaction.txt", shuffle=True)
    edge_weights = UnitEdgeWeight()
    ppi_instance = PpiInstance(graph, terminals, edge_weights)
    engine = ExpMinMaxDiverseSteinerTreeComputer(initial_fraction=init,
                                                 reduction_factor=red)
    if type(nr_of_trees) == list:
        return_dict = {}
        max_n = max(nr_of_trees)
        steiner_trees = engine(ppi_instance, n=max_n-1)
        for nt in nr_of_trees:
            t = steiner_trees.get_occurrences(include_terminals=True, first_n=nt)
            t.loc[t["#occurrences"] > nt, "#occurrences"] = nt
            t["%occurrences"] = t["#occurrences"] / max(t["#occurrences"])
            if type(threshold) == list:
                for thr in threshold:
                    return_dict[f'{thr}_{nt}'] = t[t["%occurrences"] >= thr].index.values
            else:
                return_dict[nt] = t[t["%occurrences"] >= threshold].index.values
        return return_dict
    else:
        steiner_trees = engine(ppi_instance, n=nr_of_trees)
        t = steiner_trees.get_occurrences(include_terminals=True)
        if type(threshold) == list:
            return_dict = {}
            for thr in threshold:
                return_dict[thr] = t[t["%occurrences"] >= thr].index.values
            return return_dict
        return set(t[t["%occurrences"] >= threshold].index.values)