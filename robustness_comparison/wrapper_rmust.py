import os
import random
import pandas as pd
from subprocess import call


def run_rmust(path_to_network, path_to_seeds, threshold, nr_of_trees):
    seed_name = path_to_seeds.split("/")[4].split(".")[0]
    file = open(path_to_network)
    lines = file.readlines()
    lines.pop(0)
    random.shuffle(lines)
    if type(threshold) == list:
        thr_name = "_".join(map(str, threshold))
    else:
        thr_name = threshold
    shuffled_network = f"shuffled_network_{seed_name}_{thr_name}.txt"
    f = open(shuffled_network, "w")
    f.writelines(lines)
    f.close()
    outnodes = f"outnodes_{seed_name}_{thr_name}.txt"
    open(outnodes, 'a').close()
    outedges = f"outedges_{seed_name}_{thr_name}.txt"
    open(outedges, 'a').close()
    call(
        f'java -jar robustness_comparison/MultiSteinerBackend.jar -nw {shuffled_network} -s {path_to_seeds} -on {outnodes} -oe {outedges} -hp 0.0 -m -t {nr_of_trees} -mi 2 -pd -ncd 32',
        shell=True)

    nodes = pd.read_csv(outnodes, sep="\t")
    os.remove(shuffled_network)
    os.remove(outnodes)
    os.remove(outedges)
    if type(threshold) == list:
        return_dict = {}
        for thr in threshold:
            return_dict[thr] = set(nodes[nodes["participation_number"] >= thr * 10]["node"])
        return return_dict
    nodes_out = set(nodes[nodes["participation_number"] >= threshold * 10]["node"])
    return nodes_out