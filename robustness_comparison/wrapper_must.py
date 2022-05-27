import os
import random
import pandas as pd
from subprocess import call


def run_must(path_to_network, path_to_seeds, outfile):
    seed_name = path_to_seeds.split("/")[4].split(".")[0]
    file = open(path_to_network)
    lines = file.readlines()
    lines.pop(0)
    random.shuffle(lines)
    tmpfile = f'shuffled_network{seed_name}.txt'
    f = open(tmpfile, "w")
    f.writelines(lines)
    f.close()
    if not os.path.exists('robustness_comparison/MUSTOut'):
        os.makedirs('robustness_comparison/MUSTOut')
    call(["python3", "amim_test_suite/algorithms/must/must.py", tmpfile, path_to_seeds, "0.0", "30", "0", outfile])
    output_graph = pd.read_csv(outfile, sep="\t")
    nodes = set(output_graph["source_node"]).union(set(output_graph["target_node"]))
    os.remove(tmpfile)
    os.remove(outfile)
    return nodes