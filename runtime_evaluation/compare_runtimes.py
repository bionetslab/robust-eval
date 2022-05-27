import random
import sys
import os
import subprocess
import pandas as pd
import tempfile
import time


def time_robust(path_to_network, path_to_seeds, path_to_output, nr_of_trees):
    # Run robust.
    robust = 'cd amim_test_suite/algorithms/robust/; python robust.py'
    # input seeds output initial fraction reduction factor #trees threshold
    command = f'{robust} ../../../{path_to_network} {path_to_seeds} {path_to_output} 0.25 0.9 {nr_of_trees} 0.1'
    start_time = time.perf_counter()
    subprocess.call(command, shell=True, stdout=subprocess.PIPE)
    stop_time = time.perf_counter()
    return stop_time-start_time


def time_must(path_to_network, path_to_seeds, path_to_output, nr_of_trees):
    start_time = time.perf_counter()
    subprocess.call(["python3", "amim_test_suite/algorithms/must/must.py", path_to_network, path_to_seeds, "0.0", str(nr_of_trees), "0", path_to_output])
    stop_time = time.perf_counter()
    return stop_time - start_time


def time_rmust(path_to_network, path_to_seeds, nr_of_trees):
    outnodes = tempfile.NamedTemporaryFile()
    outedges = tempfile.NamedTemporaryFile()
    start_time = time.perf_counter()
    subprocess.call(
        f'java -jar robustness_comparison/MultiSteinerBackend.jar -nw {path_to_network} -s {path_to_seeds} -on {outnodes.name} -oe {outedges.name} -hp 0.0 -m -t {nr_of_trees} -mi 2 -pd -ncd 4',
        shell=True)
    stop_time = time.perf_counter()
    return stop_time - start_time


def time_domino(path_to_network, path_to_seeds):
    domino_network = tempfile.NamedTemporaryFile()
    file = open(path_to_network)
    lines = file.readlines()
    lines.pop(0)
    with open(domino_network.name, 'w') as sifFile:
        sifFile.write("node_1\tcombined_score\tnode_2\n")
        for line in lines:
            sifFile.write("%s\tpp\t%s" % (line.split("\t")[0], line.split("\t")[1]))
    print("Slicing...")
    start_time = time.perf_counter()
    domino_slicer = tempfile.NamedTemporaryFile()
    subprocess.call(["slicer", "--network_file", domino_network.name, "--output_file", domino_slicer.name])
    print("Calling domino...")
    tempdir = tempfile.TemporaryDirectory()
    try:
        output = subprocess.check_output(["domino", "-a", path_to_seeds, "-n", domino_network.name, "-s", domino_slicer.name, "-o",
                               tempdir.name])
    except subprocess.CalledProcessError as e:
        print("The domino exception occurred")
        print(e.output)
    stop_time = time.perf_counter()
    return stop_time-start_time


def time_diamond(path_to_network, path_to_seeds, path_to_output):
    start_time = time.perf_counter()
    subprocess.call(
        ["python3", "amim_test_suite/algorithms/DIAMOnD/DIAMOnD.py", path_to_network, path_to_seeds, "200", "1",
         path_to_output])
    stop_time = time.perf_counter()
    return stop_time - start_time


def runtime_comparison():
    #os.chdir("../")
    path_to_network = "robustness_comparison/data/protein-protein-interaction.txt"
    dt = pd.read_csv(path_to_network, sep="\t")
    nodelist = list(set(dt['source_protein']).union(dt['target_protein']))
    #sample seeds
    runtime_dict = {}
    random.seed(1234)
    for i in range(25, 401, 25):
        print(i)
        for j in range(10):
            seeds = random.choices(nodelist, k=i)
            path_to_seeds = tempfile.NamedTemporaryFile()
            with open(path_to_seeds.name, 'w') as f:
                for seed in seeds:
                    f.write(seed + "\n")
            path_to_output = tempfile.NamedTemporaryFile()
            for n_trees in [5, 10, 15, 20, 25, 30]:
                #ROBUST
                runtime_robust = time_robust(path_to_network, path_to_seeds.name, path_to_output.name, n_trees)
                runtime_dict[f'ROBUST_{i}_{j}_{n_trees}'] = runtime_robust
                # R-MuST
                #runtime_rmust = time_rmust(path_to_network, path_to_seeds.name, n_trees)
                #runtime_dict[f'R-MuST_{i}_{n_trees}'] = runtime_rmust
                #MuST
                #runtime_must = time_must(path_to_network, path_to_seeds.name, path_to_output.name, n_trees)
                #runtime_dict[f'MuST_{i}_{n_trees}'] = runtime_must
            # DOMINO
            runtime_domino = time_domino(path_to_network, path_to_seeds.name)
            runtime_dict[f'DOMINO_{i}_{j}_0'] = runtime_domino
            #DIAMOnD
            runtime_diamond = time_diamond(path_to_network, path_to_seeds.name, path_to_output.name)
            runtime_dict[f'DIAMOnD_{i}_{j}_0'] = runtime_diamond
    return runtime_dict


if __name__ == '__main__':
    sys.path.append(os.getcwd())
    runtimes = runtime_comparison()
    filename = 'runtime_evaluation/runtime_results_fast.csv'
    print('Writing file...')
    with open(filename, "w") as file:
        file.write("paramters,runtime\n")
        for key, value in runtimes.items():
            file.write(f'{key},{value}\n')
        file.close()
