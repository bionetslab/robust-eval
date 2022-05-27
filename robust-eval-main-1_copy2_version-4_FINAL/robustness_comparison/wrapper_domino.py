import os
import random
from subprocess import call, check_output, CalledProcessError


def run_domino(path_to_network, path_to_seeds, outfile):
    seed_name = seed_name = path_to_seeds.split("/")[4].split(".")[0]
    file = open(path_to_network)
    lines = file.readlines()
    lines.pop(0)
    random.shuffle(lines)
    print("Writing shuffled network...")
    shuffled_network = f"shuffled_network_{seed_name}.sif"
    with open(shuffled_network, "w") as sifFile:
        sifFile.write("node_1\tcombined_score\tnode_2\n")
        for line in lines:
            sifFile.write("%s\tpp\t%s" % (line.split("\t")[0], line.split("\t")[1]))

    print("Slicing...")
    domino_slicer = f"domino_slicer_{seed_name}.sif"
    call(["slicer", "--network_file", shuffled_network, "--output_file", domino_slicer])
    print("Calling domino...")
    try:
        output = check_output(["domino", "-a", path_to_seeds, "-n", shuffled_network, "-s", domino_slicer, "-o",
                               "DOMINO/robustness"])
    except CalledProcessError as e:
        print("The domino exception occurred")
        print(e.output)

    print("Getting disease module...")
    nodes = set()
    with open(f"DOMINO/robustness/{seed_name}/modules.out") as file:
        for line in file:
            nodeList = line.split(",")
            nodeList[0] = nodeList[0].split("[")[1]
            nodeList[len(nodeList) - 1] = nodeList[len(nodeList) - 1].split("]")[0]
            for elem in nodeList:
                nodes.add(elem.strip())
    os.remove(shuffled_network)
    os.remove(f"DOMINO/robustness/{seed_name}/modules.out")
    os.remove(domino_slicer)
    os.remove(f'shuffled_network_{seed_name}.domino_slicer_{seed_name}.pkl')
    os.remove(f'shuffled_network_{seed_name}.sif.pkl')
    return nodes