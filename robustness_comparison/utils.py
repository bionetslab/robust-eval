import os
from enum import Enum
import argparse


class AlgorithmSelector(Enum):
    ROBUST = 'ROBUST'
    DOMINO = 'DOMINO'
    DIAMOND = 'DIAMOND'
    MUST = 'MUST'
    RMUST = 'RMUST'

    def __str__(self):
        return self.value


def allowed_values_floats(x):
    x = float(x)
    if not 0.0 <= x <= 1.0:
        raise argparse.ArgumentTypeError("Choose a value between 0 and 1")
    return x


def get_parser():
    parser = argparse.ArgumentParser('tests the mean overlap between different runs with shuffled input networks')
    parser.add_argument('--algorithm', type=AlgorithmSelector, choices=list(AlgorithmSelector), required=True)
    parser.add_argument('--threshold', type=allowed_values_floats, nargs='+',
                        help='Threshold for ROBUST and R-MuST, default: 0.1-0.9 in 0.1 steps')
    parser.add_argument('--init', default=0.25, type=allowed_values_floats,
                        help='Initial fraction value for ROBUST, default: 0.25')
    parser.add_argument('--red', default=0.9, type=allowed_values_floats,
                        help='Reduction factor value for ROBUST, default: 0.9')
    parser.add_argument('--threads', default=4, type=int, help='Run on multiple seed files at once with x threads')
    parser.add_argument('--trees', default=10, type=int, help='Number of trees to be returned.', nargs='+')
    return parser


def prep_parameters(algorithm=AlgorithmSelector.ROBUST, threshold=0.5, init=0.25, red=0.9, nr_of_trees=20):
    list_tuples = []

    # for robustness tests
    # sort seeds by size
    path = "robustness_comparison/data/all-seeds"
    dirpath = os.path.abspath(path)
    all_files = (os.path.join(basedir, filename) for basedir, dirs, files in os.walk(dirpath) for filename in files)
    all_seeds = sorted(all_files, key = os.path.getsize)
    all_seeds = [os.path.basename(file) for file in all_seeds][0:900]

    #for hyperparameter tests
    #all_seeds = ['0000088.txt', '0000248.txt', '0000328.txt', '0000507.txt', '0001341.txt', '0001744.txt', '0001751.txt', '0002305.txt', '0002462.txt', '0002476.txt', '0002728.txt', '0002926.txt', '0002927.txt', '0003225.txt', '0003435.txt', '0003436.txt', '0003608.txt', '0003783.txt', '0003999.txt', '0004614.txt', '0004946.txt', '0004985.txt', '0005044.txt', '0005090.txt', '0005279.txt', '0005292.txt', '0005480.txt', '0005721.txt', '0006176.txt', '0006375.txt', '0006406.txt', '0006566.txt', '0006573.txt', '0006639.txt', '0007042.txt', '0007043.txt', '0007269.txt', '0007309.txt', '0007650.txt', '0007863.txt', '0008199.txt', '0008364.txt', '0008748.txt', '0009046.txt', '0009133.txt', '0009279.txt', '0009292.txt', '0009443.txt', '0009480.txt', '0009640.txt', '0009641.txt', '0009655.txt', '0009656.txt', '0009668.txt', '0009696.txt', '0010264.txt', '0010311.txt', '0010847.txt', '0011153.txt', '0011190.txt', '0011191.txt', '0011556.txt', '0011740.txt', '0011812.txt', '0011966.txt', '0013433.txt', '0013792.txt', '0015137.txt', '0015240.txt', '0015254.txt', '0015452.txt', '0015917.txt', '0016364.txt', '0016576.txt', '0016761.txt', '0016985.txt', '0017052.txt', '0017053.txt', '0017319.txt', '0018177.txt', '0018188.txt', '0018214.txt', '0018229.txt', '0018598.txt', '0018604.txt', '0018638.txt', '0018770.txt', '0018943.txt', '0018956.txt', '0019269.txt', '0019308.txt', '0019719.txt', '0020562.txt', '0020563.txt', '0020761.txt', '0021084.txt', '0021085.txt', '0021245.txt', '0023122.txt', '0024498.txt']
    #os.listdir("robustness_comparison/data/all-seeds")[0:100]

    if type(threshold) == list and len(threshold) == 1:
        threshold = threshold[0]
    if type(nr_of_trees) == list and len(nr_of_trees) == 1:
        nr_of_trees = nr_of_trees[0]

    if type(threshold) == list and algorithm in (AlgorithmSelector.ROBUST, AlgorithmSelector.RMUST):
        print(f'Testing {algorithm} with multiple thresholds...')
        if type(nr_of_trees) == list:
            print(f'Testing {algorithm} with multiple numbers of trees...')
            for seed_file in all_seeds:
                path_to_seeds = f"robustness_comparison/data/all-seeds/{seed_file}"
                if algorithm == AlgorithmSelector.ROBUST:
                    robust_out = f"robustness_comparison/{algorithm}Out/ROBUST_{seed_file.split('.')[0]}_init{init}_red{red}.out"
                    list_tuples.append((algorithm.value, path_to_seeds, robust_out, threshold, init, red, nr_of_trees))
                else:
                    raise ValueError("#trees can only be specified as a list for ROBUST")
        else:
            for seed_file in all_seeds:
                path_to_seeds = f"robustness_comparison/data/all-seeds/{seed_file}"
                if algorithm == AlgorithmSelector.ROBUST:
                    robust_out = f"robustness_comparison/{algorithm}Out/ROBUST_{seed_file.split('.')[0]}_init{init}_red{red}.out"
                    list_tuples.append((algorithm.value, path_to_seeds, robust_out, threshold, init, red))
                else:
                    rmust_out = f"robustness_comparison/{algorithm}Out/RMUST_{seed_file.split('.')[0]}.out"
                    list_tuples.append((algorithm.value, path_to_seeds, rmust_out, threshold, nr_of_trees))
    elif algorithm == AlgorithmSelector.ROBUST:
        print(f'Running ROBUST with initial fraction={init} and reduction factor={red}')
        for seed_file in all_seeds:
            path_to_seeds = f"robustness_comparison/data/all-seeds/{seed_file}"
            robust_out = f"robustness_comparison/{algorithm}Out/ROBUST_{seed_file.split('.')[0]}_thr{threshold}_init{init}_red{red}.out"
            list_tuples.append((algorithm.value, path_to_seeds, robust_out, threshold, init, red, nr_of_trees))
    elif algorithm in [AlgorithmSelector.DIAMOND, AlgorithmSelector.DOMINO, AlgorithmSelector.MUST]:
        print(f'Running {algorithm}...')
        for seed_file in all_seeds:
            path_to_seeds = f"robustness_comparison/data/all-seeds/{seed_file}"
            outfile = f"robustness_comparison/{algorithm}Out/{algorithm}_{seed_file.split('.')[0]}.out"
            list_tuples.append((algorithm.value, path_to_seeds, outfile))
    elif algorithm == AlgorithmSelector.RMUST:
        print(f'Running R-MuST with threshold={threshold}')
        for seed_file in all_seeds:
            path_to_seeds = f"robustness_comparison/data/all-seeds/{seed_file}"
            rmust_out = f"robustness_comparison/{algorithm}Out/RMUST_{seed_file.split('.')[0]}_thr{threshold}.out"
            list_tuples.append((algorithm.value, path_to_seeds, rmust_out, threshold, nr_of_trees))
    return list_tuples
