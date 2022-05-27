from testsuite.test_runner import TestRunner
import testsuite.utils as utils
import argparse
import os
import multiprocessing as mp
import itertools as itt
import traceback


def get_parser():
    parser = argparse.ArgumentParser('tests the one-network-fits-all hypothesis')
    subparsers = parser.add_subparsers(dest='mode', required=True, help='run tests sequentially or in parallel')
    sequential = subparsers.add_parser('sequential', help='runs the tests sequentially')
    sequential.add_argument('--network', type=utils.GGINetworkSelector, choices=list(utils.GGINetworkSelector), required=True)
    sequential.add_argument('--edge_cost', type=str, choices=['UNIFORM', 'ADDITIVE', 'EXPONENTIAL'], default='UNIFORM')
    sequential.add_argument('--node_namespace', type=str, choices=['ENTREZ_GENE_ID', 'GENE_SYMBOL', 'UNIPROT_PROTEIN_ID'], default='GENE_SYMBOL')
    sequential.add_argument('--normalize', type=str, default='BAIT_USAGE', choices=['BAIT_USAGE', 'STUDY_ATTENTION', 'CUSTOM'], help='Specifies edge weight function used by ROBUST')
    sequential.add_argument('--lambda_', type=float, default=0.5, help='Hyper-parameter lambda used by PAIR_FREQ, BAIT_USAGE and CUSTOM edge weights. Should be set to value between 0 and 1.')
    sequential.add_argument('--generator', type=utils.NetworkGeneratorSelector, choices=list(utils.NetworkGeneratorSelector), required=True)
    sequential.add_argument('--method', type=utils.AlgorithmSelector, choices=list(utils.AlgorithmSelector), required=True)
    sequential.add_argument('--condition', type=utils.ConditionSelector, choices=list(utils.ConditionSelector))
    sequential.add_argument('--verbose', action='store_true', help='print progress to stdout')
    
    parallel = subparsers.add_parser('parallel', help='runs the tests in parallel')
    parallel.add_argument('--methods', type=utils.AlgorithmSelector, choices=list(utils.AlgorithmSelector), nargs='+', default=list(utils.AlgorithmSelector))
    parallel.add_argument('--networks', type=utils.GGINetworkSelector, choices=list(utils.GGINetworkSelector), nargs='+', default=list(utils.GGINetworkSelector))
    parallel.add_argument('--edge_cost', type=str, choices=['UNIFORM', 'ADDITIVE', 'EXPONENTIAL'], nargs='+', help='Hyper-parameter lambda used by PAIR_FREQ, BAIT_USAGE and CUSTOM edge weights. Should be set to value between 0 and 1.', default='UNIFORM')    
    parallel.add_argument('--node_namespace', type=str, choices=['ENTREZ_GENE_ID', 'GENE_SYMBOL', 'UNIPROT_PROTEIN_ID'], nargs='+', help='Hyper-parameter lambda used by PAIR_FREQ, BAIT_USAGE and CUSTOM edge weights. Should be set to value between 0 and 1.', default='GENE_SYMBOL')    
    parallel.add_argument('--normalize', type=str, default='BAIT_USAGE', choices=['BAIT_USAGE', 'STUDY_ATTENTION', 'CUSTOM'], nargs='+', help='Specifies edge weight function used by ROBUST')
    parallel.add_argument('--lambda', type=float, nargs='+', default=0.5, help='Hyper-parameter lambda used by PAIR_FREQ, STUDY_ATTENTION and CUSTOM edge weights. Should be set to value between 0 and 1.')
    parallel.add_argument('--generators', type=utils.NetworkGeneratorSelector, choices=list(utils.NetworkGeneratorSelector), nargs='+', default=list(utils.NetworkGeneratorSelector))
    parallel.add_argument('--conditions', type=utils.ConditionSelector, choices=list(utils.ConditionSelector), nargs='+',default=list(utils.ConditionSelector))
    parallel.add_argument('--verbose', action='store_true', help='print progress to stdout')

    return parser


def run_tests(ggi_network_selector, edge_cost_selector, node_namespace_selector, normalize_selector, lambda_selector, network_generator_selector, algorithm_selector, condition_selectors, verbose=False):
    try:
        if verbose:
            print('loading data ...')
        test_runner = TestRunner()
        if verbose:
            print('running the tests ...')
        test_runner.run_all(ggi_network_selector, edge_cost_selector, node_namespace_selector, normalize_selector, lambda_selector, network_generator_selector, algorithm_selector, condition_selectors, verbose)
        if verbose:
            print('saving the results')
        test_runner.save_results()
        return 0
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    os.chdir('amim_test_suite/testsuite')
    args = get_parser().parse_args()
    if args.mode == 'sequential':
        print('running tests sequentially ...')
        run_tests(args.network, args.edge_cost, args.node_namespace, args.normalize, args.lambda_, args.generator, args.method, [args.condition], args.verbose)

    elif args.mode == 'parallel':
        print('running tests in parallel ...')
        pool = mp.Pool(len(args.networks) * len(args.generators))
        exit_codes = pool.starmap(run_tests, list(itt.product(args.networks, args.edge_cost, args.node_namespace, args.normalize, args.lambda_, args.generators, args.methods,
                                                              [args.conditions], [args.verbose])))
        print(f'exit codes: {exit_codes}')



