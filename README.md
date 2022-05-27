# robust-eval
Evaluation of the ROBUST disease module mining method

This repository contains different directories: 

1. The AMIM test suite including the custom wrappers for ROBUST and MuST in the `amim_test_suite` directory
2. The robustness tests for DOMINO, DIAMOnD, MuST, R-MuST and ROBUST in the `robustness_comparison` directory
3. The script for the runtime evaluation are located in the `runtime_evaluation` directory
4. The R scripts for the plots shown in the paper in the `robust_evaluation` directory

# AMIM test suite
The AMIM test suite can be run for ROBUST by calling
```bash
python amim_test_suite/run_tests.py parallel --networks BIOGRID HPRD STRING APID IID_BRAIN IID_LUNG --generators ORIGINAL REWIRED EXPECTED_DEGREE SHUFFLED SCALE_FREE UNIFORM --method ROBUST --verbose
```
Other parameters:
```bash
--methods {KPM,DIAMOND,GXNA,CLUSTEX2,PINNACLEZ,GIGA,GF,COSINE,HOTNET,NETCORE,CUSTOM,DOMINO,ROBUST,MUST}
--networks {BIOGRID,HPRD,STRING,APID,IID,IID_BRAIN,IID_LUNG,CUSTOM}
--generators {ORIGINAL,EXPECTED_DEGREE,SHUFFLED,SCALE_FREE,UNIFORM,REWIRED}
--conditions {GSE112680,GSE30219,GSE75214,GSE75214_cd,GSE3790}
--verbose
```
# Robustness tests
The robustness tests can be run for ROBUST by calling
```bash
python robustness_comparison/robustness_mean_jaccard.py --algorithm ROBUST --threshold 0.1 --init 0.25 --red 0.9 --threads 32 --trees 30
```
Other parameters:
```bash
--algorithm {ROBUST,DOMINO,DIAMOND,MUST,RMUST}
--threshold THRESHOLD [THRESHOLD ...]   Threshold for ROBUST and R-MuST, default: 0.1-0.9 in 0.1 steps
--init INIT                             Initial fraction value for ROBUST, default: 0.25
--red RED                               Reduction factor value for ROBUST, default: 0.9
--threads THREADS                       Run on multiple seed files at once with x threads   
 --trees TREES [TREES ...]              Number of trees to be returned.                      
```

# Runtime evaluation 
The runtime was evaluated by calling
```bash
python runtime_evaluation/compare_runtimes.py
```