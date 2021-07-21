#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --job-name=robustness1
#SBATCH --output=%x.%j.txt
#SBATCH --error=%x.%j.err
#SBATCH --mem=100G

python robustness_comparison/robustness_mean_jaccard.py --algorithm ROBUST --threshold 0.5 --init 0.25 --red 0.9 --threads 25 --trees 20
python robustness_comparison/robustness_mean_jaccard.py --algorithm DIAMOND --threads 25
python robustness_comparison/robustness_mean_jaccard.py --algorithm DOMINO --threads 25
