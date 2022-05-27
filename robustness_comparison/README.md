<h2>How to run the robustness tests</h2>


0. Activate your environment / install all the packages needed and 
   make sure that robust-eval is your working directory. Make sure you also pulled the submodule of robust with `git submodule update --init`.
1. Switch to the robust-eval directory
2. Decide which parameters you need: 

- `--algorithm: ` must be one of ROBUST, DOMINO, DIAMOND, MUST, RMUST
- `--threads:` Number of threads that can should be used so that the evaluation can be run at the same time 
for different seed files.
- `--threshold:` only applicable for ROBUST or RMUST: Both algorithms compute 10 Steiner Trees. 
By specifying a threshold between 0.1 and 0.9, you filter the results such that only nodes that 
are contained in >=threshold*100% of the trees are returned. It's possible to specify multiple 
  thresholds at once separated by a space. 
- `--init:` only applicable for ROBUST: Initial fraction value.
- `--red:` only applicable for ROBUST: Reduction factor value

3. For e.g. ROBUST, run:

`python robustness_comparison/robustness_mean_jaccard.py --algorithm ROBUST --threads 4 --threshold 0.6 0.7 --init 0.25 --red 0.9 `
