This repository contains scripts regularly used for ancestral sequence resurrection.

PAML-DATA:
the script 'asr-pipeline_paml.py' can be used to process data calculated with PAML's codeml(1).

IQTREE-DATA:
the script asr-pipeline_iqtree.py can be used to process data calculated with the IQTREE(2) program

(1):Ziheng Yang; PAML 4: Phylogenetic Analysis by Maximum Likelihood, Molecular Biology and Evolution,
    Volume 24, Issue 8, 1 August 2007, Pages 1586–1591, https://doi.org/10.1093/molbev/msm088
(2):L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective
    stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274.
    https://doi.org/10.1093/molbev/msu300


# Detailled instructions on script usage
### PAML-DATASET

#### Preparation
1. Make sure you have the following information handy:
    - initial tree with internal node numbers (can be extracted from PAML rst file on line 15: ```head -n 15 rst | tail -n 1 >> tree_with_nodes.tree```)
    - length of the full alignment used for PAML calculations (# of positions)

### Script Usage
2. Run the script to generate full length fasta sequence for the desired node with the following command:
    ```
    ./asr-pipeline_paml.py -r <rst_file>  -n <node_number> -a <alignment_length>
    ```
    This will generate three files:
    - full length fasta sequence
    -  plots with posterior probability per reside as well as histogram of posterior probabilites
    -  file with all posterior probabilities, one per line (can be used to plot posterior probabilites on structure)

3. Use the generated full length sequence and determine which areas to cut (based on child and parent cladse). Store positions to cut in the following format:  per line two numbers separated by comma. The left number is included and right number exluded from trim (e.g. 1,50 will cut residue 1 through 49)

4. Re-run script, this time specifying the file with the positions to be trimmed:
    ```
    ./asr-pipeline_paml.py -r <rst_file>  -n <node_number> -a <alignment_length> -c <file_with_pos_to_cut>
    ```
    This will generate six files:
    - full length fasta sequence
    - the trimmed fasta sequence
    - plots with posterior probability per reside as well as histogram of posterior probabilites for full length
    - file with all posterior probabilities, one per line (can be used to plot posterior probabilites on structure) for full length trimmed fasta sequence
    - plots with posterior probability per residue as well as histogram of posterior probabilities for trimmed sequence
    - file with all posterior probabilities, one per line (can be used to plot posterior probabilites on structure) for trimmed sequence

The pdf files with plots will be named according to node as well as average posterior probability, rounded to two decimals.
