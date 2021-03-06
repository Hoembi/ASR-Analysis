This repository contains scripts regularly used for ancestral sequence resurrection.

PAML-DATA:
the script 'asr-pipeline_paml.py' can be used to process data calculated with PAML's codeml(1).

IQTREE-DATA:
the script asr-pipeline_iqtree.py can be used to process data calculated with the IQTREE(2) program

(1):Ziheng Yang; PAML 4: Phylogenetic Analysis by Maximum Likelihood, Molecular Biology and Evolution, Volume 24, Issue 8, 1 August 2007, Pages 1586–1591, https://doi.org/10.1093/molbev/msm088

(2):L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274.
    https://doi.org/10.1093/molbev/msu300


# Detailled instructions on script usage
### PAML-DATASET

### Preparation
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

3. Use the generated full length sequence and determine which areas to cut (based on child and parent cladse). Store positions to cut in the following format:  per line two numbers separated by comma defining the stretch to be trimmed (e.g. 1,20 will cut residues from 1 to 20 and the resulting sequence will start at residue 21, if only one position is to be cut define it as follows: 1,1 --> position 1 will be cut and trimmed sequence starts at 2)

4. Re-run script, this time specifying the file with the positions to be trimmed:
    ```
    ./asr-pipeline_paml.py -r <rst_file>  -n <node_number> -a <alignment_length> -c <file_with_pos_to_cut>
    ```
    This will generate seven files:
    - full length fasta sequence
    - the trimmed fasta sequence
    - the sequence file with removed positons replaced with "-" so that easy double checking of the alignment can be done
    - plots with posterior probability per reside as well as histogram of posterior probabilites for full length
    - file with all posterior probabilities, one per line (can be used to plot posterior probabilites on structure) for full length trimmed fasta sequence
    - plots with posterior probability per residue as well as histogram of posterior probabilities for trimmed sequence
    - file with all posterior probabilities, one per line (can be used to plot posterior probabilites on structure) for trimmed sequence

The pdf files with plots will be named according to node as well as average posterior probability, rounded to two decimals.

### Plotting posterior probabilities on an available structure
5.  Run the ```loadBfacts_recolor.py``` script in pymol. This will add the command loadBfacts_recolor to PYMOL. The script is an adapted version from Pietro Gatti-Lafranconi's script (Gatti-Lafranconi, P.. Pymol script: loadBfacts.py. (2014). doi:10.6084/m9.figshare.1176991.v1).
    Usage of the script:
```
    Replaces B-factors with a list of values contained in a plain txt file

    usage: loadBfactor mol, [startaa, [source, [visual]]]

    mol = any object selection (within one single object though)
    startaa = number of first amino acid in 'new B-factors' file (default=1)
    source = name of the file containing new B-factor values (default=newBfactors.txt)
    visual = redraws and colors the structre and displays bar with min/max values (default=Y)
```
    Be aware that you might have to adjust the Bfactor file to match the sequence stretch in the structure file.
