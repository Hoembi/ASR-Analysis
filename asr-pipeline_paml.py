#! /usr/bin/env python

############################################################
#   Script to process data from PAML  and calculate        #
#   ancestral sequences with their respective posteriour   #
#   probabilities                                          #
#                                                          #
#                                                          #
# 01/30/2019: initial release                              #
#                                                          #
#                   Marc Hoemberger, 2019                  #
#         email: marc.hoemberger[at]gmail.com              #
############################################################


def read_data(rst_file, cut_file, node_number, alignment_length):
    """
    Read in the PAML rst file that contains information about ancestral nodes as
    well as the file that defines positions to be cut.
    """
    import re
    import pandas

    column_name = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]

    temp_list = []
    countdown = 0
    with open(rst_file) as infile:
        for line in infile:
            if countdown > 0 and countdown < 4:
                countdown = countdown + 1
            elif countdown >= 4 and countdown < (alignment_length + 4):
                temp_list.append(line.split()[3:])
                countdown = countdown + 1
            elif countdown > (alignment_length + 4):
                countdown = 0
            elif "Prob distribution at node " + str(node_number) + ", by site" in line:
                countdown = 1

    node_pp = []
    for item in temp_list:
        res_pp = []
        for item2 in item:
            item2 = item2.replace("(", "").replace(")", "")
            item2 = float(re.sub("[A-Z]", "", item2))
            res_pp.append(item2)
        node_pp.append(res_pp)

    position = []
    for i in range(len(node_pp)):
        position.append(i + 1)

    input_data = pandas.DataFrame(node_pp)
    input_data.columns = column_name
    input_data.insert(0, "Position", position)
    input_data = input_data.set_index("Position")

    if cut_file != None:
        first = []
        second = []
        with open(cut_file, "r") as infile:
            for i in infile:
                a = i.strip("\n").split(",")
                first.append(int(a[0]))
                second.append(int(a[1]))
            exclusions = list(zip(first, second))
            to_remove = []
            for i in range(len(exclusions)):
                a = exclusions[i][0]
                b = exclusions[i][1]
                to_remove = to_remove + list(range(a, b))
    else:
        to_remove = None

    return input_data, to_remove


def generate_fastas(input_data, to_remove, node_number):
    """
    Generates both full lenght and trimmed fasta files
    """
    fasta = ""
    for i, row in input_data.iterrows():
        fasta = fasta + row.idxmax(0)

    with open("ANC-N" + str(node_number) + "_full_length.fas", "w") as outfile:
        outfile.write(">ANC-N" + str(node_number) + "\n")
        outfile.write(fasta)

    if to_remove != None:
        fasta2 = ""
        for i, row in input_data.iterrows():
            if i not in to_remove:
                fasta2 = fasta2 + row.idxmax(0)
            else:
                pass
        with open("ANC-N" + str(node_number) + "_trimmed.fas", "w") as outfile:
            outfile.write(">ANC-N" + str(node_number) + "\n")
            outfile.write(fasta2)


def calc_statistics(input_data, to_remove, node_number):
    """
    Calculates the mean posterior probabilities per residue position for both
    full length and trimmed sequences.
    """
    import numpy

    print("average pp for full length - no trimming done: ")
    full_pp = []
    for item, row in input_data.iterrows():
        full_pp.append(row.max())
    print(numpy.mean(full_pp))

    if to_remove != None:
        cut_data = input_data.drop(to_remove, axis=0, inplace=False)
        print("average pp for trimmed sequence:")
        cut_pp = []
        for item, row in cut_data.iterrows():
            cut_pp.append(row.max())
        print(numpy.mean(cut_pp))
    else:
        cut_pp = None

    return (full_pp, cut_pp)


def plot_statistics(full_pp, cut_pp, node_number):
    """
    Plot graphs for posterior probabilities per residue site as histogramm
    as well as line plot.
    """
    from matplotlib import pyplot as plt
    import matplotlib as mpl

    mpl.rcParams.update({"font.size": 12})
    import seaborn as sb
    import numpy

    ys_full = []
    for item in full_pp:
        ys_full.append(item)

    fig1 = plt.figure(figsize=(6, 4))
    sub1 = fig1.add_subplot(121)
    sub1.set_title("PP Node " + str(node_number))
    sub1.set_xlabel("residue number for full length")
    sub1.set_ylabel("posterior probabilities (%)")
    sub1 = plt.plot(ys_full, markersize=1, marker=".", mec="black", lw=0.5)
    sub1 = fig1.add_subplot(122)
    sub1.yaxis.set_label_position("right")
    sub1.yaxis.tick_right()
    sub1.set_title("Node " + str(node_number))
    sub1.set_xlabel("PP")
    sub1.set_ylabel("count")
    sub1 = sb.distplot(ys_full, kde=False, bins=25)
    fig1.savefig(
        "ANC-N"
        + str(node_number)
        + "_full_length_"
        + str(round(numpy.mean(ys_full), 2))
        + ".pdf",
        transparent=True,
    )
    if cut_pp != None:
        ys_cut = []
        for item in cut_pp:
            ys_cut.append(item)
        fig2 = plt.figure(figsize=(6, 4))
        sub2 = fig2.add_subplot(121)
        sub2.set_title("PP Node " + str(node_number))
        sub2.set_xlabel("residue number after trimming")
        sub2.set_ylabel("posterior probabilities (%)")
        sub2 = plt.plot(ys_cut, markersize=1, marker=".", mec="black", lw=0.5)
        sub2 = fig2.add_subplot(122)
        sub2.yaxis.set_label_position("right")
        sub2.yaxis.tick_right()
        sub2.set_title("Node " + str(node_number))
        sub2.set_xlabel("PP")
        sub2.set_ylabel("count")
        sub2 = sb.distplot(ys_cut, kde=False, bins=25)
        fig2.savefig(
            "ANC-N"
            + str(node_number)
            + "_trimmed_"
            + str(round(numpy.mean(ys_cut), 2))
            + ".pdf",
            transparent=True,
        )


def color_bfactor(cut_pp, full_pp, node_number):
    """
    Generate text file with posterior probabilities for plotting on structure
    as B-factors with PYMOL script.
    """
    with open("ANC-N" + str(node_number) + "_Bfactors_full-length.dat", "w") as outfile:
        for i in range(len(full_pp)):
            outfile.write(str(full_pp[i]) + "\n")
    if cut_pp != None:
        with open("ANC-N" + str(node_number) + "_Bfactors_trimmed.dat", "w") as outfile:
            for i in range(len(cut_pp)):
                outfile.write(str(cut_pp[i]) + "\n")


def execute_pipeline(rst_file, cut_file, node_number, length, bfactor):
    input_data, to_remove = read_data(rst_file, cut_file, node_number, length)
    generate_fastas(input_data, to_remove, node_number)
    full_pp, cut_pp = calc_statistics(input_data, to_remove, node_number)
    plot_statistics(full_pp, cut_pp, node_number)
    if bfactor == "yes":
        color_bfactor(cut_pp, full_pp, node_number)


def main():
    import sys
    import getopt

    help = (
        "\n*******************************************************************************************************************************\n"
        "**         ASR-Pipeline script -- python script to calculate ancestral sequences with their respective posterior             ** \n"
        "**                                probabilities (PP). Gives sequences for ancestor in fasta file, plots with PP and          ** \n"
        "**                                B-factor file that can be used to plot PP on the structure.                                **\n"
        "*******************************************************************************************************************************\n\n"
        "      Usage: ./asr_pipeline_paml.py -r <rst-file> -c <file with positions to cut> -n <node number> -b <b-factor desired> \n\n"
        "         -r (--rst_file)        : name of the input file directly from PAML, usually named rst \n"
        "         -n (--node_number)     : node number of the ancestral node that is to be calculated \n"
        "         -a (--alignment_length): length of the input alignment (# of positions) \n"
        "         -c (--cut_file)        : [OPTIONAL] file that define positions to trim. Format:  per line two numbers separated by comma. \n"
        "                                  The left number is included and right number exluded from trim (e.g. 1,50 will cut residue 1 through 49)\n"
        "         -b (--bfactor)         : [OPTIONAL] write PP into B-factor file for plotting on the structure (yes or no; Default: yes) \n"
    )

    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hr:c:n:a:b:",
            [
                "help",
                "rst_file=",
                "cut_file=",
                "node_number=",
                "alignment_length",
                "bfactor=",
            ],
        )
    except getopt.error:
        print(help)
        sys.exit(2)

    if len(opts) < 3:
        print(help)
        sys.exit(0)

    bfactor = "yes"
    cut_file = None
    for option, argument in opts:
        if option in ("-h", "--help"):
            print(help)
            sys.exit(0)
        elif option in ("-r", "--rst_file"):
            input_data = argument
        elif option in ("-c", "--cut_file"):
            cut_file = argument
        elif option in ("-n", "--node_number"):
            node_number = int(argument)
        elif option in ("-a", "--alignment_length"):
            alignment_length = int(argument)
        elif option in ("-b", "--bfactor"):
            bfactor = argument.lower()

    execute_pipeline(input_data, cut_file, node_number, alignment_length, bfactor)


if __name__ == "__main__":
    main()
