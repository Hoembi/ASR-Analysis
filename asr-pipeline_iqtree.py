#! /usr/bin/env python

############################################################
#   Script to process data from IQ-Tree  and calculate     #
#   ancestral sequences with their respective posteriour   #
#   probabilities                                          #
#                                                          #
#                                                          #
# 01/30/2019: initial release                              #
#                                                          #
#                   Marc Hoemberger,2019                   #
#         email: marc.hoemberger[at]gmail.com              #
############################################################


def read_data(node_number, state_file, cut_file):
    """
    Read in the IQ-TREE state file that contains information about ancestral nodes as
    well as the file that defines positions to be cut.
    """
    with open(state_file, "r") as infile:
        input_data = []
        for line in infile:
            if line[0] == "#" or line[5] == "S":
                pass
            elif line.split("\t")[0] == "Node" + str(node_number):
                temp = []
                for item in line.split("\t"):
                    temp.append(item)
                input_data.append(temp)

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

    return to_remove, input_data


def generate_fastas(input_data, to_remove, node_number):
    """
    Generates both full lenght and trimmed fasta files
    """
    fasta = ""
    for i in range(len(input_data)):
        fasta = fasta + input_data[i][2]

    with open("ANC-N" + str(node_number) + "_full_length.fas", "w") as outfile:
        outfile.write(">ANC-N" + str(node_number) + "\n")
        outfile.write(fasta)

    if to_remove != None:
        fasta2 = ""
        for item in input_data:
            if int(item[1]) not in to_remove:
                fasta2 = fasta2 + item[2]

        with open("ANC-N" + str(node_number) + "_trimmed.fas", "w") as outfile:
            outfile.write(">ANC-N" + str(node_number) + "_trimmed" + "\n")
            outfile.write(fasta2)


def calc_statistics(input_data, to_remove, node_number):
    """
    Calculates the mean posterior probabilities per residue position for both
    full length and trimmed sequences.
    """
    import numpy
    from future.builtins import range

    stat_all = []
    for item in input_data:
        stat = []
        for j in range(3, 23):
            stat.append(float(item[j]))
        stat_all.append(stat)

    full_pp = []
    for item in stat_all:
        full_pp.append(max(item))

    print("average pp for full length - no trimming done: ")
    print(numpy.mean(full_pp))

    full_statistics = []
    for i in list(range(len(full_pp))):
        full_statistics.append([i + 1, full_pp[i]])

    if to_remove != None:
        cut_statistics = []
        cut_pp = []
        for i in list(range(len(full_pp))):
            if (i + 1) not in to_remove:
                cut_pp.append(full_pp[i])
                cut_statistics.append([i + 1, full_pp[i]])

        print("average pp for trimmed sequence:")
        print(numpy.mean(cut_pp))

    else:
        cut_statistics = None
        cut_pp = None



    return full_statistics, cut_statistics, cut_pp, full_pp

def plot_statistics(full_statistics, cut_statistics, node_number):
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
    for item in full_statistics:
        ys_full.append(item[1])

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
        "Node"
        + str(node_number)
        + "_full_"
        + str(round(numpy.mean(ys_full), 2))
        + ".pdf",
        transparent=True,
    )
    if cut_statistics != None:
        ys_cut = []
        for item in cut_statistics:
            ys_cut.append(item[1])
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
            "Node"
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

def execute_pipeline(state_file, cut_file, node_number, bfactor):
    to_remove, input_data = read_data(
        node_number, state_file, cut_file
    )
    generate_fastas(input_data, to_remove, node_number)
    full_statistics, cut_statistics, cut_pp, full_pp = calc_statistics(
        input_data, to_remove, node_number
    )
    plot_statistics(full_statistics, cut_statistics, node_number)

    if bfactor == 'yes':
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
        "      Usage: ./asr_pipeline_iqtree.py -r <rst-file> -c <file with positions to cut> -n <node number> -b <b-factor desired> \n\n"
        "         -r (--state_file)        : name of the input file directly from PAML, usually named rst \n"
        "         -n (--node_number)     : node number of the ancestral node that is to be calculated \n"
        "         -c (--cut_file)        : [OPTIONAL] file that define positions to trim. Format:  two, comma separated numbers per line .\n"
        "                                  The left number is included and right number excluded from trim (e.g. 1,50 will cut residue 1 through 49)\n"
        "         -b (--bfactor)         : [OPTIONAL] write PP into B-factor file for plotting on the structure (yes or no; Default: yes) \n"
    )

    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hr:c:n:b:",
            [
                "help",
                "rst_file=",
                "cut_file=",
                "node_number=",
                "bfactor=",
            ],
        )
    except getopt.error:
        print(help)
        sys.exit(2)

    if len(opts) < 2:
        print(help)
        sys.exit(0)

    bfactor = "yes"
    cut_file = None
    for option, argument in opts:
        if option in ("-h", "--help"):
            print(help)
            sys.exit(0)
        elif option in ("-r", "--state_file"):
            input_data = argument
        elif option in ("-c", "--cut_file"):
            cut_file = argument
        elif option in ("-n", "--node_number"):
            node_number = int(argument)
        elif option in ("-b", "--bfactor"):
            bfactor = argument.lower()

    execute_pipeline(input_data, cut_file, node_number, bfactor)


if __name__ == "__main__":
    main()
