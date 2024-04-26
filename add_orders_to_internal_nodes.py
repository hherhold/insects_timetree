#!/usr/bin/env python

# add orders to the internal nodes of the tree


import sys
import os
import re
import argparse
import taxidTools
import ete3

# parse command line arguments. Get the input tree file and the output file.
# We also need the directory where the taxdump files are located.
# Input tree should have a -i flag.

parser = argparse.ArgumentParser(description='Add orders to internal nodes of a tree')
parser.add_argument('-i', '--input', help='The input tree file')
parser.add_argument('-o', '--output', help='The output file')
parser.add_argument('-t', '--taxdump', help='The directory where the taxdump files are located')
args = parser.parse_args()

def main():
    print("Loading NCBI taxonomy data...", end="")
    # Flush the output buffer.
    sys.stdout.flush()
    tax = taxidTools.Taxonomy.from_taxdump(f"{args.taxdump}/nodes.dmp", \
                                           f"{args.taxdump}/rankedlineage.dmp")
    print("done.")
    lin = tax.getAncestry(tax.getTaxid('Syrphidae'))
    lin.filter(['order'])
    print(lin[0].name)

if __name__ == "__main__":
    main()

    