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

# We need to keep track of all the orders we've seen.
orders = set()

def main():
    # First make sure the input file exists. There's no point in loading the
    # taxonomy DB (which is slow) if the input file doesn't exist.
    if not os.path.exists(args.input):
        print(f"Input file {args.input} does not exist.")
        exit(1)

    print("Loading NCBI taxonomy data...", end="")
    sys.stdout.flush()
    tax = taxidTools.Taxonomy.from_taxdump(f"{args.taxdump}/nodes.dmp", \
                                           f"{args.taxdump}/rankedlineage.dmp")
    print("done.")
    #lin = tax.getAncestry(tax.getTaxid('Machilidae'))
    #lin.filter(['order'])
    #print(lin[0].name)
    #exit(0)

    # There are a handful of families that are not in the taxonomy database. Let's make
    # a dictionary of these with their orders. I found these online; Pemphigidae is a
    # family of aphids, and Ascalaphidae is a family of lacewings. Pemphigidae is apparently
    # no longer used (at least according to iNaturalist).
    missing_families = {
        "Ascalaphidae": "Lepidoptera",
        "Pemphigidae": "Hemiptera",
        "Anobiidae": "Coleoptera",
        "Xylophagidae": "Diptera"
    }

    # Load the tree
    t = ete3.Tree(args.input, format=1)
    leaf_count = 0
    for n in t.traverse():
        if n.is_leaf():
            leaf_count += 1
            family_name = n.name.split("-")[0]
            #print(f"Leaf node {n.name}, cleaned family name {family_name}")

            order_name = "NA"

            try:
                family_taxid = tax.getTaxid(family_name)
                lin = tax.getAncestry(tax.getTaxid(family_name))
                lin.filter(['order'])
                order_name = lin[0].name
            except KeyError:
                print(f"Family {family_name} not found in the taxonomy database, using dictionary")
                order_name = missing_families[family_name]

            # Myida is not an insect. What family is this?
            if order_name == "Myida":
                print(f"Family {family_name} is Mydia. What is this?")

            #print("Adding order to leaf node", n.name, lin[0].name)
            n.add_features(order=order_name)

            # Keep track of all the orders we've seen.
            orders.add(order_name)

    print(f"Added orders to {leaf_count} leaf nodes")

    # Let's dump all the orders we've seen to stdout.
    #print("Orders seen:")
    #for o in orders:
    #    print(o)

    print("Searching for monophyletic groups...")
    # For each order in the set of orders we've seen, find the monophyletic group
    # (or groups) and assign the order to the internal node. Also set the name
    # to be the order name.
    for order in orders:
        print("Searching for", order)
        for node in t.get_monophyletic(values=[order], target_attr="order"):
            print(f"Setting {node.name} to {order}")
            node.add_features(order=order)
            node.name = order

    # Finally, all the internal nodes that are just numbers should have their
    # names set to blank, otherwise the tree will look weird with lots of numbers.
    for n in t.traverse():
        if not n.is_leaf() and re.match(r"\d+", n.name):
            n.name = ""

    # All done - write the tree to the output file.
    t.write(outfile=args.output, format=1)

if __name__ == "__main__":
    main()

    