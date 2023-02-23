#!/usr/bin/python

import sys
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import random
import statistics as stats
from scipy.stats import wilcoxon
import openpyxl
from pathlib import Path

G = nx.Graph()
#  give filename
file2 = "STRING_db_list.txt"

# create dictionary
string_db = {}
with open(file2) as IN:
    for line in IN:
        ensembl = line.rstrip().split("\t")[0]
        gene_name = line.rstrip().split("\t")[1]

        string_db[ensembl] = gene_name

#  give filename
file3 = "string_network_scores.txt"

with open(file3) as IN:
    for line in IN:
        node1 = line.rstrip().split("\t")[0]
        node2 = line.rstrip().split("\t")[1]
        score = int(line.rstrip().split("\t")[2])

        node1_translated = string_db[node1]
        node2_translated = string_db[node2]

        G.add_edge(node1_translated, node2_translated, weight=score)


# +
# create list
file2 = "common_nodes_hPPIN_STRING_PANCAN.txt"
combined_list = []
with open(file2) as IN:
    for line in IN:
        combined_list.append(line.rstrip().split("\t")[0])
        
removal_list = [node for node in G.nodes if node not in combined_list]
G.remove_nodes_from(removal_list)
# -

path_lengths = dict(nx.all_pairs_shortest_path_length(G))

with open("STRING_sp_lengths_total.txt", 'w') as f: 
    for key, value in path_lengths.items(): 
        f.write('%s\t%s\n' % (key, value))
