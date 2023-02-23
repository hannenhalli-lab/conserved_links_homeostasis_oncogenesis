#!/usr/bin/python

# using ast.literal_eval()
import sys
import ast
import os
import time
import pandas as pd
import pickle

# Step 2
with open(sys.argv[1], 'rb') as config_dictionary_file:
 
    # Step 3
    paths = pickle.load(config_dictionary_file)


t0 = time.time()
current_geneset_list = []
with open(sys.argv[2]) as IN:
    for line in IN:
        current_geneset_list.append(line.rstrip().split("\t")[0])
geneset_name = sys.argv[3]

genes = []
fractions7 = []
fractions6 = []
fractions5 = []
fractions4 = []
fractions3 = []
fractions2 = []
fractions1 = []
fractions0 = []
for gene,sp_paths in paths.items():
    genes.append(gene)
    sp_keys = list(sp_paths.keys())
    gene_subset = { key:value for key,value in sp_paths.items() if key in current_geneset_list}
    frac_less7 = len([key for key,value in gene_subset.items() if value < 7])/len(sp_keys)
    fractions7.append(frac_less7)
    frac_less6 = len([key for key,value in gene_subset.items() if value < 6])/len(sp_keys)
    fractions6.append(frac_less6)
    frac_less5 = len([key for key,value in gene_subset.items() if value < 5])/len(sp_keys)
    fractions5.append(frac_less5)
    frac_less4 = len([key for key,value in gene_subset.items() if value < 4])/len(sp_keys)
    fractions4.append(frac_less4)
    frac_less3 = len([key for key,value in gene_subset.items() if value < 3])/len(sp_keys)
    fractions3.append(frac_less3)
    frac_less2 = len([key for key,value in gene_subset.items() if value < 2])/len(sp_keys)
    fractions2.append(frac_less2)
    frac1 = len([key for key,value in gene_subset.items() if value == 1])/len(sp_keys)
    fractions1.append(frac1)
    frac0 = len([key for key,value in gene_subset.items() if value == 0])/len(sp_keys)
    fractions0.append(frac0)

fractions_lst = [genes,fractions7,fractions6,fractions5,fractions4,fractions3,fractions2,fractions1,fractions0]
df = pd.DataFrame(fractions_lst, dtype = float)

print(df.head)

df_transposed = df.T
df_transposed.columns = ["gene_"+geneset_name,
                         geneset_name+"_Fraction_sp_less7",
                         geneset_name+"_Fraction_sp_less6",
                         geneset_name+"_Fraction_sp_less5",
                         geneset_name+"_Fraction_sp_less4",
                         geneset_name+"_Fraction_sp_less3",
                         geneset_name+"_Fraction_sp_less2",
                         geneset_name+"_Fraction_sp_1",
                         geneset_name+"_Fraction_sp_0"]

file_out = sys.argv[4] + sys.argv[5] + "_fractions_" + geneset_name + ".csv"
df_transposed.to_csv(file_out, index=False, header=True) 

t1 = time.time()
total = t1-t0
print(total)
