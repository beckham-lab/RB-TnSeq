#!/usr/bin/python3
import numpy as np
import sys, os
import pandas as pd

# This function reads through all the old locus tags (PP_xxx format) in # Locus_Tag column of the _Summary output file from the 8_Fitness_Compare.py 
#	function and matches it with the appropriate gene name, new locus tag, and descriptor for qualified genes, as found in the updated KT2440 annotation 
#	file (GCF_000007565.2_ASM756v2_feature_table). If no match is found, the function fills in blank spaces at the corresponding location in the new output file.
#
# Note, this script muct obviously be modified for different organisms.
#
# Usage: 
# python3 PATH/9_Summary_annotate.py <PATH/C1_v_C2_Summary.csv> <PATH/feature_table.txt> <output_filename> 
#
#
#	Input Files:
#		<PATH/C1_v_C2_Summary.csv> The statistics summary file generated by the 8_fitness_compare.py function
#		<PATH/feature_table.txt> The annotated organisms feature table for KT2440, with updated locus tag IDs and gene descriptions
#
#
# 	Output File:
#		<output_filename>.csv is a csv file that modifies the _summary.csv file generated from the 6_Fitness_Compare function to include four identifier
#			columns: ‘old_locus_tag’, ‘new_locus_tag’, ‘gene_name’, and ‘description’
#			Note, the gene name is in the E. coli format (ilvA, ect.), if one exists. 
#			Gene_name and Description columns are left blank if none are provided in the feature table for the organism.


#input arguments:
Summary = sys.argv[1]
genesTable = sys.argv[2]
out_File = sys.argv[3]

# import gene numbers and and all the rest from the Summary Table
Summary_locus = np.loadtxt(Summary, delimiter=',',skiprows=1, usecols=(0), dtype=str)
Summary_rest = np.loadtxt(Summary, delimiter=',',skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10,11,12), dtype=float)

#read the labels into an array
pdSum = pd.read_csv(Summary, header=None)
npSum = np.array(pdSum, dtype=str)
Summary_rest_labels = npSum[0,1:]

#identify new header string array
first_columns = np.array(['old_locus_tag', 'new_locus_tag', 'gene_name','description'])
all_labels = np.concatenate((first_columns,Summary_rest_labels),0)

# import gene numbers and descriptions from table of genes (same as genes.GC file used in BarSeqR.pl)
featurelocus = np.loadtxt(genesTable, delimiter='\t', skiprows=1, usecols=(18), dtype=str)
Newfeaturelocus = np.loadtxt(genesTable, delimiter='\t', skiprows=1, usecols=(15), dtype=str)
featureGene = np.loadtxt(genesTable, delimiter='\t', skiprows=1, usecols=(14), dtype=str)
featureDesc = np.loadtxt(genesTable, delimiter='\t', skiprows=1, usecols=(13), dtype=str)


# read through all the (old) locus tags in Summary_locus and match the with the appropriate gene name, new locus tag, and descriptor for qualified genes
# if no match is found, fill in with blank spaces after the locus tag ID

old_tag = []
fullGene =[]
newLocus = []
fullDesc = []

for i in range(0,len(Summary_locus)):
	if Summary_locus[i] in featurelocus:
		for j in range(0,len(featurelocus)):
			if featurelocus[j]==Summary_locus[i]:
				old_tag.append(Summary_locus[i])
				fullGene.append(featureGene[j])
				newLocus.append(Newfeaturelocus[j])
				fullDesc.append(featureDesc[j])
				break
	else:
		old_tag.append(Summary_locus[i])
		fullGene.append('')
		newLocus.append('')
		fullDesc.append('')

# build new numpy array
Arr_begin = np.transpose(np.array([old_tag, newLocus, fullGene, fullDesc]))

New_array = np.hstack((Arr_begin,Summary_rest))

Full_Array = np.vstack((all_labels,New_array))

# write file with everything

np.savetxt((out_File), Full_Array, delimiter=',', fmt='%s')
