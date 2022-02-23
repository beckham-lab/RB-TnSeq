#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys, os

#	This function takes the _<30_Unused_BL_genes.csv files from the three replicates and merges them into one file, eliminating any duplicates from the list.
# 		This is an essential step prior to generating the replicates table, since if a gene doesn't have a fitness value for all replicates, this file is used to trim 
#		it from the dataset.

#	Usage: 
#		python3 PATH/6_Combined_30_Count_Replicates_List.py <PATH_input_output_Files> <output_filename> <replicate_A_<30_Unused_BL_genes.csv> <replicate_B_<30_Unused_BL_genes.csv> <replicate_C_<30_Unused_BL_genes.csv>
#
#		NOTE, the '<' must be escaped in Bash scripting using a '\'
#
#   Input_Files:
#		<replicate_A_<30_Unused_BL_genes.csv> The genes from replicate A that didn't meet the >30 counts in a gene cutoff condition
#		<replicate_B_<30_Unused_BL_genes.csv> The genes from replicate B that didn't meet the >30 counts in a gene cutoff condition
#		<replicate_C_<30_Unused_BL_genes.csv> The genes from replicate C that didn't meet the >30 counts in a gene cutoff condition
#	
#	Output file:
#		<output_filename> is a csv file with a list of all the unused (<30 count) genes from each of the baseline replicates, where duplicate names are deleted from the list.


# sysnames
In_out_path = sys.argv[1]
out_file = sys.argv[2]
RepA_less_30 = sys.argv[3]
RepB_less_30 = sys.argv[4]
RepC_less_30 = sys.argv[5]

#make sure you are in input/output directory
os.chdir(In_out_path)

#Import the lists from csv files and combine the resultant locsId lists into one array each
pdA = pd.read_csv(RepA_less_30)
A = pdA.loc[:, 'locusId']

pdB = pd.read_csv(RepB_less_30)
B = pdB.loc[:, 'locusId']

pdC = pd.read_csv(RepC_less_30)
C = pdC.loc[:, 'locusId']


#Append all lists together
All = []
All = np.append(All, A)
All = np.append(All, B)
All = np.append(All, C)

#Make a list of strings from the set of unique strings in the list

All_trimmed=list(set(All))
All_trimmed.sort()

#Save as csv in the input/output directory with the given output file name

np.savetxt(out_file, All_trimmed, delimiter=',', fmt='%s')




