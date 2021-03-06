#!/usr/bin/python3
from numpy import *
import numpy as np
import pandas as pd
import sys, os

# This function loads your all.poolcount file, splits into groups based upon biological replicate, groups insertions according to gene locus, and eliminates strains that have 
# fewer than three insertions in the baseline sample, according to the baseline used. 

# NOTE- This function breaks if using less than or more than 3 biological replicates. Modify accordingly.

# Usage:
#   python3 PATH/4_BarSeqProc_loadExps.py <PATH/all.poolcount_file> <baseline_condition> <Test Conditions>

# Inputs:
#   all.poolcount file is file generated by BarSeqR.pl function 
#   baseline condition chosen as a reference (1, 2, 3, etc.)
#   Test Conditions Chosen '4,5,6,etc.''   NOTE. This is the set number as it apears in the .poolcount file. Not the neccisarily particular condition index.
#		So third set in would be '3'. Seperate condition set numbers with commas e.g. '2,3,4,5'

# Output files:
#   <baseline condition><replicates>_counts-in-genes.csv—Files made from extracting rows from the all.poolcount file where the barcoded transposon insertion
#		falls within a gene and breaking apart by replicates A-C 
#   <baseline condition><replicates>_in3genes.csv—File made from extracting rows from the <baseline condition><replicate>_counts-in-genes.csv file, where >= 3
#		counts for each distinct transposon insertion exists in the T=0 baseline condition. This breaks apart into sets based upon replicate A,B, or C


if len(sys.argv) < 3:
	print("Run this in the directory where you want your output files.")
	print("Usage: " + sys.argv[0] + " <path to all.poolcount> <baseline condition> <replicate set>")
	print("Arg 1: e.g., /path/to/all.poolcount/")
	print("Arg 2: baseline condition (1,2,3,etc.)")
	print("Arg 3: Test Conditions to compare to the baseline")
	sys.exit(0)

if not os.path.isfile(sys.argv[1]):
	print (sys.argv[1] + " does not exist... aborting...")
	sys.exit(0)

# sysnames
poolFile = sys.argv[1]
baseline = sys.argv[2]
baseline = np.array(baseline)
Test_Conditions = sys.argv[3]
Test_Conditions = np.array(Test_Conditions.split(',')) #needed to change the string to an array and ignore the commas in the argv input


#################################################################################
##### Split all.poolcount data into the three biological replicates #####
################################################################################# 

# Use Pandas to make import all.poolcount as a dataframe
Counts_Table = pd.read_csv(poolFile, sep="\t")

# Determine the total number of condition columns in the master dataframe
Total_Columns = list(range(7,Counts_Table.shape[1],1))

# Dataframe For Replicate A, B, and C
drop_indexA = list(np.delete(Total_Columns, slice(0,len(Total_Columns), 3))) #Indexes to drop to create replicate dataframes
drop_indexB = list(np.delete(Total_Columns, slice(1,len(Total_Columns), 3))) 
drop_indexC = list(np.delete(Total_Columns, slice(2,len(Total_Columns), 3)))

drop_colsA = [j for i,j in enumerate(Counts_Table.columns) if i in drop_indexA] #fetch the column names at the indexes you wish to drop from the datasets
drop_colsB = [j for i,j in enumerate(Counts_Table.columns) if i in drop_indexB]
drop_colsC = [j for i,j in enumerate(Counts_Table.columns) if i in drop_indexC] 

Counts_Table_A = Counts_Table.drop(drop_colsA, axis=1) # Drop the other two replicate columns from the master dataframe
Counts_Table_B = Counts_Table.drop(drop_colsB, axis=1)
Counts_Table_C = Counts_Table.drop(drop_colsC, axis=1)


################################################################################
# Determine the Baseline Condition (to be placed in the first column of count  #
# data in each dataframe) and the enrichment condition indexes that will be    #
# compared to the designated baseline                                          #
################################################################################


BL_Enrichment_Col_Indexes = str(int(baseline)-1) #Taking away 1 accounts for python indexing
for i in range(0,len(Test_Conditions)):
    BL_Enrichment_Col_Indexes = np.append(BL_Enrichment_Col_Indexes,str(int(Test_Conditions[i])-1)) #Taking away 1 accounts for python indexing
BL_Enrichment_Col_Indexes = BL_Enrichment_Col_Indexes.astype(int)

# First determine the total number of condition columns in the replicate dataframes (all should be the same)
Total_Conditions_Columns = list(range(7,Counts_Table_A.shape[1],1))
drop_index_Conditions = list(np.delete(Total_Conditions_Columns, BL_Enrichment_Col_Indexes)) #Indexes of conditions to drop from each Replicate dataframe

# Eliminate the condition columns not needed for the current analysis each Replicate Dataframe
drop_Cond_colsA = [j for i,j in enumerate(Counts_Table_A.columns) if i in drop_index_Conditions] 
Table_A = Counts_Table_A.drop(drop_Cond_colsA, axis=1)

drop_Cond_colsB = [j for i,j in enumerate(Counts_Table_B.columns) if i in drop_index_Conditions] 
Table_B = Counts_Table_B.drop(drop_Cond_colsB, axis=1)

drop_Cond_colsC = [j for i,j in enumerate(Counts_Table_C.columns) if i in drop_index_Conditions] 
Table_C = Counts_Table_C.drop(drop_Cond_colsC, axis=1)


###################################################################################
# Trim Replicate Tables to only include rows where an insertion is within a gene  #
# and then trim to only include rows where the mutant count for each unqiue       #
# transposon is at least 3 in the baseline                                        #
###################################################################################

#select only rows that are within genes (ingenes); write output to csv    
ingenes_A = Table_A.dropna(subset=['locusId'])
ingenes_B = Table_B.dropna(subset=['locusId'])
ingenes_C = Table_C.dropna(subset=['locusId'])

ingenes_A.to_csv(('BL_'+str(baseline)+'A_counts-in-genes.csv'), index=False, header=True)
ingenes_B.to_csv(('BL_'+str(baseline)+'B_counts-in-genes.csv'), index=False, header=True)
ingenes_C.to_csv(('BL_'+str(baseline)+'C_counts-in-genes.csv'), index=False, header=True)
      
# select only rows that have >= 3 reads/strain in T0 (in3genes); write output to csv
colname_A = ingenes_A.columns[7] # Need to identify the column names for the baseline columns in the replicate tables
colname_B = ingenes_B.columns[7]
colname_C = ingenes_C.columns[7]

in3genes_A = ingenes_A[ingenes_A[colname_A] >= 3]
in3genes_B = ingenes_B[ingenes_B[colname_B] >= 3]
in3genes_C = ingenes_C[ingenes_C[colname_C] >= 3]

in3genes_A.to_csv(('BL_'+str(baseline)+'A_in3genes.csv'), index=False, header=True)
in3genes_B.to_csv(('BL_'+str(baseline)+'B_in3genes.csv'), index=False, header=True)
in3genes_C.to_csv(('BL_'+str(baseline)+'C_in3genes.csv'), index=False, header=True)

