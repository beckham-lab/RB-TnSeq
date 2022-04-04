#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys, os
import matplotlib.pyplot as plt
import seaborn as sns
from distutils.util import strtobool

#This allows you to change text in illustrator for the PDFs exported in this function
plt.rcParams['pdf.fonttype'] = 42
#This changes the sans-serif font output for plots to Arial
plt.rcParams.update({'font.sans-serif':'Arial'})

#####
#Generate a heatmap showing fitness values for the genes listed in the first column of the input file. The conditions to include in the heatmap are listed in the second column of the input file (The glucose only baseline is added automatically to the heatmap)
##### 

#Example call:
 #
 #	#python3 10_heatmap.py <PATH/Included_Genes_List_File.csv> <PATH to Annoted genes file>
 #
 #	Note, all the optional inputs below can be declared in the command line to change the heatmap output as desired.
 #
 #
 #	Input File:
 #		Included_Genes_List_File -- A file where the first column is the a list of genes the user desires to include in the heatmap output. These become row in the heatmap
 #									The second column is a list of conditions files to include in the heatmap output as columns. Note, this abstracts from the Annotated_Summary Files
 #									Pay particular notcie to the code at Lines 71 and 123, since if your input file isn't named appropriately, this function will break
 #									Note, replicates are added individually as a column
 #										There is no functionality to use the mean of a particular condition and present it as a column in the heatmap output 
 #
 #
 #	Output File:
 #		heatmap_'Genes_n_Conditions_Filename'-- A heatmap file where columns are the input conditions (and the parallel medium comparison condition) and rows are the input genes for analysis
 #		
 #

def main(Genes_n_Conditions_File,Input_Fitfile_Dir_Path, font_scale=1.0, method='average', metric='euclidean', robust='True', figheight=8, figwidth=8, cmap='RdBu', center=0, col_cluster='False', row_cluster='True', linewidths=0 , linecolor='black'):

	#set appropriate **kwargs to ints, floats, bools, or arrays
	font_scale=float(font_scale)
	figheight = int(figheight)
	figwidth = int(figwidth)
	figsize=np.array([figwidth,figheight])
	center=float(center)
	robust=strtobool(robust)
	col_cluster=strtobool(col_cluster)
	row_cluster=strtobool(row_cluster)
	linewidths = int(linewidths)

	# load the .csv file with pandas, reading all rows and columns, including headers and convert to np array 
	Whole_File = pd.read_csv(Genes_n_Conditions_File, header=None)
	print(Whole_File)
	Whole_File = Whole_File.fillna('')
	Whole_File = np.array(Whole_File[:])
	Genes = np.array(Whole_File[1:,0],dtype=str)
	Conditions = np.array(Whole_File[1:,1],dtype=str)

#######################################################################################################################################################################################
# For each gene, navigate to each file for the conditions listed in the input file and add the fitness values into an array
#	the array format is conditions and replicates across columns and different genes down rows
#######################################################################################################################################################################################

	#Start building the numpy fitness value array. For this, iterate through each condition file and extract fitness values for the designated genes
	# for the first dfitness file, also extract the glucose alone fitness values
	fitness_array =[]
	condition_counter = 1 	#This is a counter to keep track of whether this is the first fitness file or not
	for i in range(0, len(Conditions)):
		if condition_counter == 1: 
			current_file_name= Input_Fitfile_Dir_Path+'/M9_Glucose_v_M9_'+Conditions[i]+'_Annotated_Summary.csv'
			current_file = pd.read_csv(current_file_name, header=None)
			current_file = np.array(current_file[:])
			genes_list = np.array(current_file[1:,0],dtype=str)
			Description_List = np.array(current_file[1:,3],dtype=str)
			Gene_Names = np.array(current_file[1:,2],dtype=str)
			glucose_fitness = np.array(current_file[1:,4:7], dtype=float)
			condition_fitness = np.array(current_file[1:,8:11], dtype=float)
			Rep_labels = np.append(np.array(current_file[0,4:7], dtype=str),np.array(current_file[0,8:11],dtype=str))
			
			#Build an array for the Condition Labels, extracted from the input file
			Condition_Labels = Rep_labels
			
			#add the matching fitness values for this condition (and glucose) into seperate arrays for each replicate for all the genes listed in the heatmap input file
			GluA = []
			GluB = []
			GluC = []
			CondA = []
			CondB = []
			CondC = []
			NEW_Description=[]
			for j in range(0, len(Genes)):
				if len(Genes[j]) < 1:  #This happens if there are less genes than conditions in the file
					break
				if Genes[j] not in genes_list: #Not all genes may exist in all datasets, if this happens fill with fitness_array with blank cells and move to next gene
					GluA = np.append(GluA, '')
					GluB = np.append(GluB, '')
					GluC = np.append(GluC, '')
					CondA = np.append(CondA, '')
					CondB = np.append(CondB, '')
					CondC = np.append(CondC, '') 
					NEW_Description = np.append(NEW_Description, Genes[j])  #This accounts for if the locus ID and associated description is not present in the fitness files
					continue
				a = np.where(genes_list == Genes[j])
				GluA = np.append(GluA, glucose_fitness[a,0])
				GluB = np.append(GluB, glucose_fitness[a,1])
				GluC = np.append(GluC, glucose_fitness[a,2])
				CondA = np.append(CondA, condition_fitness[a,0])
				CondB = np.append(CondB, condition_fitness[a,1])
				CondC = np.append(CondC, condition_fitness[a,2])
				
				#The following .replace lines are my work around for in the description contains a forward slash in the string and also replaces spaces with underscores
				Gene_Description = str(np.char.replace(Description_List[a],' ','_')) 
				Gene_Description = str(np.char.replace(Gene_Description,'/','--'))
				NEW_Description=np.append(NEW_Description,np.char.replace([Genes[j]+'('+str(Gene_Names[a])+')'+': '+Gene_Description],'--','/')) #The forward slash is added back in for the descriptions that contain them
				
				# I was having issues with ' and [] symbols showing up in the merged strings, so this was my crude work-around
				NEW_Description=np.char.replace(NEW_Description,'[','')
				NEW_Description=np.char.replace(NEW_Description,']','')
				NEW_Description=np.char.replace(NEW_Description,'\'','')
			fitness_array = np.column_stack((GluA,GluB,GluC,CondA,CondB,CondC))
			condition_counter = condition_counter+1
		else:
			if len(Conditions[i]) < 1:  #This happens if there are less conditions than genes in the file
				break
			current_file_name= Input_Fitfile_Dir_Path+'/M9_Glucose_v_M9_'+Conditions[i]+'_Annotated_Summary.csv'
			current_file = pd.read_csv(current_file_name, header=None)
			current_file = np.array(current_file[:])
			genes_list = np.array(current_file[1:,0],dtype=str)
			Description_List = np.array(current_file[1:,3],dtype=str)
			Gene_Names = np.array(current_file[1:,2],dtype=str)
			condition_fitness = np.array(current_file[1:,8:11], dtype=float)
			Rep_labels=[]
			Rep_labels = np.append(Rep_labels,np.array(current_file[0,8:11],dtype=str))
			Condition_Labels = np.append(Condition_Labels,Rep_labels)
			CondA = []
			CondB = []
			CondC = []
			NEW_Description=[]
			for j in range(0, len(Genes)):
				if len(Genes[j]) < 1:  #This happens if there are less genes than conditions in the file
					break
				if Genes[j] not in genes_list: #Not all genes may exist in all datasets, if this happens fill with fitness_array with blank cells and move to next gene
					CondA = np.append(CondA, '')
					CondB = np.append(CondB, '')
					CondC = np.append(CondC, '')
					NEW_Description = np.append(NEW_Description, Genes[j]) # This accounts for if the locus ID and associated description is not present in the fitness files
					continue
				a = np.where(genes_list == Genes[j])
				CondA = np.append(CondA, condition_fitness[a,0])
				CondB = np.append(CondB, condition_fitness[a,1])
				CondC = np.append(CondC, condition_fitness[a,2])

				#The following .replace lines are my work around for in the description contains a forward slash in the string and also replaces spaces with underscores
				Gene_Description = str(np.char.replace(Description_List[a],' ','_')) 
				Gene_Description = str(np.char.replace(Gene_Description,'/','--'))
				NEW_Description=np.append(NEW_Description,np.char.replace([Genes[j]+'('+str(Gene_Names[a])+')'+': '+Gene_Description],'--','/')) #The forward slash is added back in for the descriptions that contain them

				# I was having issues with ' and [] symbols showing up in the merged strings, so this was my crude work-around
				NEW_Description=np.char.replace(NEW_Description,'[','')
				NEW_Description=np.char.replace(NEW_Description,']','')
				NEW_Description=np.char.replace(NEW_Description,'\'','')

			fitness_array = np.column_stack((fitness_array,CondA,CondB,CondC))
			condition_counter = condition_counter+1
	
	#Make a dataframe for all the data with 'Genes' as the row label and 'Condition' as the column label

	Full_array = pd.DataFrame(fitness_array,columns=Condition_Labels,index=NEW_Description)
	Full_array.index.name='Genes'
	
	#Deterimine file name from input file name and save csv table for heatmap
	output_file = 'heatmap_'+Genes_n_Conditions_File.split('/')[-1]
	Full_array.to_csv(output_file)
	
####################################################################################################################################################################################################################
# GENERATE HEATMAP USING SEABORN
####################################################################################################################################################################################################################

	Heatmap_File = pd.read_csv(output_file, na_values="NaN")
	Heatmap_File = Heatmap_File.dropna(how='any') # drops all rows with NaNs, since these can't exist in a heatmap
	Heatmap_File= Heatmap_File.set_index('Genes')
	Heatmap_File.columns.name='Conditions'

	#build the seaborn heatmap (cluster map)
	sns.set(font_scale=font_scale)
	clustered_HM = sns.clustermap(Heatmap_File, method=method, metric=metric, robust=robust, figsize=figsize, cmap=cmap, center=center, col_cluster=col_cluster, row_cluster=row_cluster, linewidths=linewidths , linecolor=linecolor)

	output_file = output_file.split('.csv')[0]+'.pdf'
	plt.savefig(output_file, dpi = 300)

	#, pivot_kws={'index':'Genes', 'columns':'Conditions'}


if __name__=='__main__':
	main(sys.argv[1],
		sys.argv[2],   #required compare file name
		**dict(arg.split('=') for arg in sys.argv[3:])) # take in optional kwargs

