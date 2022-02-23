#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys, os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
from distutils.util import strtobool

#This allows you to change text in illustrator for the PDFs exported in this function
plt.rcParams['pdf.fonttype'] = 42
#This changes the sans-serif font output for plots to Arial
plt.rcParams.update({'font.sans-serif':'Arial'})

#####
#Generate a scatterplot using the comparison summary data from 6A_Gene_Annotate and color any values where the y-condition (test) differs from the x-condition (reference) 
#####

# Usage:
#	python3 PATH/11_2D_Graph.py <PATH/Annotated_Summary.csv> <kwags>
#
# Input File:
#	Annotated_Summary.csv- File generated as output from the 9_Summary_annotate.py function
#
# Output File:
#	'Annotated_Summary_identifier'_<q_val_cutoff>.pdf -- PDF of a scatterplot named according to the name of the Annotated_Summary file used as input,
#		where the q value cutoff used to flag significant genes is also provided in the title

def main(Compare_File, q_val_cutoff=0.1, Show_Sig_Labels='False', Additional_Genes_Labeled='', force_text_val=.2, force_points_val=10, label_font=6, label_color='blue', 
	arrow_color='blue', nonsig_marker = '.', sig_marker='*', diff_between=1.5, high_marker_color='cyan', low_marker_color='r', marker_size='18', 
	yequalsx_line='False', yequalsx_line_style = '-', yequalsx_line_color= 'blue', fitness_diff_line_style = '--', fitness_diff_line_color = 'red'):
	#Break apart the Additional_Genes_Labeled and arrange into a numpy array
	Additional_GeneLabels =[]
	if len(Additional_Genes_Labeled)>0:
		Additional_GeneLabels = np.array(Additional_Genes_Labeled.split(" "), dtype=str)

	#set appropriate **kwargs to ints or floats
	q_val_cutoff=float(q_val_cutoff)
	force_text_val=float(force_text_val)
	force_points_val=int(force_points_val)
	label_font=int(label_font)
	marker_size=int(marker_size)
	diff_between=float(diff_between)
	yequalsx_line=strtobool(yequalsx_line)
	Show_Sig_Labels=strtobool(Show_Sig_Labels)

	#Derive labels for the graph from the Compare_File name
	noPath_label=Compare_File.split("/")[-1]
	Cond1Label=noPath_label.split("_v_")[0]
	long_Cond2Label=noPath_label.split("_v_")[1]
	Cond2Label=long_Cond2Label.split("_Annotated")[0]

	# load the .csv files with pandas, reading all rows and columns, including headers and convert to np array 
	Whole_File = pd.read_csv(Compare_File, header=None)
	Whole_File = np.array(Whole_File[:])

	#Extract out gene IDs (locus tags), mean fitness values, and q-scores
	GenesIDs = np.array(Whole_File[1:,0],dtype=str)
	meansA = np.array(Whole_File[1:,7],dtype=float)
	meansB = np.array(Whole_File[1:,11],dtype=float)
	sort_qVal = np.array(Whole_File[1:,15],dtype=float)

	#plot points where Cond2-Cond1 > 1 in blue
	#plot points where Cond2-Cond1 < -1 in red
	#plot points where 1 > Cond2-Cond1 > -1
	#This populates seperate arrays form performing this during plotting
	Difference = []

	high_test_x_noSig =[]
	high_test_y_noSig =[]
	high_test_x_Sig =[]
	high_test_y_Sig =[]
	high_test_sig_labels =[]

	low_test_x_noSig =[]
	low_test_y_noSig =[]
	low_test_x_Sig =[]
	low_test_y_Sig =[]
	low_test_sig_labels =[]


	equal_test_x_noSig =[]
	equal_test_y_noSig =[]
	equal_test_x_Sig =[]
	equal_test_y_Sig =[]

	for i in range(0, len(GenesIDs)):
		Difference = np.append(Difference, meansB[i] - meansA[i])

		if Difference[i] > diff_between and sort_qVal[i] > q_val_cutoff:
			high_test_x_noSig = np.append(high_test_x_noSig,meansA[i])
			high_test_y_noSig = np.append(high_test_y_noSig,meansB[i])

		if Difference[i] > diff_between and sort_qVal[i] < q_val_cutoff:
			high_test_x_Sig = np.append(high_test_x_Sig,meansA[i])
			high_test_y_Sig = np.append(high_test_y_Sig,meansB[i])
			high_test_sig_labels = np.append(high_test_sig_labels,GenesIDs[i])

		if Difference[i] < -diff_between and sort_qVal[i] > q_val_cutoff:
			low_test_x_noSig = np.append(low_test_x_noSig,meansA[i])
			low_test_y_noSig = np.append(low_test_y_noSig,meansB[i])

		if Difference[i] < -diff_between and sort_qVal[i] < q_val_cutoff:
			low_test_x_Sig = np.append(low_test_x_Sig,meansA[i])
			low_test_y_Sig = np.append(low_test_y_Sig,meansB[i])
			low_test_sig_labels= np.append(low_test_sig_labels,GenesIDs[i])

		if Difference[i] < diff_between and Difference[i] > -diff_between and sort_qVal[i] > q_val_cutoff:
			equal_test_x_noSig = np.append(equal_test_x_noSig,meansA[i])
			equal_test_y_noSig = np.append(equal_test_y_noSig,meansB[i])

		if Difference[i] < diff_between and Difference[i] > -diff_between and sort_qVal[i] < q_val_cutoff:
			equal_test_x_Sig = np.append(equal_test_x_Sig,meansA[i])
			equal_test_y_Sig = np.append(equal_test_y_Sig,meansB[i])

	#plot points where Cond2-Cond1 > diff_between in cyan and not significant with circles
	plt.scatter(high_test_x_noSig, high_test_y_noSig, color=high_marker_color, marker=nonsig_marker, s=marker_size, lw=0)
	#plot points where Cond2-Cond1 > diff_between in cyan and not significant with stars
	plt.scatter(high_test_x_Sig, high_test_y_Sig, color=high_marker_color, marker=sig_marker, s=marker_size, lw=0)
	#plot points where Cond2-Cond1 < -diff_between in red and not significant with circles
	plt.scatter(low_test_x_noSig, low_test_y_noSig, color=low_marker_color, marker=nonsig_marker, s=marker_size, lw=0)
	#plot points where Cond2-Cond1 < -diff_between in red and not significant with stars
	plt.scatter(low_test_x_Sig, low_test_y_Sig, color=low_marker_color, marker=sig_marker, s=marker_size, lw=0)
	#plot points where diff_between > Cond2-Cond1 > -diff_between in black and not significant with circles
	plt.scatter(equal_test_x_noSig, equal_test_y_noSig, color='k', marker=nonsig_marker, s=marker_size, lw=0)
	#plot points where diff_between > Cond2-Cond1 > -diff_between in black and not significant with stars
	plt.scatter(equal_test_x_Sig, equal_test_y_Sig, color='k', marker=sig_marker, s=marker_size, lw=0)

	#Fill in the Title and the axes titles
	plt.title(Cond1Label+' vs '+Cond2Label, fontsize=12)
	plt.xlabel(Cond1Label+' Fitness', fontsize=10)
	plt.ylabel(Cond2Label+' Fitness', fontsize=10)

	#Make a horizontal line at y=o and a vertical line at x=0
	plt.axhline(y=0, color='dimgrey', zorder=0, linestyle='-')
	plt.axvline(x=0, color='dimgrey', zorder=0, linestyle='-')

	#Fine y=x line
	xmin= np.min(plt.xlim())
	xmax= np.max(plt.xlim())

	xline= np.linspace(xmin,xmax,5000)
	yline= xline
	
	#plot a y=x line if asked for it
	if yequalsx_line:
		plt.plot(xline, yline, yequalsx_line_color, alpha=0.5, zorder=1,linestyle=yequalsx_line_style)

	#plot lines where fitness differs from y=x by fitness of y=x+diff_between and y=x-diff_between
	y_lower = np.linspace(xmin-diff_between,xmax-diff_between,5000)
	y_upper = np.linspace(xmin+diff_between,xmax+diff_between,5000)

	plt.plot(xline, y_lower, fitness_diff_line_color, alpha=0.5, zorder=1,linestyle=fitness_diff_line_style)
	plt.plot(xline, y_upper, fitness_diff_line_color, alpha=0.5, zorder=1,linestyle=fitness_diff_line_style)


	####################################################################################################################################################################################################################################################

	##### Plot additional labels with the labels for genes that satisfy the statistical and fold-change cutoffs ######
	LabelIDs = []
	LabelmeansA = []
	LabelmeansB = []
	GeneLabels = []
	if Show_Sig_Labels:
		Sig_GeneLabels=np.append(high_test_sig_labels,low_test_sig_labels)
	else:
		Sig_GeneLabels =[]
	GeneLabels=np.append(Additional_GeneLabels,Sig_GeneLabels)
	if len(GeneLabels)>0:
		for n in range(0,len(GenesIDs)):
			if GenesIDs[n] in GeneLabels: 
				LabelIDs = np.append(LabelIDs, GenesIDs[n])
				LabelmeansA = np.append(LabelmeansA, meansA[n])
				LabelmeansB = np.append(LabelmeansB, meansB[n])

	#plot the labels
	texts=[]
	for i in range(0,len(LabelIDs)):
		texts.append(plt.text(LabelmeansA[i], LabelmeansB[i], LabelIDs[i], weight='bold', c=label_color, fontsize=label_font))

	#Use arrows to point from labels to the corresponding point
	adjust_text(texts,force_text=force_text_val,force_points=force_points_val,arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0', color=arrow_color))
	####################################################################################################################################################################################################################################################

	#Add a legend for the scatterplot
	legend = [str('= q-value < '+str(q_val_cutoff)), str('|x-y| > '+str(diff_between))]

	#Add a legend for significant genes on the scatterplot
	black_star = plt.scatter([], [], color='k', marker=sig_marker, s=marker_size, lw=0)
	diff_line = Line2D([0,1],[0,1],linestyle=fitness_diff_line_style, color=fitness_diff_line_color)
	plt.legend(handles=[black_star, diff_line], labels=legend, loc='lower right', fontsize=8)

	#Save the scatterplot as a PDF
	plt.savefig(Cond1Label+'_v_'+Cond2Label+'_'+str(q_val_cutoff)+'.pdf', dpi = 300)

if __name__=='__main__':
	main(sys.argv[1],   #required compare file name
		**dict(arg.split('=') for arg in sys.argv[2:])) # take in optional kwargs



	