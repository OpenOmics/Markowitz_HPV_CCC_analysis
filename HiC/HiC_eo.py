"""
# HiC_eo.py
# Used to create chromosome-association preference heatmaps for
# Figures 1B and 1C in Warburton et al. (2025).
#
# Input:  SM159.transValidPairs (trans-only allValidPairs from HiC-Pro)
# Output: SM159.transValidPairs.summary.txt
#         SM159_chrAssnPref_noYM.pdf
#
# To run for SM163, change "SM159" to "SM163" in the filenames below.

Based upon: Moquin SA et al. 2018. The Epstein-Barr Virus Episome 
Maneuvers between Nuclear Chromatin Compartments during Reactivation. 
Journal of Virology 92(3): e01413-17.

Expected = [ 	( (chrA/all) x (chrB/(all-chrA)) ) + 
				( (chrB/all) x (chrA/(all-chrB)) ) ] x totalPairs 
“chrA” represents the number of single-end interchromosomal reads 
		containing chromosome A, 
“chrB” represents the number of single-end interchromosomal reads 
		containing chromosome B, 
“totalPairs” represents the total number of interchromosomal 
		paired-end reads,
“all” represents the total number of interchromosomal single-end 
		reads, which is equal to 2 × totalPairs.
New definitions:
“totalPairs” are all read pairs that are trans interactions
“chrA” are the number of these trans interactions where one of the 
			reads aligns to chromosome A
"""

########################################
# MODULES

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

########################################
# FUNCTIONS

def expected_interactions(chrA, chrB, totalPairs):
    """ calculates expected interaction value for two chromosomes """
    partA = chrA / (2 * totalPairs)
    partB = chrB / ( (2 * totalPairs) - chrA )
    partC = chrB / (2 * totalPairs)
    partD = chrA / ( (2 * totalPairs) - chrB )
    expected  = ( ( partA * partB ) + ( partC * partD ) ) * totalPairs
    return(expected) 
    
def chr_assn_pref(chr1, chr2, countsData, totalPairs):
	chrA = countsData[chr1]
	chrB = countsData[chr2]
	expected = expected_interactions(chrA, chrB, totalPairs)
	if chr1 + "." + chr2 in countsData.keys():
		actual = countsData[ chr1 + "." + chr2 ]
	else:
		actual = countsData[ chr2 + "." + chr1 ]
	return (actual / expected)
	
def get_totalPairs(countsData):
	totalPairs = sum([value for key,value in countsData.items() 
		if "." in key])
	return(totalPairs)

def prep_heatmap(countsData):
	totalPairs = get_totalPairs(countsData)
	#chrs = [key for key in countsData.keys() if "." not in key]
	#chrs.sort()
	chrs = ['HPV31REF','chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
	        'chr9','chr10','chr11','chr12','chr13','chr14','chr15',
	        'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
	df = pd.DataFrame(index=chrs, columns=chrs)
	df = df.fillna(1)
	for i in range(len(chrs)):
		for j in range(i+1,len(chrs)):
			chr1 = chrs[i]
			chr2 = chrs[j]
			score = chr_assn_pref(chr1, chr2, countsData, totalPairs)
			df.loc[chr1,chr2] = score
			df.loc[chr2,chr1] = score
	return(df)

########################################
### Part 1
# Input: Files of only trans valid pairs files
# Output: Two column table where the first column the name of chromosomes
# or interacting chromosomes. The second column is counts:
# counts of SE aligning to that chromosome or of PE aligning to that 
# chromosome pair.

TransValidPairsFile = "SM159.transValidPairs"
outputFile = "SM159.transValidPairs.summary.txt"

ind={}
paired={}

f = open(TransValidPairsFile, 'r')

for i,line in enumerate(f): # use this to read one line at a time
    chrA = line.split('\t')[1]
    chrB = line.split('\t')[4]
    ind.setdefault(chrA,[]).append(1)
    ind.setdefault(chrB,[]).append(1)
    paired.setdefault(chrA + "." +  chrB,[]).append(1)
    if (i % 1000000 ) == 0:
        print(i)

f.close()

ind2 = [ key + "\t" + str(sum(value)) for (key,value) in ind.items() ]
paired2 = [ key + "\t" + str(sum(value)) for (key,value) in paired.items() ]

f = open(outputFile, 'w')
f.write( '\n'.join(ind2) + '\n' + '\n'.join(paired2) )
f.close()

### Part 2
# Input: file from part 1

inputFile = "SM159.transValidPairs.summary.txt"

f = open(inputFile, 'r')
countsData = f.readlines()
f.close()

countsData = [ line.strip().split('\t') for line in countsData ]
countsData = { line[0]:int(line[1]) for line in countsData }

countsData2 = { key:value for (key,value) in countsData.items() 
	if "chrY" not in key }
countsData3 = { key:value for (key,value) in countsData2.items() 
	if "chrM" not in key }

df = prep_heatmap(countsData)

ax = sns.heatmap(df,cmap="seismic",center=1)
plt.savefig("SM159_chrAssnPref_noYM.pdf")
