# K36-K9Quant.py
# Jonathan M. Galazka, Andrew D. Klocko
# 2015
# Usage: python ./K9Quant.py
#
# Takes a whole genome observed/expected HiC map and outputs the number of intra- and inter-
# chromosomal interactions above a given threshold.
# Edited by Andrew Reckard for Andrew D. Klocko
# Modified by OluwadareLab to Accept Input Argument for .h5 file converted to a N x N Square Matrix

import sys
import numpy as np
import h5py
import argparse

## epi1 = first epigenetic mark to assess for enrichment
## epi2 = second epigenetic marks to assess for enrichment


# Process first-epi-mark track into 1 (enriched) or 0 (not enriched)
'''Will need to specify our first-epi-mark_path'''
epi1_path = 'SRR6058127_WT-N150-H3K9me3_ChIPseq-nc14_10kbins_norm-RPKM-correct.txt'
epi1_array = np.transpose(np.loadtxt(epi1_path))
epi1_median = np.median(epi1_array)
epi1_median_plus = epi1_median + (epi1_median * 1)
for i in range(0, epi1_array.size):
	if epi1_array[i] < epi1_median_plus:
		epi1_array[i] = 0
	else:
		epi1_array[i] = 1

# Process second-epi-mark track into 1 (enriched) or 0 (not enriched)
'''Will need to specify our second-epi-mark_path'''
epi2_path = 'SRR6058127_WT-N150-H3K9me3_ChIPseq-nc14_10kbins_norm-RPKM-correct.txt'
epi2_array = np.transpose(np.loadtxt(epi2_path))
epi2_median = np.median(epi2_array)
epi2_median_plus = epi2_median + (epi2_median * 1)
for i in range(0, epi2_array.size):
	if epi2_array[i] < epi2_median_plus:
		epi2_array[i] = 0
	else:
		epi2_array[i] = 1

# Load HiC observed/expected dataset
'''Will need to specify our dataset (and resolution?)'''
#dataset = 'N150'
resolution = 10000

# get chromosome starts and ends
res_string = str(resolution)
'''Specify our chr_starts_path'''
chr_starts_path = 'ChromosomeStartPoints_nc14-10kb.txt'
chr_ends_path = 'ChromosomeEndPoints_nc14-10kb.txt'
chr_starts_array = np.loadtxt(chr_starts_path, delimiter=' ')
chr_ends_array = np.loadtxt(chr_ends_path, delimiter=' ')
chr_starts = np.transpose(chr_starts_array.astype(int)[0:7])
chr_ends = np.transpose(chr_ends_array.astype(int)[0:7])


'''Specify our dataset_path'''
# dataset_path = 'GSM1825702_NMF39_WT_wholegenome_HiCarray_10k_observed_expected.txt'
# array = np.loadtxt(dataset_path, delimiter=' ')


'''This is the new code, the test.txt is the N X N matrix generated from the homer file generated with hicConvert. 
	This N X N matrix generation process happens in the dataconvert.py file.
'''
parser=argparse.ArgumentParser()
parser.add_argument("contact_matrix",help="File name of a white space delimited square contact matrix")
args=parser.parse_args()
if args.contact_matrix:
	dataset_path = args.contact_matrix  # this is the n x n matrix output in a txt file
	
array = np.genfromtxt(dataset_path, delimiter=' ') #[:,:-1]


lg_array = array[chr_starts[0]:chr_ends[6],chr_starts[0]:chr_ends[6]] # get array of 7 chromosomes
array_x_dim, array_y_dim = lg_array.shape # get the dimensions of this array

print('The size of your array is:')
print(array_x_dim, array_y_dim)


# count number of links above threshold
total_inter_links = 0
inter_links_first2second = 0
total_intra_links = 0
intra_links_first2second = 0
for i in range(0, array_y_dim):

	if(i >= chr_starts[0] and i <= chr_ends[0]):
		lg_string1 = 'LGI'
		start1 = (i - chr_starts[0]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[1] and i <= chr_ends[1]):
		lg_string1 = 'LGII'
		start1 = (i - chr_starts[1]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[2] and i <= chr_ends[2]):
		lg_string1 = 'LGIII'
		start1 = (i - chr_starts[2]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[3] and i <= chr_ends[3]):
		lg_string1 = 'LGIV'
		start1 = (i - chr_starts[3]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[4] and i <= chr_ends[4]):
		lg_string1 = 'LGV'
		start1 = (i - chr_starts[4]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[5] and i <= chr_ends[5]):
		lg_string1 = 'LGVI'
		start1 = (i - chr_starts[5]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[6] and i <= chr_ends[6]):
		lg_string1 = 'LGVII'
		start1 = (i - chr_starts[6]) * resolution
		end1 = start1 + resolution
	
	for j in range(0, array_x_dim):
	
		if(j >= chr_starts[0] and j <= chr_ends[0]):
			lg_string2 = 'LGI'
			start2 = (j - chr_starts[0]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[1] and j <= chr_ends[1]):
			lg_string2 = 'LGII'
			start2 = (j - chr_starts[1]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[2] and j <= chr_ends[2]):
			lg_string2 = 'LGIII'
			start2 = (j - chr_starts[2]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[3] and j <= chr_ends[3]):
			lg_string2 = 'LGIV'
			start2 = (j - chr_starts[3]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[4] and j <= chr_ends[4]):
			lg_string2 = 'LGV'
			start2 = (j - chr_starts[4]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[5] and j <= chr_ends[5]):
			lg_string2 = 'LGVI'
			start2 = (j - chr_starts[5]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[6] and j <= chr_ends[6]):
			lg_string2 = 'LGVII'
			start2 = (j - chr_starts[6]) * resolution
			end2 = start2 + resolution

		value = float(lg_array[i,j])
		
		# filter out nan, -inf, inf, zeros or non-K9-enriched regions
		if(
				np.isnan(value) == True or
				np.isposinf(value) == True or
				np.isneginf(value) or
				value == 0.00 or
				epi1_array[i] == 0):
		
			pass
		
		# determine if region passes thresholds, Thresholds are altered from histone to histone?
		else:
			###uncomment this out if we need to log transform the array - but obs/exp may already be log transformed
			
			value = np.log2(value)

			# for calculating inter links
			## change the value based on the signal strength of the matrix - before, used 3.5 as the value because we had an obs-v-exp matrix, and that displayed the log2 difference
			## will try 1000; but we are taking the log2 of it, so perhaps 3.5 is fine
			

			if(value >= 2.25 and lg_string1 != lg_string2):
				total_inter_links = total_inter_links + 1
			
			if(value >= 2.25 and epi2_array[j] == 1 and lg_string1 != lg_string2):
			
				inter_links_first2second = inter_links_first2second + 1
				
			##for calculating intra links
			## change the value based on the signal strength of the matrix - before, used 2.25 as the value because we had an obs-v-exp matrix, and that displayed the log2 difference
			##will try 1000; but we are taking the log2 of it, so perhaps 3.5 is fine
			
			if(value >= 1.75 and lg_string1 == lg_string2):
			
				total_intra_links = total_intra_links + 1
				
			if(value >= 1.75 and epi2_array[j] == 1 and lg_string1 == lg_string2):
			
				intra_links_first2second = intra_links_first2second + 1
				
# output				
print('\n')
print('Here are the link quantification results of your strain:')
print('Inter links: ' + str(total_inter_links))
print('Inter links first-mark to second-mark: ' + str(inter_links_first2second))
print('Inter links first-mark to not second-mark: ' + str(total_inter_links - inter_links_first2second))
print('\n')
print('Intra links: ' + str(total_intra_links))
print('Intra links first-mark to second-mark: ' + str(intra_links_first2second))
print('Intra links first-mark to not second-mark: ' + str(total_intra_links - intra_links_first2second))

					
					
					
					
					
					
					
					
					
					
