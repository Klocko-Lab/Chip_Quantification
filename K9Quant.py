# K9Quant.py
# Jonathan M. Galazka, Andrew D. Klocko
# 2015
# Usage: python ./K9Quant.py
#
# Takes a whole genome observed/expected HiC map and outputs the number of intra- and inter-
# chromosomal interactions above a given threshold.
# Edited by Andrew Reckard for Andrew D. Klocko

import sys
import numpy as np

# Process K9 track into 1 (enriched) or 0 (not enriched)
'''Will need to specify our ChIPk9_path'''
k9_path = 'SRR6058127_WT-N150-H3K9me3_ChIPseq-all_to-nc14_bin10kp.txt'
k9_array = np.transpose(np.loadtxt(k9_path))
k9_median = np.median(k9_array)
k9_median_plus = k9_median + (k9_median * 1)
for i in range(0, k9_array.size):
	if k9_array[i] < k9_median_plus:
		k9_array[i] = 0
	else:
		k9_array[i] = 1

# Load HiC observed/expected dataset
'''Will need to specify our dataset (and resolution?)'''
dataset = 'NMF39'
resolution = 10000

# get chromosome starts and ends
res_string = str(resolution)
'''Specify our chr_starts_path'''
chr_starts_path = 'ChromosomeStarts.txt'
chr_starts_array = np.loadtxt(chr_starts_path, delimiter=' ')
chr_starts = np.transpose(chr_starts_array.astype(int)[0:7])
chr_ends = np.transpose(chr_starts_array.astype(int)[1:8])


'''Specify our dataset_path'''
dataset_path = 'GSM1825702_NMF39_WT_wholegenome_HiCarray_10k_observed_expected.txt'
array = np.loadtxt(dataset_path, delimiter=' ')


lg_array = array[chr_starts[0]:chr_ends[6],chr_starts[0]:chr_ends[6]] # get array of 7 chromosomes
array_x_dim, array_y_dim = lg_array.shape # get the dimensions of this array


# count number of links above threshold
total_inter_links = 0
inter_links_k92k9 = 0
total_intra_links = 0
intra_links_k92k9 = 0
for i in range(0, array_y_dim):

	if(i >= chr_starts[0] and i < chr_ends[0]):
		lg_string1 = 'LGI'
		start1 = (i - chr_starts[0]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[1] and i < chr_ends[1]):
		lg_string1 = 'LGII'
		start1 = (i - chr_starts[1]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[2] and i < chr_ends[2]):
		lg_string1 = 'LGIII'
		start1 = (i - chr_starts[2]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[3] and i < chr_ends[3]):
		lg_string1 = 'LGIV'
		start1 = (i - chr_starts[3]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[4] and i < chr_ends[4]):
		lg_string1 = 'LGV'
		start1 = (i - chr_starts[4]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[5] and i < chr_ends[5]):
		lg_string1 = 'LGVI'
		start1 = (i - chr_starts[5]) * resolution
		end1 = start1 + resolution
		
	elif(i >= chr_starts[6] and i < chr_ends[6]):
		lg_string1 = 'LGVII'
		start1 = (i - chr_starts[6]) * resolution
		end1 = start1 + resolution
	
	for j in range(0, array_x_dim):
	
		if(j >= chr_starts[0] and j < chr_ends[0]):
			lg_string2 = 'LGI'
			start2 = (j - chr_starts[0]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[1] and j < chr_ends[1]):
			lg_string2 = 'LGII'
			start2 = (j - chr_starts[1]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[2] and j < chr_ends[2]):
			lg_string2 = 'LGIII'
			start2 = (j - chr_starts[2]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[3] and j < chr_ends[3]):
			lg_string2 = 'LGIV'
			start2 = (j - chr_starts[3]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[4] and j < chr_ends[4]):
			lg_string2 = 'LGV'
			start2 = (j - chr_starts[4]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[5] and j < chr_ends[5]):
			lg_string2 = 'LGVI'
			start2 = (j - chr_starts[5]) * resolution
			end2 = start2 + resolution
			
		elif(j >= chr_starts[6] and j < chr_ends[6]):
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
				k9_array[i] == 0):
		
			pass
		
		# determine if region passes thresholds, Thresholds are altered from histone to histone?
		else:

			value = np.log2(value)

			if(value >= 3.5 and lg_string1 != lg_string2):
	
				total_inter_links = total_inter_links + 1
				
			if(value >= 3.5 and k9_array[j] == 1 and lg_string1 != lg_string2):
			
				inter_links_k92k9 = inter_links_k92k9 + 1

			if(value >= 2.25 and lg_string1 == lg_string2):
	
				total_intra_links = total_intra_links + 1
				
			if(value >= 2.25 and k9_array[j] == 1 and lg_string1 == lg_string2):
			
				intra_links_k92k9 = intra_links_k92k9 + 1
				
# output				
print('\n')
print(dataset + ' link quantification:')				
print('Inter links: ' + str(total_inter_links))
print('Inter links K9 to K9: ' + str(inter_links_k92k9))
print('Inter links K9 to not K9: ' + str(total_inter_links - inter_links_k92k9))
print('\n')
print('Intra links: ' + str(total_intra_links))
print('Intra links K9 to K9: ' + str(intra_links_k92k9))
print('Intra links K9 to not K9: ' + str(total_intra_links - intra_links_k92k9))

					
					
					
					
					
					
					
					
					
					