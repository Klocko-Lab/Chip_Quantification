#!/bin/bash
# Script created by Oluwadarelab for Klocko's Lab Project, 2021
#
# Please make sure the following Requirements are installed:
#	HiCExplorer library: You should have HiCexplorer library installed because we used the hicConvertFormat function
#	gzip installed: We used gunzip to uncompress too 
#Execute script = ./process.sh
#Outputs Generated: 
#		HiC-Explorer Output compressed (.gz) file
#     	A N by N Matrix


#SPECIFY INPUT AND OUTPUT FILE NAMES
h5input=$1
hiCExploreroutput="hiCExploreroutput"   #Optional
NxN_Matrix="NxN_Matrix.txt"				#Optional

# Genrate homer output using hicConvertFormat as compressed file
echo "h5 file input filename: $1"
hicConvertFormat --matrices $h5input  --inputFormat h5 -o "${hiCExploreroutput}.gz"  --outputFormat homer

#Uncompress output file
gunzip -f "${hiCExploreroutput}.gz"

#Convert HicExplorer .homer Output to  N by N Square Matrix
python dataconvert.py $hiCExploreroutput $NxN_Matrix

#Run Script from Klocko Lab with Square Matrix as Input (We changed some portions of your code to accept input argument )
python epigenetic-mark-Quant.py $NxN_Matrix

echo "Processing Completed Successfully"
