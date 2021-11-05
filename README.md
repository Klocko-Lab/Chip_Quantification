# Chip_Quantification
Programs for quantifying ChIP-seq data

The program process-AK.sh is the shell file that utilizes the dataconvert.py and epigenetic-mark-Quant-v2.py to 
quantify the number of contacts between regions marked by various epigenetic marks. The shell requires an input 
.h5 matrix to run and compare regions of enriched ChIP-seq data to contact regions outlined by the HiC. The ChIP-seq
files are hard coded into the epigenetic-mark-Quant-v2.py program, and need to be edited directly when looking at a
new ChIP-seq data set.

The shell command needs to be converted to executable before running.

Generic Example Shell command
./process.sh name_of_h5_file.h5

This will output a .gz compressed, homer format file and an NxN matrix text file.

This process utilizes the following python libraries to run.

-Hicexplorer library
-gunzip
