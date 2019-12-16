# Demultiplex

*******************************************************************************

Program	   : demultiplex.py

Description: determine the level of index hopping and barcode errors using 
the demultiplexing process 

Parameters : none

Returned	 : 1. read quality figures and statistics 
2. demultiplexed FASTQ files

******************************************************************************* 

Vocab

BARCODES - the unique sequences that were attached to your each invidivual 
samples’ genetic material before the samples got all mixed together. 
https://astrobiomike.github.io/amplicon/demultiplexing

DEMULTIPLEXING - refers to the step in processing where you’d use the barcode 
information in order to know which sequences came from which samples after 
they had all be sequenced together.
https://astrobiomike.github.io/amplicon/demultiplexing

INDEX SWAPPING - a specific type of misassignment that results in the incorrect 
assignment of libraries from the expected index to a different index 
(in the multiplexed pool). 
https://www.illumina.com/science/education/minimizing-index-hopping.html
