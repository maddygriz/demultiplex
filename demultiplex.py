#!/usr/bin/env python

#*******************************************************************************
#Program	: demultiplex.py
#Description: determine the level of index hopping and barcode errors using 
#               the demultiplexing process 
#Parameters	: none
#Returned	: read quality figures and statistics
#             demultiplexed FASTQ files
#******************************************************************************* 

#Vocab

#BARCODES - the unique sequences that were attached to your each invidivual 
#samples’ genetic material before the samples got all mixed together. 
#https://astrobiomike.github.io/amplicon/demultiplexing

#DEMULTIPLEXING - refers to the step in processing where you’d use the barcode 
#information in order to know which sequences came from which samples after 
#they had all be sequenced together.
#https://astrobiomike.github.io/amplicon/demultiplexing

#INDEX SWAPPING - a specific type of misassignment that results in the incorrect 
#assignment of libraries from the expected index to a different index 
#(in the multiplexed pool). 
#https://www.illumina.com/science/education/minimizing-index-hopping.html


#package install
import argparse
import gzip

#*******************************************************************************
#Function	: get_arguments
#Description: get the arguments from the command line for generalization 
#				of the program
#Parameters	: none
#Returned	: parse_args - the arguments needed for the program
#******************************************************************************* 

def get_arguments():
	parser = argparse.ArgumentParser (description = "graphing fastq scores")
	parser.add_argument("-r", "--read1", help="file for first read", required=True, type=str)
	parser.add_argument("-R", "--read2", help="file for second read", required=True, type=str)
	parser.add_argument("-i", "--index1", help="file for first index", required=True, type=str)
	parser.add_argument("-I", "--index2", help="file for second index", required=True, type=str)
	parser.add_argument("-o", "--out", help="output file for ending summary", required=True, type=str)
	return parser.parse_args()
args = get_arguments()

#*******************************************************************************
#Function	: matching
#Description: to compare two indexes and determine if they are the same
#				(index2 reverse compliment of index1)
#Parameters	: index1  - index from first read
#			  index2  - index from second read
#Returned	: matched - bool
#			    TRUE  - if matched
#				FALSE - if not matched
#******************************************************************************* 

def matching(index1, index2):
    matched = False
    inverse = ""
    
    for x in index2:
        if x == "A":
            inverse =  "T" + inverse 
        elif x == "T":
            inverse =  "A" + inverse
        elif x == "G":
            inverse = "C" + inverse
        elif x == "C":
            inverse = "G" + inverse   
        else:
            inverse = "N" + inverse

#assuming N base calls are normally on the first base
#first 100 reads in both files only have Ns on first base
#this allows one base pair difference between matching indexes
#keeps down the rate of false positives 

    if containsN(index1):
        index1 = index1[1:]
        inverse = inverse[1:]
        
    elif containsN(index2):
        index1 = index1[:-1]
        inverse = inverse[:-1]
        
    elif containsN(index1) and containsN(index2):
        index1 = index1[1:-1]
        inverse = inverse[1:-1]
        
    return (inverse == index1)

	
#*******************************************************************************
#Function	: writeToBadFile
#Description: write reads that don't meet standard to file
#			  Reason concatenated to the + index line in file
#Parameters	: outputFile - fastq file that the reads are being written to
#			  read       - the read needing to be put in the output file
#			  reason     - reason for read failing, concat onto + index line in []
#				options:
#					QS - low quality score
#					N  - index contains Ns
#					IH - index hopping
#Returned	: none
#*******************************************************************************   	
	
def writeToBadFile(outputFile, read, reason):
    N = 1
    IH = 2
    QS = 4
    E = 8 #error in index
    
    if (reason & N):
        read[2] = read[2] + " " + "(N)"
        
    if (reason & IH):
        read[2] = read[2] + " " + "(IH)"
        
    if (reason & QS):
        read[2] = read[2] + " " + "(QS)"
        
    if (reason & E):
        read[2] = read[2] + " " + "(E)"
        
    with open (outputFile, "a") as out:
        out.write(read[0])
        out.write("\n")
        out.write(read[1])
        out.write("\n")
        out.write(read[2])
        out.write("\n")
        out.write(read[3])
        out.write("\n")	
	
#*******************************************************************************
#Function	: writeToGoodFile
#Description: write reads that meet standards to file
#Parameters	: outputFile - fastq file that the reads are being written to
#			  read       - the read needing to be put in the output file
#Returned	: none
#*******************************************************************************
	
def writeToGoodFile(outputFile, read):
    with open (outputFile, "a") as out:
        out.write(read[0])
        out.write("\n")
        out.write(read[1])
        out.write("\n")
        out.write(read[2])
        out.write("\n")
        out.write(read[3])
        out.write("\n")	

#*******************************************************************************
#Function	: convertPhred
#Description: convert letter to a phred score 
#Parameters	: letter - quality score from fastq file
#Returned	: score  - phred score value (int)
#*******************************************************************************

def convertPhred(letter):
    return ord(letter)-33

#*******************************************************************************
#Function	: meanQS
#Description: determine the average quality score in the read
#Parameters	: read - the read that contains the quality scores
#Returned	: mean - mean of quality scores (float)
#*******************************************************************************	

def meanQS(read):
    total = 0
    bp = 0
    for base in read[3]:
        total = total + convertPhred(base)
        bp = bp + 1
    return total/bp
	
#*******************************************************************************
#Function	: compareMeans 
#Description: compare mean quality scores to designated cutoff
#Parameters	: mean - mean of quality scores of specific read
#Returned	: above   - bool
#				TRUE  - if above or equal to cutoff
#				FALSE - if below cutoff
#*******************************************************************************	  	
	
def compareMeans(read):
    cutoffQS = 32
    return meanQS(read) >= cutoffQS

#*******************************************************************************
#Function	: containsN
#Description: determine if the index contains an N
#Parameters	: index - the index of the read
#Returned	: N 	  - bool
#				TRUE  - if index contains NameError
#				FALSE - if index does not contain N
#*******************************************************************************	  

def containsN(index):
    return "N" in index
	
#*******************************************************************************
#Function	: readInLine
#Description: read in one line of the file and strip the white space
#Parameters	: file - the file to be read
#Returned	: line of file with no white space
#*******************************************************************************	 
 
def readInLine(file):
    line = file.readline()
    return line.strip() 
	

#initialize constants
fastqLines = 4
end = False
N = 0
IH = 0
QS = 0
reads = 0


#open files
R1 = gzip.open(args.read1, "r")
R2 = gzip.open(args.read2, "r")
I1 = gzip.open(args.index1, "r")
I2 = gzip.open(args.index2, "r")

#initialize empty arrays
read1   = [None] * 4
indexR1 = [None] * 4
indexR2 = [None] * 4
read2   = [None] * 4

#initialize dict with zero count for every index
indices = {
    'GTAGCGTA': 0, 'CGATCGAT': 0, 'GATCAAGG': 0, 'AACAGCGA': 0,
    'TAGCCATG': 0, 'CGGTAATC': 0, 'CTCTGGAT': 0, 'TACCGGAT': 0, 
    'CTAGCTCA': 0, 'CACTTCAC': 0, 'GCTACTCT': 0, 'ACGATCAG': 0, 
    'TATGGCAC': 0, 'TGTTCCGT': 0, 'GTCCTAAG': 0, 'TCGACAAG': 0,
    'TCTTCGAC': 0, 'ATCATGCG': 0, 'ATCGTGGT': 0, 'TCGAGAGT': 0,
    'TCGGATTC': 0, 'GATCTTGC': 0, 'AGAGTCCA': 0, 'AGGATAGC': 0,
    'SNP': 0
}

#run until end of file
while end == False:
    #initialize bit flag
    flag = 0
    
    #add corresponding index sequence to the "+" line of the read fastq file
    for x in range(fastqLines):
        read1[x] = readInLine(R1)
        if read1[x] == "":
            end = True
            break
        indexR1[x] = readInLine(I1)
        indexR2[x] = readInLine(I2)
        read2[x] = readInLine(R2)
        
        if x == 2:
            read1[2] = read1[2] + " " + indexR1[1]
            index1 = indexR1[1]
            read2[2] = read2[2] + " " + indexR2[1]
            index2 = indexR2[1]
    if end == True:
        break
    
    #increment read counter
    reads = reads + 1 
    
    #add errors to bit flag 
    if (containsN(index1)) or (containsN(index2)):
        flag = flag + 1
        N = N + 1
    
    if not(matching(index1, index2)):
        flag = flag + 2
        IH = IH + 1
    
    if not(compareMeans(index1)) or not(compareMeans(index2)):
        flag = flag + 4
        QS = QS + 1
    
    #if no error and index is valid, write file to "valid" location
    #else write file to "rejects" location
    if flag == 0:
        if index1 in indices:
            indices[index1] = indices[index1] + 1
            writeToGoodFile('/projects/bgmp/maddyg/Bi624/demultiplexing/outputFiles/%s_R1.fastq' % index1, read1)
            writeToGoodFile('/projects/bgmp/maddyg/Bi624/demultiplexing/outputFiles/%s_R2.fastq' % index1, read2)
        else:
            flag = flag + 8
            indices["SNP"] = indices["SNP"] + 1
            writeToBadFile('/projects/bgmp/maddyg/Bi624/demultiplexing/outputFiles/rejects_R1.fastq', read1, flag)
            writeToBadFile('/projects/bgmp/maddyg/Bi624/demultiplexing/outputFiles/rejects_R2.fastq', read1, flag)
    else:
        writeToBadFile('/projects/bgmp/maddyg/Bi624/demultiplexing/outputFiles/rejects_R1.fastq', read1, flag)
        writeToBadFile('/projects/bgmp/maddyg/Bi624/demultiplexing/outputFiles/rejects_R2.fastq', read2, flag)

#write output file
with open(args.out, "w") as out:
    out.write("Reads: ")
    out.write(str(reads))
    out.write("\n")
    out.write("N's  : ")
    out.write(str(N))
    out.write("\n")
    out.write("IH   : ")
    out.write(str(IH))
    out.write("\n")
    out.write("QS   : ")
    out.write(str(QS))
    out.write("\n")
    out.write("\n")
    out.write("% of reads for each index")
    out.write("\n")
    out.write(str(indices))

#close files
R1.close()
R2.close()
I1.close()
I2.close()
