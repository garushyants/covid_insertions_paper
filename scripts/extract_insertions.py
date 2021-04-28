#!/usr/bin/env python

import sys
from Bio import AlignIO
import optparse
from bisect import bisect_left

#Options
parser=optparse.OptionParser()
parser.add_option('-i', '--infile', default = "/home/garushyants/sars_cov2_nih/GISAID_long_insertions_Feb12/GISAID_long_insertions_aligned_mod.fasta",help='input alignment in fasta format', type='str')
parser.add_option('-r', '--reference', default= 'NC_045512.2', help='reference sequence ID', type='str')
parser.add_option('-o', '--outfile', default= 'out.csv', help='name of the output file with insertions', type='str')

options, args=parser.parse_args()
alignmentFile = options.infile
refseq_name = options.reference
outfile = options.outfile

AnnotationFile = "NC_045512.2.gff3"
allowedLetters = ["a","t","g","c"]

##
#Read annotation
AnnotDict = dict()

with open(AnnotationFile, 'r') as annotf:
	for line in annotf:
		if not line.startswith("#"):
			line =line.rstrip("\n")
			if not line == "":
				parts = line.split("\t")
				if parts[2] == "gene":
					start = int(parts[3])
					stop = int(parts[4])
					name = parts[8].split(";")[2]
					AnnotDict[start] = [stop, name]

##
#Manipulations with alignment
#read_alignment
alignment = AlignIO.read(alignmentFile, "fasta")
alilength = alignment.get_alignment_length()

##
def interval_extract(list): 
    length = len(list) 
    i = 0
    while (i< length): 
        low = list[i] 
        while i <length-1 and list[i]+1 == list[i + 1]: 
            i += 1
        high = list[i] 
        if (high - low >= 1): 
            yield [low, high] 
        elif (high - low == 1): 
            yield [low, ] 
            yield [high, ] 
        else: 
            yield [low, ] 
        i += 1
##

#Get information about reference sequence
i=0
refnum = -1
RefNumDict = dict()
GapPosArray = list()

for record in alignment:
	#print(record.id+"\t"+str(i))
	if refseq_name in record.id:
		refnum = i
		seqstr = str(record.seq).lower()
		#find coordinates of gaps and write correct positions:
		alipos = 0
		refpos = 0
		for letter in enumerate(seqstr):
			if letter[1] in allowedLetters:
				refpos +=1
			elif letter[1] == "-":
				#print letter[1]
				GapPosArray.append(letter[0])
			RefNumDict[letter[0]] = refpos	
	i+=1
print(refnum)

##Get intervals with insertions
print(len(GapPosArray))
GapPosIntervals = list(interval_extract(GapPosArray))	
print(GapPosIntervals)	

##Find insertions
with open(outfile,'w') as ouf:
	for interval in GapPosIntervals:
		#print(interval)
		if len(interval)>1:
			if (RefNumDict[interval[0]] > 99) and (RefNumDict[interval[1]] < 29803): #we remove first and last 100 nucleotides because they are not well aligned
				alignmentslice = alignment[:,interval[0]:interval[1]+1]
				#print(alignmentslice)
				for line in alignmentslice:
					#print(line.seq)
					alistr = str(line.seq).lower()
					alistr_nogap =alistr.replace("-","")
					if (alistr_nogap.isalpha()) and (float(alistr_nogap.count('n'))/float(len(alistr_nogap)) < 0.5):
						annot = ""
						#get annotation for the gene
						RefCoord = RefNumDict[interval[0]]
						annotlist = sorted(list(AnnotDict))
						index = bisect_left(annotlist,RefCoord)
						if RefCoord < AnnotDict[annotlist[index-1]][0]:
							annot = AnnotDict[annotlist[index-1]][1].replace("Name=","")
						#
						print(line.id+"\t"+annot+"\t"+str(RefCoord)+"\t"+str(len(alistr_nogap))+"\t"+alistr_nogap)
						ouf.write(line.id+"\t"+annot+"\t"+str(RefCoord)+"\t"+str(len(alistr_nogap))+"\t"+alistr_nogap+"\n")
