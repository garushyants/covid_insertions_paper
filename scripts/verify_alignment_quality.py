#!/usr/bin/env python3

import sys
from Bio import AlignIO
import optparse
from bisect import bisect_left

#Options
parser=optparse.OptionParser()
parser.add_option('-i', '--inali', default = "./msa0617_genomes_with_inserts_oneref.fasta",help='input alignment in fasta format', type='str')
parser.add_option('-r', '--reference', default= 'hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31|Asia', help='reference sequence ID', type='str')
parser.add_option('-n', '--listofinserts', default= '../data/GISAID_msa0617_insertions_filtered.csv', help='reference sequence ID', type='str')
parser.add_option('-o', '--outpath', default= '../data/Short_alignments/', help='name of the output file with insertions', type='str')

options, args=parser.parse_args()
alignmentFile = options.inali
listofinserts = options.listofinserts
refseq_name = options.reference
outpath = options.outpath

################
##Read inserts description

insertsDict =dict()
insertsMetadataDict = dict()
GeneDict = dict()

with open(listofinserts,'r') as insf:
	next(insf)
	for ll in insf:
		genome,gene,position,length,seq = ll.split("\t")
		#
		GeneDict[position] = gene
		#
		pnum =int(position)
		if not pnum in insertsDict:
			insertsDict[pnum] = list()
		insertsDict[pnum].append(int(length))
		#
		if not position+"_"+seq in insertsMetadataDict:
			insertsMetadataDict[position+"_"+seq] =list()
		insertsMetadataDict[position+"_"+seq].append(genome)	
		

##
#Manipulations with alignment
#read_alignment
alignment = AlignIO.read(alignmentFile, "fasta")
alilength = alignment.get_alignment_length()

#Get information about reference sequence
i=0
refnum = -1
RefNumDict = dict()
GapPosArray = list()

allowedLetters = ["a","t","g","c"]

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
			RefNumDict[int(refpos)] = int(letter[0])
	i+=1
print(refnum)
#print(RefNumDict)

#Get short alignment
border = 30
for p in insertsDict:
	max_length = max(insertsDict[p])
	print(max_length)
	start = RefNumDict[p-1] - border
	end = RefNumDict[p] + max_length + border
	alignmentslice = alignment[:,start:end]
	print(alignmentslice)
	output_handle = open(outpath+str(p)+".fasta",'w')
	AlignIO.write(alignmentslice, output_handle, "fasta")
