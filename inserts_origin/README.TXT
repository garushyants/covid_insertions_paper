Search for sources of insertions
--------------------------------

Input file with insertion = ins_ext.seq
format:
>name
sequence, one line, insertion is surrounded by spaces

example:

>MN908947.3
GCTAGTTATCAGACTCAGACTAATT CTCCTCGGCGGG CACGTAGTGTAGCTAGTCAATCC



Running the program

olig2a database_name number_of_mismatches strand O/R

number_of_mismatches = maximum allowed number of mismatches
strand = D (the direct strand only), C (direct and complimetary strands)
O/R = O (analysis of the insertion), R (analysis of a shuffled insertion)
example:

olig2a vir.fas 2 D O

olig2a vir.fas 2 D R

vir.fas - database file

OUTput file:
out.txt




Batch file to run the system

RUNP database_name number_of_mismatches strand
number_of_mismatches = maximum allowed number of mismatches
strand = D (the direct strand only), C (direct and complimetary strands)

An example of a batch file:
RUNP.bat vir 2 D



OUTput files of RUNP.bat (example for the vir.fas database):

**** outvir.txt  - output of the best candidate sources of the insertion
example:
>  21479 2 d >Bat_SARS_CoV_HKU31 0 29592 + 29592                         
CTCCTCGGCGGG  - insertion
CTTCTCGTCGGG  - putative source sequence)
(Position of the best source sequence CTTCTCGTCGGG, number of mismatches,
description of a SARS-COV genome)


**** weightvir.txt - statistics of the actual inserion (O) and shuffled version (R)
P value at the end of the file
example:
O   2
R   2
.......
R   2
R   2

P =  0.80

****additional output file:
fnavir.txt 
Extended context of insertions (-+ 150 bp). This file can be used as input for 
FASTA3 or BLASTN program

