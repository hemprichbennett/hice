#This script takes all of your input sequences, filters them by length and collapses them by haplotypes, then excludes whatever exclusion threshold you have set

import re
import sys
import string
import gc


del sys.argv[0]

if not len(sys.argv) ==3:
 	sys.exit("Please specify three arguments, 1: the fasta, 2: the groups file, 3: the outfile name")
 
accessionPattern = re.compile(">.+") 
 
inFasta = sys.argv[0]
fastaData = open(inFasta, mode = 'r') 

inGroups = sys.argv[1]
groupsData = open(inGroups, mode = 'r')

outGroups = sys.argv[2]
outData = open(outGroups, mode = 'w')

num_lines = int(sum(1 for line in fastaData))/2

namesList= [None] * int(num_lines)

i = 1

#gc.disable()
#for line in fastaData:

#Make a big list of all of the sequence names that survived our previous QC
with open(inFasta) as input_file:
	for i, line in enumerate(input_file):
		#print(line)
		line = line.rstrip()
		#print(line)
		if line.startswith('>'):
			#print(line)
			line = line.replace('>','')
			#print(line)
			namesList.append(line)

			
			#i +=1

print('namesList filled')


#Check the groups file against the namesList, outputting only sequences which were found in our fasta file
for line in groupsData:
	line = line.rstrip()
	lineList = line.split('\t')
	#print(lineList[0])
	if lineList[0] in namesList:
		#print(lineList[0])
		outData.write(line + '\n')

quit()
