#This script takes your groupfile and searches for duplicates, making sure theres only one of each entry

import re
import sys
import string
import gc


print('Script started')

del sys.argv[0]

if not len(sys.argv) ==2:
 	sys.exit("You have specified " + str(len(sys.argv))+ "arguments. Please specify two arguments, 1: the groups file input, 2: the groups file output")
 
accessionPattern = re.compile(">.+") 
 
inGroups = sys.argv[0]
inData = open(inGroups, mode = 'r') 

outGroups = sys.argv[1]
outData = open(outGroups, mode = 'w')

writtenList = []

#Check the groups file against the namesList, outputting only sequences which were found in our fasta file
for line in inData:
	line = line.rstrip()
	
	#print(lineList[0])
	if line not in writtenList:
		outData.write(line + '\n')
		writtenList.append(line)

quit()