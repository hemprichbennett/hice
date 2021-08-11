

import re
import sys
import string
import gc


print('Script started')

del sys.argv[0]

if not len(sys.argv) ==2:
 	sys.exit("You have specified " + str(len(sys.argv))+ "arguments. Please specify two arguments, 1: the fasta, 2: the groups file")
 
accessionPattern = re.compile(">.+") 
 
inFasta = sys.argv[0]
fastaData = open(inFasta, mode = 'r') 

print('fasta loaded')
inGroups = sys.argv[1]
groupsData = open(inGroups, mode = 'r')

print('groups opened')
outGroups = str(inGroups) + '.' + str(inFasta) + '.groups'
outData = open(outGroups, mode = 'w')

print('out opened')

num_lines = int(sum(1 for line in fastaData))/2
print('num_lines known')
namesList= [None] * int(num_lines)

print('namesList made')

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

writtenList = []

#Check the groups file against the namesList, outputting only sequences which were found in our fasta file
for line in groupsData:
	line = line.rstrip()
	lineList = line.split('\t')
	#print(lineList[0])
	if lineList[0] in namesList:
		#print(lineList[0])
		outData.write(line + '\n')
		writtenList.append(line)

quit()
