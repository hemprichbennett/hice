#This script takes all of your input sequences, filters them by length and collapses them by haplotypes, then excludes whatever exclusion threshold you have set

import re
import sys
import string

del sys.argv[0]

if not len(sys.argv) ==6:
 	sys.exit("Please specify six arguments, 1: the input file, 2: output file, 3: output metadata file, 4: minimum sequence length and 5: maximum sequence length, minimum copies")
 
accessionPattern = re.compile(">.+") 
 
inFile = sys.argv[0]
inData = open(inFile, mode = 'r') 

outFile = sys.argv[1]
outData = open(outFile, mode = 'w')

outMetaFile = sys.argv[2]
outMetaData = open(outMetaFile, mode = 'w')
outMetaData.write('sample, sequence_length\n')

minLength = sys.argv[3]
minLength = int(minLength)

maxLength = sys.argv[4]
maxLength = int(maxLength)

minCopies = sys.argv[5]
minCopies = int(minCopies)

n = 0
m = 0
u = 0
s = 0
seqDict ={}

for line in inData:
	line = line.rstrip()
	if accessionPattern.match(line):
		n +=1
		#name = line
		lineList = line.split('_')
		name = lineList[1]
		name = name.replace('[', '.')
		name = name.replace(']', '.')
		name = name.replace('-', '.')
		name = name.split('.')[0]
        
		sequence = next(inData)
		sequence = sequence.replace('\n','')
		seqLength = len(sequence)
		if seqLength > minLength and seqLength < maxLength:
			m +=1
			#outData.write('>' + name + '\n' + sequence)
			#outMetaData.write(name + ',' + str(seqLength) + '\n')
			
			if sequence in seqDict:
				if name in seqDict[sequence]:
					seqDict[sequence][name] +=1
				else:
					seqDict[sequence][name] = 1	
			else:
				seqDict[sequence] = {}
				seqDict[sequence][name] = 1	
		else:
			u +=1
		outMetaData.write(name + ',' + str(seqLength) + '\n')
#print(seqDict)

i = 0 #we need to use a counter to ensure all sequence names that we output are unique, otherwise QIIME gets angry

for sequence in seqDict:
#	seq = key
	i+=1
	for name in seqDict[sequence]:

		value = seqDict[sequence][name]
		if value < minCopies:
			s +=1
			continue
		#print(str(sequence), str(name), (value))
		outstring = '>'+ name + '_' + str(i) + '-' + str(value) + '\n' + sequence + '\n'		
		#print(outstring)
		outData.write(outstring)
print('of ' + str(n) + ' sequences detected, ' + str(m) + ' were of the correct length, '	+ str(u) + ' were excluded for length, and ' + str(s) + ' were excluded due to their number of representatives')
quit()
