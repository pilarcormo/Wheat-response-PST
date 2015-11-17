from __future__ import division
import sys
import csv

directory  = sys.argv[1] #directory where 
libraries = sys.argv[2]   #name of text file with all the library names listed 
outputfile = sys.argv[3] #name for the output csv file 


print 'The directory where the libraries are is ', directory
print 'The text file with the libary names listed is ', libraries,'.txt'
print 'Output file is ', outputfile,'.csv'


array = []
dic = {}

file_tail=open(''+str(directory)+'/'+str(libraries)+'.txt')


for library in file_tail:
	a =  library.strip('\n')
	name='2015_isolates/'+str(a)+'/top_hat/align_summary.txt'
	f = open(name) 
	array = []
	for line in f:
		if "Mapped" in line:
	  		a = float(line.split(':')[1].split('(')[1].split('%')[0])
	  		array.append(a)
	average = sum(array)/len(array)
	dic[str(library.strip('\n'))] = average

with open(''+str(directory)+'/'+str(outputfile)+'.csv', 'wb') as f:  # Just use 'w' mode in 3.x
    writer = csv.writer(f)
    for row in dic.iteritems():
        writer.writerow(row)


print 'Your out file with the percent aligned to the reference genomes can be found at '+str(directory)+'/'+str(outputfile)+'.csv'