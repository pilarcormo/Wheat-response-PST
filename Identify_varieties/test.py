sample = open("LIB12451_Reference_greater_2x.tab", "r")

array=[]
for line in sample:
	x=line.split('\t')
	array.append([x[2].strip('\n'),x[4].strip('\n')])

for i in range(len(array)):
	if array[i][0]!=array[i][1]:
		print array[i]
	else:
		print array[i]
