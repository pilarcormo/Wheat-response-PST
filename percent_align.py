f = open('2015_isolates/LIB18751/top_hat/align_summary.txt')

for line in f:
	if "Mapped" in line:
  		a = line.split(' ')[2]
  		print a