varieties = open("SNP-markers/varieties_numbercode.csv")
example_varieties = open("SNP-markers/example/info_7_varieties_nucleotides.txt")
example_code = open("SNP-markers/example/info_7_varieties.txt")


marker = {}
code_dic = {}
code_nuc = {}
decode = {}


for line in varieties:
	codes = []
	a = line.split(',')
	lent= len(a)
	marker[str(a[1])] = a[2:lent]

for line in example_code:
	codes = []
	a = line.strip('\r\n').split('\t')
	lent= len(a)
	code_dic[str(a[0])] = a[5:lent]
	
for line in example_varieties:
	codes = []
	a = line.strip('\r\n').split('\t')
	lent= len(a)
	code_nuc[str(a[0])] = a[1:lent]

n = 0

for k in code_dic.keys():
	if k in code_nuc.keys():
		for array in code_dic[k]:
			print array 
			# while n <= len(array):
			# 	print array 
				# decode[array[n]] = code_nuc[k][n]
				# n = n + 1	
				# print decode 

		# decode[code_dic[k]] = code_nuc[k]
		# print decode 

		# decode[str(k)] = code_dic.values()













