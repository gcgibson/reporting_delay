data = []
with open ("province-biweek_with_delays.csv") as f:
	for line in f.readlines():
		data.append(line.split(','))

count = 0 

for dat in data:
	print (dat)
	if dat[2] == "10":
		count +=1

print (count)
