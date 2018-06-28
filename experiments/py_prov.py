lines = []
with open("provinces.txt") as f:
	for line in f.readlines():
		lines.append(line.replace("\n",""))



import numpy as np

for i in np.arange(0,len(lines),2):
	print (lines[i])
