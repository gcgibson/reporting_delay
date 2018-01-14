import sys
import numpy as np
np.set_printoptions(threshold=np.inf)

n_t_d = []
with open("province-biweek_with_delays.csv") as f:
	i = 0
	for line in f.readlines():
		if i > 0:
			n_t_d.append(line.replace("\n","").split(','))
		i+=1

date_to_index = {}

i = 0
for elm in n_t_d:
	date_to_index[elm[0]+elm[1]] = i

	i+=1


d_to_i = {}
i = 0
iter_ =  date_to_index.keys()
iter_.sort()
for key in iter_:

	d_to_i[key] = i
	i+=1

print (d_to_i)
n_t = np.zeros((i,i))

for elm in n_t_d:
	try:
		sick_date = d_to_i[elm[0]+elm[1]]
		report_date = d_to_i[elm[4] + elm[5]]
		n_t[sick_date][report_date] = elm[3]
	except:
		pass


D= 10

n_t_d = []
for row in range(len(n_t)):
	if len(n_t[row][row:row+D]) == D:
		n_t_d.append(n_t[row][row:row+D].tolist())

np.save("n_t_d",n_t_d)
