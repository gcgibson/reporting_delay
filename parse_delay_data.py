
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
for key in date_to_index.keys():
	d_to_i[key] = i
	i+=1


n_t = np.zeros((i,i))

for elm in n_t_d:
	try:
		sick_date = d_to_i[elm[0]+elm[1]]
		report_date = d_to_i[elm[4] + elm[5]]
		n_t[sick_date][report_date] = elm[3]
	except:
		pass


np.save("n_t_d",n_t) 


