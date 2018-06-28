import numpy as np
import matplotlib.pyplot as plt


data = [[ 347.28100104262353,153.34373192527656],
        [ 11228.117369561476,10946.19174444719],
        [ 5717.183655722671,5666.965509178719],
        [ 4208.343619339098,4156.521913956512],
        [1312.0241200924213,1204.4695259041287],
	[1568.943615868976,1516.0333648674161],
	[1125.4026182610548,1089.4167834789414],
	[632.2973126044274,747.9118780492763],
	[605.5663508720593,641.4845399247652],
	[31442.662605386042,30782.027513819696]]

columns = ('Ignore Delay', 'Delay Adjusted')
rows = ['%d' % x for x in (2006,2007,2008,2009,2010,2011,2012,2013,2014,2015)]

fig, ax = plt.subplots()

# Hide axes
ax.xaxis.set_visible(False) 
ax.yaxis.set_visible(False)

# Table from Ed Smith answer
clust_data = np.random.random((10,3))
ax.table(cellText=data,rowLabels=rows,colLabels=columns,loc='center')
plt.show()
