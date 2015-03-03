# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import os
import re

# para generar graficos superpuestos de cant de poros vs tiempo de carpeta de poros

FOLDERS = [
	('../salida/pulsativo/02-16/', '1600 V/cm'),
	('../salida/pulsativo/02-17/120kvm/', '1200 V/cm'),
]

font = { 'size': 16 }

def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

plt.clf()
fig, ax = plt.subplots()
fig.set_size_inches(12, 4.5)

for folder, label in FOLDERS:
	poros_folder = folder + '/poros/'
	files = filter(lambda x : x.endswith('.csv') and x.startswith('poros'), os.listdir(poros_folder))
	files = map(lambda x : poros_folder + x, files)

	values = []
	
	for file in files:
		reg = re.match('.*\\poros\-\d+\-(.*)\.csv', file)
		time = float(reg.groups()[0]) * 1000
		poros = count_lines(file)
		values.append((time, poros))

	values = sorted(values, key = lambda x : x[0])
	
	times = [x[0] for x in values]
	poros = [x[1] for x in values]

	times = np.array(times)
	poros = np.array(poros)

	ax.plot(times, poros, label = label)
	
ax.legend(loc = 'upper left')

plt.xlabel('Tiempo [ms]')
plt.ylabel('Poros')

#plt.show()

outfile = 'temp.png'
pylab.savefig(outfile, bbox_inches='tight')
print outfile
