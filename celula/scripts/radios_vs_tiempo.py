# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import os
import re

# para generar graficos superpuestos de radios vs tiempo de carpeta de poros

FOLDERS = [
	('../salida/pulsativo/02-16/', '1600 V/cm'),
	('../salida/pulsativo/02-17/120kvm/', '1200 V/cm'),
]

font = { 'size': 16 }

def get_radio_avg(fname):
	total = 0
	with open(fname) as f:
		for i, line in enumerate(f):
			if i == 0: continue
			total += float(line.split(',')[1])
	avg = 0
	if i > 0: avg = total / 1
	#print 1, avg
	return avg

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
		radios = get_radio_avg(file) * 1e9
		values.append((time, radios))

	values = sorted(values, key = lambda x : x[0])
	
	times = [x[0] for x in values]
	radios = [x[1] for x in values]

	times = np.array(times)
	radios = np.array(radios)

	ax.plot(times, radios, label = label)
	#plt.semilogx(times, radios, label = label)
	
ax.legend(loc = 'upper right')

#ax.set_xlim(0, 5)

plt.xlabel('Tiempo [ms]')
plt.ylabel('Radio [nm]')

#plt.show()

outfile = 'temp.png'
pylab.savefig(outfile, bbox_inches='tight')
print outfile
