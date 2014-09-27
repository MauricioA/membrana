# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import os

MESH = 128		#cant de angulos discretos
NLINEAS = 3
MAX_LIN = 10e-6
FOLDERS = (
	'C:\\Archivos\\Membrana\\celula\\salida\\25-128\\40KVm\\',
	'C:\\Archivos\\Membrana\\celula\\salida\\25-128\\80KVm\\',
)

font = { 'size': 16 }

files = ('%sitv.csv' % FOLDERS[0], '%sitv.csv' % FOLDERS[1])
iis = [int(round(((MESH-1) * i) / (NLINEAS-1), 0)) for i in range(NLINEAS)]
angles = [None for i in range(NLINEAS)]
times = ([[] for x in range(NLINEAS)], [[] for x in range(NLINEAS)])
itvs = ([[] for x in range(NLINEAS)], [[] for x in range(NLINEAS)])
first = True

for j in range(2):
	file = files[j]
	with open(file) as f:
		for line in f:
			spl = line.split(',')
			time = float(spl[0])
			angle = float(spl[1])
			itv = float(spl[2])
			if time > MAX_LIN: break

			if first or abs(time - lasttime) > 1e-12:
				i = 0
				first = False
			else:
				i += 1

			lasttime = time

			if i in iis:
				pos = iis.index(i)
				times[j][iis.index(i)].append(time)
				itvs[j][iis.index(i)].append(itv)
				angles[pos] = angle

fig, ax = plt.subplots()
fig.set_size_inches(16, 6)

extra_string = ('40KV/m', '80KV/m')
style = ('-', '--')

for i in range(2):
	for j in range(NLINEAS):
		tms = np.array(times[i][j])
		its = np.array(itvs[i][j])
		ax.plot(tms*1e6, its, label =  extra_string[i] + '  %.0f$^\circ$' % (angles[j] * 57.2957795), ls = style[i])

for label in ax.legend().get_texts():
	label.set_fontsize('small')

plt.ylabel('ITV [V]', fontdict=font)
plt.xlabel('Tiempo [$\\mu$s]', fontdict=font)
plt.xlim(0, MAX_LIN*1e6)

filename = 'C:\\Temp\\temp.png'
pylab.savefig(filename, bbox_inches='tight')
print filename