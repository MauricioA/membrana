# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import os

# para generar graficos de itv vs tiempo de itv.csv

MESH = 192		    # cant de angulos discretos
NLINEAS = 7
MAX_LIN_US = 5   # 10
MAX_LOG_SEG = 1e-9  
MAXTIME = max(MAX_LIN_US, MAX_LOG_SEG)
FOLDERS = [
	'../salida/last/02-03/120/',
]

font = { 'size': 16 }

for folder in FOLDERS:
	file = '%sitv.csv' % folder
		
	times = []
	itvss = [[] for x in range(NLINEAS)]
	iis = [int(round(((MESH-1) * i) / (NLINEAS-1), 0)) for i in range(NLINEAS)]
	angles = [0 for x in range(NLINEAS)]
	first = True
	
	with open(file) as f:
		for line in f:
			spl = line.split(',')
			time = float(spl[0])
			angle = float(spl[1])
			itv = float(spl[2])
			if time > MAXTIME: break
				
			if first or abs(time - lasttime) > 1e-12:
				i = 0
				first = False
			else:
				i += 1

			lasttime = time

			if i in iis:
				pos = iis.index(i)
				if i==0: times.append(time)
				itvss[pos].append(itv)
				angles[pos] = angle

		times = np.array(times)

		plt.clf()
		i = 0
		fig, ax = plt.subplots()
		for itvs in itvss:
				
			if not len(itvs) == len(times):
				print '%s %s' % (len(itvs), len(times))
				assert False
				
			itvs = np.array(itvs)
			ax.plot(times*1e6, itvs, label = '%.0f$^\circ$' % (angles[i] * 57.2957795))
			i += 1

		legend = ax.legend()
		for label in legend.get_texts():
			label.set_fontsize('small')

		plt.ylabel('ITV [V]', fontdict=font)
		plt.xlabel('Tiempo [$\\mu$s]', fontdict=font)
		plt.xlim(0, MAX_LIN_US)
			
		filename = folder + 'itv-lin.png'
		pylab.savefig(filename,	bbox_inches='tight')
		print filename

		#logarítmico:

		#plt.clf()
		#i = 0
		#fig, ax = plt.subplots()
		#for itvs in itvss:
		#	assert len(itvs) == len(times)
		#	itvs = np.array(itvs)
		#	plt.semilogx(times, itvs, label = '%.0f$^\circ$' % (angles[i] * 57.2957795))
		#	i += 1

		#legend = ax.legend()
		#for label in legend.get_texts():
		#	label.set_fontsize('small')

		#plt.ylabel('PTM [V]', fontdict=font)
		#plt.xlabel('Tiempo [s]', fontdict=font)
		#plt.xlim(0, MAX_LOG_SEG)

		#filename = folder + 'itv-log.png'
		#pylab.savefig(filename,	bbox_inches='tight')
		#print filename
