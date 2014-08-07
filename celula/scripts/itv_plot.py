# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import os

#para generar graficos de itv vs tiempo de itv.csv

MAXTIME = 20e-3
NLINEAS = 7
FOLDERS = [		# (ruta, radio, cant(angulos))
	#("../salida/10-64/", 10, 63),
	#("../salida/25-64/", 25, 64),
	#("../salida/50-64/", 50, 64),
	#("../salida/25-128/", 25, 128),
	#('../salida/ph7/25_64_cls/', 25, 64),
	('../salida/ph7/10_d/', 10, 64),
]

font = { 'size': 16 }

for (folder, radio, mesh) in FOLDERS:
	#for dir in os.listdir(folder):
	#	if dir[-3:] == ".7z" or dir[-4:] == ".tmp": continue
	#	file = "%s%s/itv.csv" % (folder, dir)
	
		file = '%s/itv.csv' % folder
		
		times = []
		itvss = [[] for x in range(NLINEAS)]
		iis = [int(round(((mesh-1) * i) / (NLINEAS-1), 0)) for i in range(NLINEAS)]
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

			#import pdb; pdb.set_trace()

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
			plt.xlim(0, 10)
			
			pylab.savefig(
				 #'itvs/itv-time-lin-%s-%s-%s.png' % (radio, mesh, dir), 
				 'itvs/ph7/lin-%s-%s-%s.png' % (radio, mesh, '50KVm'), 
				 bbox_inches='tight'
			)
			print "ploted"

			plt.clf()
			i = 0
			fig, ax = plt.subplots()
			for itvs in itvss:
				assert len(itvs) == len(times)
				itvs = np.array(itvs)
				plt.semilogx(times, itvs, label = '%.0f$^\circ$' % (angles[i] * 57.2957795))
				i += 1

			legend = ax.legend()
			for label in legend.get_texts():
				label.set_fontsize('small')

			plt.ylabel('ITV [V]', fontdict=font)
			plt.xlabel('Tiempo [s]', fontdict=font)
			plt.xlim(0, 20e-3)
			pylab.savefig(
				 #'itvs/itv-time-log-%s-%s-%s.png' % (radio, mesh, dir), 
				 'itvs/ph7/log-%s-%s-%s.png' % (radio, mesh, '50KVm'), 
				 bbox_inches='tight'
			)
			print "ploted"
