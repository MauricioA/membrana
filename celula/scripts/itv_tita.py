import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys
import os

#para generar graficos de itv vs angulo de los itv.csv

FOLDERS = [		# (ruta, radio, cant(angulos))
	('../salida/pulsativo/03-06/300vcm/', 25, 192),
	('../salida/pulsativo/03-06/400vcm/', 25, 192),
	('../salida/pulsativo/03-06/500vcm/', 25, 192),
	('../salida/pulsativo/03-06/800vcm/', 25, 192),
]

font = { 'size': 16 }

TIMES = [750e-9, 1.5e-6, 16e-6]
LABELS = ['750 ns', '1.5 $\\mu$s ', '16 $\\mu$s']

for (folder, radio, mesh) in FOLDERS:
	#for dir in os.listdir(folder):
	#	if dir[-3:] == ".7z" or dir[-4:] == ".tmp": continue
    #   file = "%s%s/itv.csv" % (folder, dir)
	file = "%s/itv.csv" % folder
	itime = 0
	found = False
	currTime = None

	angles = [[] for x in TIMES]
	itvs = [[] for x in TIMES]
		
	with open(file) as f:
		for line in f:
			spl = line.split(',')
			time = float(spl[0])
			angle = float(spl[1])
			itv = float(spl[2])

			if not found:
				if time >= TIMES[itime]:
					found = True
					currTime = time
			else:
				if abs(time - currTime) > 1e-12:
					found = False
					itime += 1
					if itime == len(TIMES): break

			if found:
				angles[itime].append(angle)
				itvs[itime].append(itv)

		#import pdb; pdb.set_trace()

		angles = np.array(angles)
		itvs = np.array(itvs)
		plt.clf()
		i = 0
			
		fig, ax = plt.subplots()

		for i in range(len(TIMES)):
			#ax.plot(angles[i], itvs[i], label = '%.0f$\\mu$s' % (TIMES[i]*1e6))
			ax.plot(angles[i], itvs[i], label = LABELS[i])

		legend = ax.legend()		#descomentar para poner label

		#for label in legend.get_texts():
		#	label.set_fonstize('small')

		ax.set_ylabel('PTM [V]', fontdict=font)
		ax.set_xlabel('$\\theta$ [rad]', fontdict=font)
		ax.set_xlim(0, np.pi)

		#plt.show()

		pylab.savefig(
			'%s/itv-tita-%s-%s.png' % (folder, radio, mesh), 
			bbox_inches='tight'
		)

		sys.stdout.write('.')
		sys.stdout.flush()

sys.stdout.write('\n')
