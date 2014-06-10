import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys
import os

#para generar graficos de itv vs angulo de los itv.csv

FOLDERS = [		# (ruta, radio, cant(angulos))
	("../salida/10-64/", 10, 63),
	("../salida/25-64/", 25, 64),
	("../salida/50-64/", 50, 64),
	("../salida/25-128/", 25, 128),
]

TIMES = [1e-6, 2e-6, 4e-6, 16e-6]

for (folder, radio, mesh) in FOLDERS:
	for dir in os.listdir(folder):
		if dir[-3:] == ".7z" or dir[-4:] == ".tmp": continue
		file = "%s%s/itv.csv" % (folder, dir)
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
				ax.plot(angles[i], itvs[i], label = '%.0fus' % (TIMES[i]*1e6))

			legend = ax.legend()
			#for label in legend.get_texts():
			#	label.set_fonstize('small')

			ax.set_ylabel('TMV [V]')
			ax.set_xlabel('Tita')
			ax.set_xlim(0, np.pi)

			pylab.savefig(
				 'itvs/itv-tita-%s-%s-%s.png' % (radio, mesh, dir), 
				 bbox_inches='tight'
			)

			sys.stdout.write('.')
			sys.stdout.flush()

sys.stdout.write('\n')
