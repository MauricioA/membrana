import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys
import os

#para generar un graficos de itv vs angulo con varios resultados superpuestos

FOLDERS = [		# (ruta, radio, cant(angulos), )
	("../salida/pulsativo/02-27/poisson-radio-10/", 10, 192, "10 $\mu$m"),
]

font = { 'size': 16 }

angles = []
itvs = []

for i in range(len(FOLDERS)):
	folder, radio, mesh, _ = FOLDERS[i]
	file = "%s/itv.csv" % folder
	angles_i = []
	itvs_i = []

	with open(file) as f:
		for line in f:
			spl = line.split(',')
			time = float(spl[0])
			angle = float(spl[1])
			itv = float(spl[2])

			angles_i.append(angle)
			itvs_i.append(itv)

		angles.append(np.array(angles_i))
		itvs.append(np.array(itvs_i))
	
plt.clf()
		
fig, ax = plt.subplots()

for i in range(len(FOLDERS)):
	ax.plot(angles[i], itvs[i], label = FOLDERS[i][3])

legend = ax.legend()

#for label in legend.get_texts():
#	label.set_fonstize('small')

ax.set_ylabel('PTM [V]', fontdict=font)
ax.set_xlabel('$\\theta$ [rad]', fontdict=font)
ax.set_xlim(0, np.pi)

pylab.savefig(
	'temp.png',
	bbox_inches='tight'
)
