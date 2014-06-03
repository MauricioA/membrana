import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys
import os

#para generar histogramas de poros sobre poros.csv

FOLDERS = [		# (ruta, radio, cant(angulos))
	("../salida/10-64/", 10, 63),
	("../salida/25-64/", 25, 64),
	("../salida/50-64/", 50, 64),
]

TIMES = [30e-6, 100e-6, 1e-3]

for (folder, radio, mesh) in FOLDERS:
	for dir in os.listdir(folder):
		if dir[-3:] == ".7z" or dir[-4:] == ".tmp": continue
		file = "%s%s/poros.csv" % (folder, dir)
		found = False
		currTime = None
		radios = []
		
		for TIME in TIMES:
			with open(file) as f:
				for line in f:
					spl = line.split(',')
					time  = float(spl[0])
					angle = float(spl[1])
					radio = float(spl[2])

					if not found:
						if time >= TIME:
							found = True
							currTime = time
							#import pdb; pdb.set_trace()
					else:
						if abs(time - currTime) > 1e-12:
							break;

					if found:
						radios.append(radio)

				if len(radios) == 0: continue
				radios = np.array(radios)

				plt.clf()
			
				pylab.hist(radios, 20)
				pylab.savefig('poros/hist-radios-%s-%s-%s-%s.png' % (TIME, radio, mesh, dir))

				sys.stdout.write('.')
				sys.stdout.flush()

sys.stdout.write('\n')
sys.stdout.flush()

#			fig, ax = plt.subplots()

#			for i in range(len(TIMES)):
#				ax.plot(angles[i], itvs[i], label = '%.0fus' % (TIMES[i]*1e6))

#			legend = ax.legend()
#			#for label in legend.get_texts():
#			#	label.set_fonstize('small')

#			ax.set_ylabel('TMV [V]')
#			ax.set_xlabel('Tita')
#			ax.set_xlim(0, np.pi)

#			pylab.savefig(
#				 'itvs/itv-tita-%s-%s-%s.png' % (radio, mesh, dir), 
#				 bbox_inches='tight'
#			)

#			sys.stdout.write('.')
#			sys.stdout.flush()

#sys.stdout.write('\n')
