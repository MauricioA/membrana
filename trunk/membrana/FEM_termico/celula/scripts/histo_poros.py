import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys
import os

#para generar histogramas de poros sobre poros.csv

# (ruta, radio, cant(angulos))
FOLDERS = [		
	("../salida/10-64/", 10, 63),
	("../salida/25-64/", 25, 64),
	("../salida/50-64/", 50, 64),
	("../salida/25-128/", 25, 128),
]

TIMES = [(5e-6, '5e-6'), (10e-6, '10e-6'), (100e-6, '100e-6'), (1e-3, '1e-3'), (5e-3, '5e-3')]
font = {'size': 16}

for (folder, radio, mesh) in FOLDERS:
	for dir in os.listdir(folder):
		if dir[-3:] == ".7z" or dir[-4:] == ".tmp": continue
		file = "%s%s/poros.csv" % (folder, dir)
		
		for (TIME, time_name) in TIMES:
			found = False
			currTime = None
			radios = []

			with open(file) as f:
				for line in f:
					spl = line.split(',')
					time  = float(spl[0])
					angle = float(spl[1])
					radiop= float(spl[2])

					if not found:
						if time >= TIME:
							found = True
							currTime = time
					else:
						if abs(time - currTime) > 1e-12:
							break;

					if found and radiop > 1e-3:
						radios.append(radiop)

				if len(radios) == 0: continue
				radios = np.array(radios)

				plt.clf()
			
				pylab.hist(radios * 1e3, 100)
				plt.xlabel('Radio [nm]', fontdict=font)
				pylab.savefig('poros/hist-radios-%s-%s-%s-%s.png' % (time_name, radio, mesh, dir))
				
				sys.stdout.write('.')
				sys.stdout.flush()

sys.stdout.write('\n')
sys.stdout.flush()
