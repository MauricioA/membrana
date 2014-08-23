import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys, os, pdb

#para generar curvas de nivel de transporte sobre transporte.csv

FOLDERS = [		# (ruta, radio, cant(angulos))
	#("../salida/10-64/", 10, 63),
	#("../salida/25-64/", 25, 64),
	#("../salida/50-64/", 50, 64),
	#("../salida/25-128/", 25, 128),
	("../salida/ph7/25/", 25, 64),
	("../salida/ph7/25_128/", 25, 128),
]

ESPECIES = ['H', 'OH', 'NA', 'CL']
esps_latex = ['H$^+$', 'OH$^-$', 'Na$^+$', 'Cl$^-$']
font = { 'size': 16 }

for (folder, radio, mesh) in FOLDERS:
	#for dir in os.listdir(folder):
	#	if dir[-3:] == ".7z" or dir[-4:] == ".tmp": continue
	#	file = "%s%s/transporte.csv" % (folder, dir)
	
		file = '%s/transporte.csv' % folder

		with open(file) as f:
			first = next(f).decode()
			f.seek(-1024, 2)
			last = f.readlines()[-1].decode()

		firstTime = float(first.split(',')[0])
		lastTime  = float( last.split(',')[0])
		
		firstYs = []
		lastsYs = []
		firstCs = [[] for e in range(len(ESPECIES))]
		lastsCs = [[] for e in range(len(ESPECIES))]

		with open(file) as f:
			for line in f:
				spl = line.split(',')
				time  = float(spl[0])
				x = float(spl[1])
				y = float(spl[2])

				if abs(x - 0) > 1e-9: continue

				if abs(time - firstTime) < 1e-12:
					Ys = firstYs
					Cs = firstCs
				elif abs(time - lastTime) < 1e-12:
					Ys = lastsYs
					Cs = lastsCs
				else:
					continue

				#pdb.set_trace()

				Ys.append(y)
				for esp in range(len(ESPECIES)):
					Cs[esp].append(float(spl[esp+3]))

			#pdb.set_trace()

			firstYs = np.array(firstYs)
			lastsYs = np.array(lastsYs)
			firstCs = [np.array(x) for x in firstCs]
			lastsCs = [np.array(x) for x in lastsCs]
			
			for esp in range(len(ESPECIES)):
				plt.clf()
				fig, ax = plt.subplots()
				plt.semilogy(firstYs, firstCs[esp], 'bo')
				plt.semilogy(lastsYs, lastsCs[esp], 'ro')
				plt.xlabel('Y [$\\mu$m]', fontdict=font)
				plt.ylabel('%s [at $\\mu$m$^{-3}$]' % esps_latex[esp], fontdict=font)

				max = np.amax(firstYs)
				plt.axvline(x=max/5 * 2, color="gray", ls='--')
				plt.axvline(x=max/5 * 3, color="gray", ls='--')
				
				filename = 'curvas/ph7/%s-%s-%sKvm-%s.png' % (radio, mesh, 50, ESPECIES[esp])
				pylab.savefig(filename)
				sys.stdout.write('plotted %s\n' % filename)
				sys.stdout.flush()

sys.stdout.write('\n')
sys.stdout.flush()
