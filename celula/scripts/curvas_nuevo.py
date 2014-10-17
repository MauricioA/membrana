import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys, os, pdb

input_folder = '../salida/temp/conduc/50e-12/transporte/'

font = { 'size': 16 }
labels = { 'H': 'H$^+$', 'OH': 'OH$^-$', 'Na': 'Na$^+$', 'Cl': 'Cl$^-$' }
x_curva = 0
tolerancia = 1e-6
nesps = 4
H, OH, Na, Cl = [i for i in range(4)]

def getValues(archivo):
	ys = []
	concs = [[] for x in range(nesps)]
	
	with open(archivo) as f:
		for line in f.readlines()[1:]:
			x, y, h, oh, na, cl = [float(x) for x in line.split(',')]
			
			if abs(x - x_curva) < tolerancia:
				ys.append(y)
				concs[H].append(h)
				concs[OH].append(oh)
				concs[Na].append(na)
				concs[Cl].append(cl)

	return {
		'y': np.array(ys),
		'H': np.array(concs[H]),
		'OH': np.array(concs[OH]),
		'Na': np.array(concs[Na]),
		'Cl': np.array(concs[Cl]),
	}

iniciales = getValues(input_folder + 'molar-0-0.000000.csv')

out_dir = input_folder + 'curvas'

try: os.mkdir(out_dir) 
except: pass

files = [x for x in os.listdir(input_folder) if (x.startswith('molar'))]

for file in files:
	nFile = int(file.split('-')[1])
	new_values = getValues(input_folder + file)
	esps = ['H', 'OH', 'Na', 'Cl']

	for esp in esps:
		fig, ax = plt.subplots()
		
		ax.plot(iniciales['y'], iniciales[esp], 'bo')
		ax.plot(new_values['y'], new_values[esp], 'ro')

		#plt.semilogy(new_values['y'], new_values[esp], 'ro')

		#plt.xlabel('', fontdict=font)
		#plt.ylabel('', fontdict=font)

		max = np.amax(iniciales['y'])
		plt.axvline(x=max/5 * 2, color="gray", ls='--')
		plt.axvline(x=max/5 * 3, color="gray", ls='--')

		outfile = '%s_%02d.png' % (esp, nFile)
		pylab.savefig('%s/%s' % (out_dir, outfile))
		print outfile
