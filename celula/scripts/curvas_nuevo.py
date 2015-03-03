# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys, os

input_folder = '../salida/pulsativo/02-16/transporte/'

x_curva = 2.5
tolerancia = 1

#FILES = [
#	'molar.csv.038',
#]

FIRST_FILE = 'molar.csv.000'

xlim = 150

####

nesps = 4
H, OH, Na, Cl = [i for i in range(4)]
labels = { 'H': 'H$^+$', 'OH': 'OH$^-$', 'Na': 'Na$^+$', 'Cl': 'Cl$^-$' }
font = { 'size': 16 }

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

iniciales = getValues(input_folder + FIRST_FILE)

out_dir = input_folder + 'curvas'

try: os.mkdir(out_dir) 
except: pass

#files = [x for x in os.listdir(input_folder) if x in FILES]
files = [x for x in os.listdir(input_folder) if x.startswith('molar')]

for file in files:
	nFile = int(file.split('.')[2])
	new_values = getValues(input_folder + file)
	esps = ['H', 'OH', 'Na', 'Cl']

	for esp in esps:
		fig, ax = plt.subplots()
		fig.set_size_inches(12, 4.5)
		
		ax.semilogy(iniciales['y'], iniciales[esp], 'bo', label = 'Inicial')
		ax.semilogy(new_values['y'], new_values[esp], 'ro', label = 'Final')

		plt.xlabel('z', fontdict=font)
		plt.ylabel(labels[esp], fontdict=font)

		ax.set_xlim(0, xlim)

		# para poner líneas de borde célula
		#max = np.amax(iniciales['y'])
		#plt.axvline(x=max/5 * 2, color="gray", ls='--')
		#plt.axvline(x=max/5 * 3, color="gray", ls='--')
		plt.axvline(x=50, color="gray", ls='--')
		plt.axvline(x=100, color="gray", ls='--')

		ax.legend(loc = 'lower right')

		outfile = '%s-%02d.png' % (esp, nFile)
		pylab.savefig('%s/%s' % (out_dir, outfile), bbox_inches='tight')
		print outfile
