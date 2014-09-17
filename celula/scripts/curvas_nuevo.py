import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import pylab
import sys, os, pdb

input_folder = '../salida/poster64/dt25/150kvm/transporte/'
font = { 'size': 16 }
labels = { 'H': 'H$^+$', 'OH': 'OH$^-$', 'Na': 'Na$^+$', 'Cl': 'Cl$^-$' }

#retorna un dicc con np.array de Y y dos especies segun el archivo
def getValues(archivo):
	ys = []
	c1s = []
	c2s = []
	
	with open(archivo) as f:
		for line in f.readlines()[1:]:
			numbers = [float(x) for x in line.split(',')]
			x = numbers[0]
			y = numbers[1]
			c1 = numbers[2]
			c2 = numbers[3]
			
			if abs(x) < 1e-6:
				ys.append(y)
				c1s.append(c1)
				c2s.append(c2)

	if 'concent' in archivo:
		c1_esp = 'Na'
		c2_esp = 'Cl'
	elif 'pH' in archivo or 'hidro' in archivo:
		c1_esp = 'H'
		c2_esp = 'OH'
	else: 
		assert False

	return {
		'y': np.array(ys),
		c1_esp: np.array(c1s),
		c2_esp: np.array(c2s),
	}

iniciales_ph  = getValues(input_folder + 'pH-0-0.000000.csv')
iniciales_con = getValues(input_folder + 'concent-0-0.000000.csv')
iniciales_hid = getValues(input_folder + 'hidro-0-0.000000.csv')

out_dir = input_folder + 'curvas'

try: os.mkdir(out_dir) 
except: pass

files = [
	x for x in os.listdir(input_folder) 
	if (x.startswith('concent') or x.startswith('pH') or x.startswith('hidro')) and 
	not (x.startswith('concent-0-') or x.startswith('pH-0-') or x.startswith('hidro-0-'))
]

for file in files:
	new_values = getValues(input_folder + file)
	
	if file.startswith('pH'):
		dicc_inicial = iniciales_ph
		esps = ['H', 'OH']
		outfile = '-pH' + file[len('pH'):-len('.csv')] + '.png'
	elif file.startswith('concent'):
		dicc_inicial = iniciales_con
		esps = ['Na', 'Cl']
		outfile = file[len('concent'):-len('.csv')] + '.png'
	elif file.startswith('hidro'):
		dicc_inicial = iniciales_hid
		esps = ['H', 'OH']
		outfile = '-hidro' + file[len('hidro'):-len('.csv')] + '.png'
	else: 
		assert False

	for esp in esps:
		#plt.clf()
		fig, ax = plt.subplots()
		
		ax.plot(dicc_inicial['y'], dicc_inicial[esp], 'bo')
		ax.plot(new_values['y'], new_values[esp], 'ro')
		
		#plt.xlabel('???', fontdict=font)
		#plt.ylabel('???', fontdict=font)

		#max = np.amax(ys)
		#plt.axvline(x=max/5 * 2, color="gray", ls='--')
		#plt.axvline(x=max/5 * 3, color="gray", ls='--')

		#pylab.show()

		pylab.savefig(out_dir + '/' + esp + outfile)
		print esp + outfile
