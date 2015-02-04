import matplotlib.pyplot as plt
import numpy as np
import pylab
import sys
import os

#para generar histogramas de poros

folders = [		
	'../salida/last/02-03/120/',
]

font = {'size': 16}
MIN_RADIO = 1e-3

for folder in folders:
	poros_folder = folder + 'poros/'
	files = filter(lambda x : x.endswith('.csv'), os.listdir(poros_folder))

	try: os.mkdir(poros_folder + 'graficos/')
	except: pass

	for file in files:
		radios = []

		with open(poros_folder + file) as f:
			for line in f.readlines()[1:]:
				tita, radio = [float(x) for x in line.split(',')]
				if radio > MIN_RADIO:
					radios.append(radio)
			
		if len(radios) == 0: continue
		radios = np.array(radios)
		plt.clf()
		
		pylab.hist(radios * 1e3, 100)
		plt.xlabel('Radio [nm]', fontdict=font)
		filename = poros_folder + '/graficos/' + file[:-len('.csv')] + '.png'
		pylab.savefig(filename)
		print filename
