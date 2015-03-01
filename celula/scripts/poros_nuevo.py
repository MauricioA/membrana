# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pylab
import sys
import os

#para generar histogramas de poros

folders = [
	'../salida/pulsativo/02-17/longoff/',
]

#ESTOS = [
#	'poros-1-0.000100.csv',
#	'poros-101-0.010100.csv',
#	'poros-201-0.020100.csv',
#	'poros-301-0.030100.csv',
#]

font = {'size': 16}

#MIN_RADIO = 1e-3
MIN_RADIO = 0

MIN_X = 0
MAX_X = 70
MIN_Y = 0
MAX_Y = 200

dir_salida = "/graficosnew/"

for folder in folders:
	poros_folder = folder + 'poros/'
	files = filter(lambda x : x.endswith('.csv'), os.listdir(poros_folder))

	try: os.mkdir(poros_folder + dir_salida)
	except: pass

	for file in files:
		radios = []

		#if not file in ESTOS: continue

		with open(poros_folder + file) as f:
			for line in f.readlines()[1:]:
				tita, radio = [float(x) for x in line.split(',')]
				if radio > MIN_RADIO:
					radios.append(radio)
			
		if len(radios) == 0: continue
		radios = np.array(radios)
		plt.clf()
		
		#descomentar para escala automática
		fig, ax = plt.subplots()
		ax.set_xlim(MIN_X, MAX_X)
		ax.set_ylim(MIN_Y, MAX_Y)

		pylab.hist(radios * 1e3, 100)
		plt.xlabel('Radio [nm]', fontdict=font)
		filename = poros_folder + dir_salida + file[:-len('.csv')] + '.png'
		pylab.savefig(filename)
		print filename
