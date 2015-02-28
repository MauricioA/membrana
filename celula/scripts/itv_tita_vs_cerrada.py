# -*- coding: utf-8 -*-

from __future__ import unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import math
from pylab import *
import pylab
import sys
import os

folder = "../salida/pulsativo/02-27/poisson-1600Vcm/"

font = { 'size': 16 }

# real

file = "%s/itv.csv" % folder
angles_exp = None
itvs_exp = None

with open(file) as f:
	angles = []
	itvs = []

	for line in f:
		spl = line.split(',')
		time = float(spl[0])
		angle = float(spl[1])
		itv = float(spl[2])

		angles.append(angle)
		itvs.append(itv)

	angles_exp = np.array(angles)
	itvs_exp = np.array(itvs)

# formula cerrada

angles_teo = []
itvs_teo = []

lo = 0.15
li = 0.2
lm = 5e-6
d = 5e-9
R = 25e-6
num = num = 3 * lo * (3*d*R*R*li + (3*d*d*R - d**3) * (lm-li))
den = 2 * R**3 * (lm+2*lo) * (lm + .5 * li) - 2* (R-d)**3 * (lo-lm) * (li-lm)
fs = num / den
E = 160e3

tita = 0
while tita < math.pi:
	itv = fs * E * R * cos(tita)
	angles_teo.append(tita)
	itvs_teo.append(itv)
	tita += math.pi / 256
		
plt.clf()
fig, ax = plt.subplots()
fig.set_size_inches(16, 6)

ax.plot(angles_exp, itvs_exp, 'ro', label = u'simulaci\u00f3n')
ax.plot(angles_teo, itvs_teo, label = u'f\u00f3rmula', linewidth=2)

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
