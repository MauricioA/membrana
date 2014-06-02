import matplotlib.pyplot as plt
import numpy as np
import pylab
import os

#para generar graficos de itv vs tiempo de itv.csv

MAXTIME = 20e-3
FOLDERS = [
	#"../salida/10-64/",
	"../salida/25-64/",
]

ANGLES = [
	0.0245429305,
	0.564503713,
	1.05537859,
	1.59533926,
	2.08621406,
	2.62617492,
	3.11704972,
]

for folder in FOLDERS:
	for dir in os.listdir(folder):
		if dir[-3:] == ".7z" or dir[-4:] == ".tmp": continue
		file = "%s%s/itv.csv" % (folder, dir)
		
		times = []
		itvss = [[] for x in ANGLES]
		
		with open(file) as f:
			for line in f:
				spl = line.split(',')
				time = float(spl[0])
				angle = float(spl[1])
				itv = float(spl[2])
				if time > MAXTIME: break
				for i in range(len(ANGLES)):
					ang = ANGLES[i]
					if abs(angle - ang) < 1e-6:
						itvss[i].append(itv)
						if i == 0: times.append(time)

			times = np.array(times)

			plt.clf()
			for itvs in itvss:
				assert len(itvs) == len(times)
				itvs = np.array(itvs)
				plt.plot(times*1e6, itvs)

			plt.ylabel('TMV [V]')
			plt.xlabel('Tiempo [us]')
			plt.xlim(0, 10)
			pylab.savefig('itvs/itv-time-lin-%s-%s.png' % (folder, dir), bbox_inches='tight')

			plt.clf()
			for itvs in itvss:
				assert len(itvs) == len(times)
				itvs = np.array(itvs)
				plt.semilogx(times*1e3, itvs)

			plt.ylabel('TMV [V]')
			plt.xlabel('Tiempo [ms]')
			plt.xlim(0, 20e-3)
			pylab.savefig('itvs/itv-time-log-%s-%s.png' % (folder, dir), bbox_inches='tight')
			

