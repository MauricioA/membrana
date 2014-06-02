import matplotlib.pyplot as plt
import numpy as np

#para generar graficos de itv vs tiempo de itv.csv

MAXTIME = 5e-3
FILE = "../salida/25-64/40KVm/itv.csv"
ANGLES = [
	0.0245429305,
	0.662680648,
	1.20264138,
	1.79168882,
	2.3316501,
	3.11704972,
]

times = []
itvss = [[] for x in ANGLES]

with open(FILE) as f:
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

for itvs in itvss:
	assert len(itvs) == len(times)
	itvs = np.array(itvs)
	plt.plot(times, itvs)

plt.xscale('log')
plt.show()

#import matplotlib.pyplot as plt
#import numpy as np

#times = []
#itvs0 = []
#itvs1 = []
#itvs2 = []

#MAXTIME = 4e-6

#with open("itv.dat") as f:
#	i = 0
#	n = 0
#	for line in f:
#		spl = line.split()
#		if i == n:
#			n = int(spl[2])
#			time = float(spl[1])
#			if time > MAXTIME: break
#			times.append(time)
#			assert n == 63
#			i = 0
#		else:
#			if i == 0:
#				itvs0.append(float(spl[1]))
#			if i == 1:
#				itvs1.append(float(spl[1]))
#			if i == 2:
#				itvs2.append(float(spl[1]))
#			i += 1

#assert(len(times) == len(itvs0))

#times = np.array(times)
#itvs0 = np.array(itvs0)
#itvs1 = np.array(itvs1)
#itvs2 = np.array(itvs2)

#plt.plot(times, itvs0)
#plt.plot(times, itvs1)
#plt.plot(times, itvs2)
#plt.show()

1e-009,0.0245429305,0
1e-009,0.0736299341,0
1e-009,0.122718049,0
1e-009,0.171805299,0
1e-009,0.220892493,0
1e-009,0.269980961,0
1e-009,0.319068053,0
1e-009,0.368154946,0
1e-009,0.417242166,0
1e-009,0.466330301,0
1e-009,0.515417736,0
1e-009,0.564503713,0
1e-009,0.613592405,0
1e-009,0.662680648,0
1e-009,0.711766658,0
1e-009,0.760853772,0
1e-009,0.809942555,0
1e-009,0.859029668,0
1e-009,0.908115678,0
1e-009,0.957203921,0
1e-009,1.00629261,0
1e-009,1.05537859,0
1e-009,1.10446603,0
1e-009,1.15355416,0
1e-009,1.20264138,0
1e-009,1.25172827,0
1e-009,1.30081537,0
1e-009,1.34990383,0
1e-009,1.39899103,0
1e-009,1.44807828,0
1e-009,1.49716639,0
1e-009,1.5462534,0
1e-009,1.59533926,-0
1e-009,1.64442626,-0
1e-009,1.69351438,-0
1e-009,1.74260163,-0
1e-009,1.79168882,-0
1e-009,1.84077729,-0
1e-009,1.88986438,-0
1e-009,1.93895127,-0
1e-009,1.98803849,-0
1e-009,2.03712663,-0
1e-009,2.08621406,-0
1e-009,2.13530004,-0
1e-009,2.18438873,-0
1e-009,2.23347698,-0
1e-009,2.28256299,-0
1e-009,2.3316501,-0
1e-009,2.38073888,-0
1e-009,2.429826,-0
1e-009,2.47891201,-0
1e-009,2.52800025,-0
1e-009,2.57708894,-0
1e-009,2.62617492,-0
1e-009,2.67526235,-0
1e-009,2.72435049,-0
1e-009,2.77343771,-0
1e-009,2.8225246,-0
1e-009,2.87161169,-0
1e-009,2.92070016,-0
1e-009,2.96978735,-0
1e-009,3.0188746,-0
1e-009,3.06796272,-0
1e-009,3.11704972,-0
