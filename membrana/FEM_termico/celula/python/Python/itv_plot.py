import matplotlib.pyplot as plt
import numpy as np

#para generar graficos de los archivos itv.dat

times = []
itvs0 = []
itvs1 = []
itvs2 = []

MAXTIME = 4e-6

with open("itv.dat") as f:
	i = 0
	n = 0
	for line in f:
		spl = line.split()
		if i == n:
			n = int(spl[2])
			time = float(spl[1])
			if time > MAXTIME: break
			times.append(time)
			assert n == 63
			i = 0
		else:
			if i == 0:
				itvs0.append(float(spl[1]))
			if i == 1:
				itvs1.append(float(spl[1]))
			if i == 2:
				itvs2.append(float(spl[1]))
			i += 1

assert(len(times) == len(itvs0))

times = np.array(times)
itvs0 = np.array(itvs0)
itvs1 = np.array(itvs1)
itvs2 = np.array(itvs2)

plt.plot(times, itvs0)
plt.plot(times, itvs1)
plt.plot(times, itvs2)
plt.show()
