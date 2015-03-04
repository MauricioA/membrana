import numpy as np
from math import log, exp
import pylab
import matplotlib.pyplot as plt

mp = 2e-6
k = -500e-6 / log(0.1)

x = np.linspace(0, 1.5e-3, num = 300)
y = [mp * (1 - exp(- t / k)) for t in x]

fig, ax = plt.subplots()
fig.set_size_inches(12, 4.5)

plt.plot(x * 1e3, np.array(y) * 1e9, '-', linewidth=2)

plt.xlabel('Tiempo [ms]')
plt.ylabel('$\Delta_t$ [ns]')

#plt.show()

outfile = '../../latex/template/graficos/acoplado/deltat.png'
pylab.savefig(outfile, bbox_inches='tight')
