import numpy as np
import matplotlib.pyplot as plt
import pylab

E = 120e3
radio = 25e-6
font = { 'size': 16 }

x = np.linspace(0, np.pi, num=1000)
y = 1.5 * E * radio * np.cos(x)

fig, ax = plt.subplots()

plt.plot(x, y)

ax.set_ylabel('PTM [V]', fontdict=font)
ax.set_xlabel('$\\theta$ [rad]', fontdict=font)
ax.set_xlim(0, np.pi)

#plt.show()

pylab.savefig(
	'itv-cos.png', 
	bbox_inches='tight'
)
