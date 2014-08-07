import math

r = 10
h = 100
d = 5e-3

n_curva = 512
n_vertices = (n_curva+1) * 2 + 4

print "%d %d %d %d" % (n_vertices, 2, 0, 0)

dtheta = math.pi/n_curva
theta = []
x = -math.pi/2
for i in range(n_curva):
	theta.append(x)
	x += dtheta

theta.append(math.pi/2)

v = 0

for i in range(n_curva+1):
	v += 1
	print "%d %f %f" % (v, math.cos(theta[i]) * r, math.sin(theta[i]) * r + h/2)

for i in range(n_curva+1):
	v += 1
	print "%d %f %f" % (v, math.cos(theta[i]) * (r + d), math.sin(theta[i]) * (r + d) + h/2)

v += 1
print "%d %f %f" % (v, 0., 0.)

v += 1
print "%d %f %f" % (v, 50., 0.)

v += 1
print "%d %f %f" % (v, 50., 100.)

v += 1
print "%d %f %f" % (v, 0., 100.)

segs_i = []
segs_o = []
for i in range(1, n_curva+1):
	segs_i.append((i, i+1))
	segs_o.append((i + n_curva + 1, i + n_curva + 2))

segs_otros = [(v, v-1), (v-1, v-2), (v-2, v-3), (v-3, n_curva+2), (n_curva+2, 1), (1, n_curva+1), (n_curva+1, v-4), (v-4, v)]

segs_tot = segs_i + segs_o + segs_otros

print "%d, %d" % (len(segs_tot), 0)

for s in range(len(segs_tot)):
	print "%d %d %d" % (s+1, segs_tot[s][0], segs_tot[s][1])
	
print 0

#~ print segs_otros
