import math

e = 40e3	# [V/m] (100v/cm)
r = 50e-6	# [m]
d = 5e-9	

so = 0.15	# [S/m]
sm = 5e-7
si = 0.2

num = 3*so * (3*d*(r**2)*si + (3*(d**2)*r - (d**3)) * (sm - si))
den = 2*(r**3) * (sm + 2*so) * (sm + .5*si) - 2 * ((r - d)**3) * (so - sm) * (si - sm)
f = num / den

def itv_teo(tita):
	return f * e * r * math.cos(tita)

n_0 = 1.5e9
r_m = 0.8e-9
r_ast = 0.51e-9
q = (r_m / r_ast) ** 2
v_ep = 0.258
alfa = 1e9

n = 128

densidades = []		# tita, densidad, itv, neq
t = math.pi/n
while t < math.pi:
	densidades.append((t, 0, itv_teo(t), n_0 * math.exp(q * (itv_teo(t) / v_ep) ** 2)))
	t += math.pi/n

time = 0.
delta_t = 100e-9
iter = 0
while time < 1:
	densidades_nueva = []
	densidad_prom = 0
	
	for (tita, densidad, vm, neq) in densidades:
		dens = densidad + delta_t * alfa * math.exp((vm/v_ep) ** 2) * (1 - densidad  / neq)	
		densidades_nueva.append((tita, dens, vm, neq))
		densidad_prom += dens

	densidades = densidades_nueva
	densidad_prom /= len(densidades)
	if iter % 50000 == 0:
		print '%.2e %.5e' % (time, densidad_prom * (1e6)**-2)
	iter += 1
	time += delta_t
