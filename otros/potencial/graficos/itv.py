import argparse, re, math

centerX 	= 0
centerY 	= 50
radioCel 	= 10
dMemb 		= 5e-3		# en um
tolerancia 	= 1e-3

parser = argparse.ArgumentParser()
# parser.add_argument('fem')
parser.add_argument('datos')
parser.add_argument('salida')
args = parser.parse_args()

# fem = open(args.fem)
# lines = fem.readlines()
# fem.close()

# i = 0
# while True:
# 	line = lines[i]
# 	if line[0] != '*':
# 		break
# 	i += 1

# nPuntos = int(lines[i])
# i += 1
# puntos = [None]

# p = 0
# while p < nPuntos:
# 	line = lines[i]
# 	reg = re.match('(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)', line)
# 	i += 1
# 	if reg.groups():
# 		assert float(reg.groups()[0]) == p+1
# 		puntos.append((
# 			float(reg.groups()[1]), 
# 			float(reg.groups()[2])
# 		))
# 		p += 1

# assert nPuntos+1 == len(puntos)

# while True:
# 	line = lines[i]
# 	if line[0] != '*':
# 		break
# 	i += 1

# i += 1
# nElemsGrupo1 = int(re.match('\d+\s+(\d+)', lines[i]).groups()[0])
# i += 1
# nElemsMemb   = int(re.match('\d+\s+(\d+)', lines[i]).groups()[0])
# puntosMemb = set()
# i += 1

# while True:
# 	line = lines[i]
# 	if line[0] != '*':
# 		break
# 	i += 1

# i += nElemsGrupo1
# e = 0

# while p < nElemsMemb: 
# 	line = lines[i]
# 	# groups = re.match('(\d+)\s+(\d+)\s+(\d+)\s(\d+)', line).groups()	#tri
# 	groups = re.match('(\d+)\s+(\d+)\s+(\d+)\s(\d+)\s(\d+)', line).groups()	#quad
# 	# for j in range(1, 4):	#tri
# 	assert int(groups[0] == e+1)
# 	for j in range(0, 4):	#quad
# 		puntosMemb.add(puntos[int(groups[j+1])])
# 	e += 1
	
resultados = open(args.datos)
lines = resultados.readlines()
resultados.close()
valores = []

internos = []	# (grados, valor)
externos = []

for line in lines[1:]:
	reg = re.match('\s*(\-?\d+\.?\d*E?e?\+?\-?\d*)\,\s*(\-?\d+\.?\d*E?e?\+?\-?\d*)\,\s*(\-?\d+\.?\d*E?e?\+?\-?\d*)', line)
	groups = reg.groups()
	x = float(groups[0])
	y = float(groups[1])
	v = float(groups[2])

	radio = math.sqrt((x - centerX) ** 2 + (y - centerY) ** 2)
	tita = math.degrees(math.acos((y - centerY) / radio))

	if abs(radio - radioCel) < tolerancia:
		internos.append((tita, v))
	elif abs(radio - (radioCel + dMemb)) < tolerancia:
		externos.append((tita, v))

internos = sorted(internos, key = lambda ang : ang[0])
externos = sorted(externos, key = lambda ang : ang[0])


# 	for punto in puntosMemb:
# 		if abs(punto[0] - x) < tolerancia and abs(punto[1] - y) < tolerancia:
# 			valores.append((x, y, v))
# 			break


# for valor in valores:
# 	x = valor[0]
# 	y = valor[1]
# 	v = valor[2]
# 	radio = math.sqrt((x - centerX) ** 2 + (y - centerY) ** 2)
# 	tita = math.degrees(math.acos((y - centerY) / radio))
	
# 	if abs(radio - radioCel) < tolerancia:
# 		internos.append((tita, v))
# 	elif abs(radio - (radioCel + dMemb)) < tolerancia:
# 		externos.append((tita, v))
	
# internos = sorted(internos, key = lambda ang : ang[0])
# externos = sorted(externos, key = lambda ang : ang[0])

# print len(internos)
# print len(externos)

f_salida = open(args.salida, "w")

f_salida.write("tita, itv [V]\n")

v_promedio = 0

for i in range(len(internos)):
	tita1 = internos[i][0]
	tita2 = externos[i][0]
	assert abs(tita1 - tita2) < 1/60.
	tita = (tita1 + tita2) / 2
	valor = externos[i][1] - internos[i][1]
	f_salida.write("%f, %e\n" % (tita, valor))
	v_promedio += valor
	
f_salida.close()

v_promedio /= len(internos)

print "promedio = %e" % v_promedio
