import argparse, re, math

centerX 	= 0
centerY 	= 50
radioCel 	= 10
dMemb 		= 25e-3		# en um
tolerancia 	= 10e-3

parser = argparse.ArgumentParser()
parser.add_argument('fem')
parser.add_argument('datos')
parser.add_argument('salida')
args = parser.parse_args()

fem = open(args.fem)
lines = fem.readlines()
fem.close()

nPuntos = int(lines[1])
puntos = [None]

for i in range(2, nPuntos+2):
	line = lines[i]
	reg = re.match('(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)', line)
	puntos.append((
		float(reg.groups()[1]), 
		float(reg.groups()[2])
	))
	
nElemsGrupo1 = int(re.match('\d+\s+(\d+)', lines[nPuntos+4]).groups()[0])
nElemsMemb   = int(re.match('\d+\s+(\d+)', lines[nPuntos+5]).groups()[0])
puntosMemb = set()

for i in range(nPuntos + nElemsGrupo1 + 8, nPuntos + nElemsGrupo1 + 8 + nElemsMemb): 
	line = lines[i]
	groups = re.match('(\d+)\s+(\d+)\s+(\d+)\s(\d+)', line).groups()
	for j in range(1, 4):
		puntosMemb.add(puntos[int(groups[j])])
	
resultados = open(args.datos)
lines = resultados.readlines()
resultados.close()
valores = []

for line in lines[1:]:
	reg = re.match('\s+(\d+\.?\d*E?\+?\-?\d*)\,\s+(\d+\.?\d*E?\+?\-?\d*)\,\s+(\d+\.?\d*E?\+?\-?\d*)', line)
	groups = reg.groups()
	x = float(groups[0])
	y = float(groups[1])
	v = float(groups[2])
	for punto in puntosMemb:
		if abs(punto[0] - x) < tolerancia and abs(punto[1] - y) < tolerancia:
			valores.append((x, y, v))
			break

internos = []	# (grados, valor)
externos = []

for valor in valores:
	x = valor[0]
	y = valor[1]
	v = valor[2]
	radio = math.sqrt((x - centerX) ** 2 + (y - centerY) ** 2)
	tita = 180 - 90 - math.degrees(math.asin((y - centerY) / radio))
	
	if abs(radio - radioCel) < tolerancia:
		internos.append((tita, v))
	elif abs(radio - (radioCel + dMemb)) < tolerancia:
		externos.append((tita, v))
	
internos = sorted(internos, key = lambda ang : ang[0])
externos = sorted(externos, key = lambda ang : ang[0])

print len(internos)
print len(externos)

f_salida = open(args.salida, "w")

f_salida.write("tita, itv [V]\n")

v_promedio = 0

for i in range(len(internos)):
	tita = (internos[i][0] + externos[i][0]) / 2
	valor = externos[i][1] - internos[i][1]
	f_salida.write("%f, %e\n" % (tita, valor))
	v_promedio += valor
	
f_salida.close()

v_promedio /= len(internos)

print "promedio = %e" % v_promedio
