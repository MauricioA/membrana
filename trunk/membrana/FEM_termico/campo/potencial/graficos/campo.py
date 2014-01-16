import argparse, re, math

centerX 	= 0
centerY 	= 50
radioCel 	= 10
dMemb 		= 25e-3	   # en um
tolerancia 	= 1e-3

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
valores = []	# (grados, valor)

minDist = radioCel - tolerancia
maxDist = radioCel + dMemb + tolerancia

for line in lines[1:]:
	reg = re.match('\s+(\d+\.?\d*E?\+?\-?\d*)\,\s+(\d+\.?\d*E?\+?\-?\d*)\,\s+(\d+\.?\d*E?\+?\-?\d*)', line)
	groups = reg.groups()
	x = float(groups[0])
	y = float(groups[1])
	c = float(groups[2])
	dist = math.sqrt((x - centerX)**2 + (y - centerY)**2)
	if dist > minDist and dist < maxDist:
		radio = math.sqrt((x - centerX) ** 2 + (y - centerY) ** 2)
		tita = math.degrees(math.acos((y - centerY) / radio))
		valores.append((tita, c))
	
valores = sorted(valores, key = lambda ang : ang[0])

print len(valores)

f_salida = open(args.salida, "w")

f_salida.write("tita, campo[?]\n")

for valor in valores:
	f_salida.write("%f, %e\n" % (valor[0], valor[1]))
	
f_salida.close()
