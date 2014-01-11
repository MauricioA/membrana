import argparse, re, math

parser = argparse.ArgumentParser()
parser.add_argument('fem')
parser.add_argument('results')
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
	
voltajes = open(args.results)
lines = voltajes.readlines()
voltajes.close()
tolerancia = 10e-3
voltajes = []

for line in lines[1:]:
	reg = re.match('\s+(\d+\.?\d*E?\+?\-?\d*)\,\s+(\d+\.?\d*E?\+?\-?\d*)\,\s+(\d+\.?\d*E?\+?\-?\d*)', line)
	groups = reg.groups()
	x = float(groups[0])
	y = float(groups[1])
	v = float(groups[2])
	for punto in puntosMemb:
		if abs(punto[0] - x) < tolerancia and abs(punto[1] - y) < tolerancia:
			voltajes.append((x, y, v))
			break

print len(voltajes)

centerX = 0
centerY = 5
angulos = []	# (grados, voltaje)

for voltaje in voltajes:
	x = voltaje[0]
	y = voltaje[1]
	v = voltaje[2]
	radio = math.sqrt((x - centerX) ** 2 + (y - centerY) ** 2)
	if y < 5:
		tita = math.degrees(math.asin((x - centerX) / radio))
	else:
		tita = math.degrees(math.asin((y - centerY) / radio)) + 90
	angulos.append((tita, v))
	
angulos = sorted(angulos, key = lambda ang : ang[0])

print 'tita [grados], phi [V]'
i = 0;
while i < len(angulos):
	ang = angulos[i]
	print str(ang[0]) + ', ' + str(ang[1])
	i += 1
