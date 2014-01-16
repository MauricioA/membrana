import argparse, re, math

centerX 	= 0
centerY 	= 50
radioCel 	= 10
dMemb 		= 25e-3		# en um
tolerancia 	= 1e-3
#~ sigma_memb  = 500e-15

#~ 25 4.864341e-06
#~ 50 2.183236e-06
#~ 75 2.304688e-06

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

internos = []	# (grados, voltaje)
externos = []

for voltaje in voltajes:
	x = voltaje[0]
	y = voltaje[1]
	v = voltaje[2]
	radio = math.sqrt((x - centerX) ** 2 + (y - centerY) ** 2)
	if y < 5:
		tita = 180 - math.degrees(math.asin((x - centerX) / radio))
	else:
		tita = 180 - math.degrees(math.asin((y - centerY) / radio)) - 90
	
	if abs(radio - radioCel) < tolerancia:
		internos.append((tita, v))
	elif abs(radio - (radioCel + dMemb)) < tolerancia:
		externos.append((tita, v))
	
internos = sorted(internos, key = lambda ang : ang[0])
externos = sorted(externos, key = lambda ang : ang[0])

print len(internos)
print len(externos)

volt_intra 	= open("itv.csv", "w")
#~ flujo_intra = open("flujo.csv", "w")

volt_intra.write("tita, itv [V]\n")
#~ flujo.write("angulo, flujo\n")

v_trans_promedio = 0

for i in range(len(internos)):
	tita = (internos[i][0] + externos[i][0]) / 2
	voltaje = externos[i][1] - internos[i][1]
	#~ flujo = sigma_memb * voltaje / dMemb
	volt_intra.write("%f, %e\n" % (tita, voltaje))
	#~ flujo.write("%f, %e\n" % (tita, flujo))
	v_trans_promedio += voltaje
	
volt_intra.close()
#~ flujo.close()

v_trans_promedio /= len(internos)

print "V promedio = %e" % v_trans_promedio
