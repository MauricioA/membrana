import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument('msh')
args = parser.parse_args()

centerX = 0
centerY = 50
radio = 10
memb = 5e-3
toler = 1e-3

msh = open(args.msh)
lines = msh.readlines()
msh.close()

nPuntos = int(lines[0])
puntos = [None]

for i in range(1, nPuntos+1):
	line = lines[i]
	reg = re.match('\d+\s+(\d+\.?\d*)\s+(\d+\.?\d*)', line)
	puntos.append((
		float(reg.groups()[0]) / 10, 
		float(reg.groups()[1]) / 10
	))

nElems = int(lines[nPuntos+1])

externos = [None]
membrana = [None]
internos = [None]

for i in range(nPuntos+2, nPuntos+2+nElems):
	line = lines[i]
	reg = re.match('\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+', line)
	ipuntos = [int(reg.groups()[j]) for j in range(0, 4)]
	nodosElem = [puntos[j] for j in ipuntos]
	dominio = membrana
	for nodo in nodosElem:
		dist = ((nodo[0] - centerX) ** 2 + (nodo[1] - centerY) ** 2) ** .5
		if dist > (radio + memb + toler):
			dominio = externos
		elif dist < (radio - toler):
			dominio = internos
	dominio.append(ipuntos)

print nPuntos

for i in range(1, nPuntos+1):
	print '%s\t%s\t%s' % (i, puntos[i][0], puntos[i][1])

doms = [externos, membrana, internos]

for i in range(len(doms)):
	print '%s\t%s' % (i+1, len(doms[i])-1)

for dom in doms:
	for i in range(1, len(dom)):
		print '%s\t%s\t%s\t%s\t%s' % (i, dom[i][0], dom[i][1], dom[i][2], dom[i][3])
