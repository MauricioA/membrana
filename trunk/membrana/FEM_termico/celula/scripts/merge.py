import argparse, re

#Para convertir del merged generado por Automesh2d al formato usado y cambiar escala

ENTRADA = "50_128_merged.msh"
SALIDA  = "50_128.msh"
ESCALA  = 100
RADIO   = 50
ANCHO   = 250
ALTO    = 500
MEMB    = 5e-3

centerX = 0
centerY = ALTO / 2
toler   = 1e-3

msh = open(ENTRADA)
lines = msh.readlines()
msh.close()

nPuntos = int(lines[0])
puntos = [None]

for i in range(1, nPuntos+1):
	line = lines[i]
	reg = re.match('\d+\s+(\d+\.?\d*)\s+(\d+\.?\d*)', line)
	puntos.append((
		float(reg.groups()[0]) / ESCALA, 
		float(reg.groups()[1]) / ESCALA
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
		if dist > (RADIO + MEMB + toler):
			dominio = externos
		elif dist < (RADIO - toler):
			dominio = internos

	dominio.append(ipuntos)

out = open(SALIDA, "w")

out.write("%s\n" % nPuntos)

for i in range(1, nPuntos+1):
	out.write('%s %s %s\n' % (i, puntos[i][0], puntos[i][1]))

out.write('%s\n' % len(doms))

doms = [externos, membrana, internos]

for i in range(len(doms)):
	out.write('%s %s\n' % (i+1, len(doms[i])-1))

for dom in doms:
	for i in range(1, len(dom)):
		out.write('%s %s %s %s %s\n' % (i, dom[i][0], dom[i][1], dom[i][2], dom[i][3]))

out.close()