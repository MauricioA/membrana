import re
#usar merge.py en vez de esto!

SCALE_X = 1. / 100	#CAMBIAR POR LO QUE CORRESPONDA!
SCALE_Y = 1. / 100
ENTRADA = "merged.msh"
SALIDA  = "scaled.msh"

fem = open(ENTRADA)
lines = fem.readlines()
fem.close()

out = open(SALIDA, "w")

i = 0
nPuntos = int(lines[i])
i += 1
out.write("%d\n" % nPuntos)

for p in range(nPuntos):
    reg = re.match('(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)', lines[i])
    nP = int(reg.groups()[0])
    x = float(reg.groups()[1])
    y = float(reg.groups()[2])

    assert nP == p+1
    out.write("%d %e %e\n" % (nP, x * SCALE_X, y * SCALE_Y))
    i += 1

while i < len(lines):
    out.write(lines[i])
    i += 1

out.close()
