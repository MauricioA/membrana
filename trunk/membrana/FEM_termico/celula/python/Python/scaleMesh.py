import re

SCALE_X = 1e-6
SCALE_Y = 1e-6
ENTRADA = "grande.msh"
SALIDA = "salida.msh"

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
