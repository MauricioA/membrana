ENTRADA = "start.txt"
SALIDA = "start.0.txt"

fi = open(ENTRADA)
fo = open(SALIDA, "w")
line = fi.readline()
o = []

while line != "":
    a = line.split()
    if abs(float(a[1]) - 0) < 1e-9:
        o.append(line)
    line = fi.readline()

o = sorted(o, key = lambda x: float(x.split()[2]))

for x in o:
    fo.write(x)

fi.close()
fo.close()
