ARCHIVO = "densidad.dat"

area = 0
lines = 0

with open(ARCHIVO) as f:
    f.readline()
    for line in f:
        t, a = [float(x) for x in line.split()]
        area += a
        lines += 1

print area
print lines
