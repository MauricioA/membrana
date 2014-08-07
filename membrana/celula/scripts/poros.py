ENTRADA = "../../salida/poros.dat"
NPOROS  = "../../salida/nporos.dat"

time = 0
pasoPoro = 1e-9 * 500
fNPoros = open(NPOROS, "w")
depol = 0.024947
db = 0.872679
hb = 2.268914
h = 3.116646
maximos = []
f = open(ENTRADA)
lines = f.readlines()
ln = 0

while ln < len(lines):
    time += pasoPoro
    linea = lines[ln].split()
    ln += 1
    assert linea[0] == "paso"
    n = int(linea[2])
    fNPoros.write(str(time) + " " + str(n) + "\n")
    maxDep, maxdb, maxhb, maxh = (0, 0, 0, 0)
    for i in range(n):
        angulo, radio = [float(x) for x in lines[ln].split()]
        if abs(angulo - depol) < 1e-6 and radio > maxDep:
            maxDep = radio
        if abs(angulo - db) < 1e-6 and radio > maxdb:
            maxdb = radio
        if abs(angulo - hb) < 1e-6 and radio > maxhb:
            maxhb = radio
        if abs(angulo - h) < 1e-6 and radio > maxh:
            h = radio
        ln += 1
    m = (time, maxDep, maxdb, maxhb, maxh)
    maximos.append(m)
    print time
        
f.close()
fNPoros.close()

MAXIMO = "../../salida/maximo.dat"
f = open(MAXIMO, "w")

for m in maximos:
    f.write(str(m[0]) + " " + str(m[1]) + " " + str(m[2]) + " " + str(m[3]) + " " + str(m[4]) + "\n")

f.close()
