ENTRADA = "../../salida/historia_12.dat"
SALIDA  = "../../salida/historia_12_resumen.dat"
FAC = 256

fi = open(ENTRADA)
fo = open(SALIDA, "w")
deb = 0
nNodes = int(fi.readline())
deb += 1
fo.write(str(nNodes) + "\n")
line = fi.readline()
deb += 1
n = FAC
last = line

while line != "":
    a = line.split()
    if len(a) != 2 or a[0] != "paso:":
        print "ERROR! " + str(deb)
        while len(a) != 2 or a[0] != "paso:":
            line = fi.readline()
            deb += 1
            a = line.split()

    time = float(a[1])
    do = n == FAC

    if do:
        fo.write(line)
        n = -1

    for i in range(nNodes):
        line = fi.readline()
        deb += 1
        if do:
            fo.write(line)

    line = fi.readline()
    deb += 1
    n += 1

fi.close()
fo.close()
