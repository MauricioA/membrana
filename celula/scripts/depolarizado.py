ENTRADA = "../../salida/poros.dat"
SALIDA = "../../salida/depolarizado.dat"
ANGULO = 3.116646e+000
TOLER = 10e-6

radios = []

with open(ENTRADA) as f:
    lines  = f.readlines()
    i = 0
    while i < len(lines):
        split = lines[i].split()
        t = float(split[0])
        n = int(split[1])
        rs = [t]
        for j in range(n):
            i += 1
            split = lines[i].split()
            if abs(float(split[0]) - ANGULO) < TOLER:
                rs.append(float(split[1]))
        radios.append((rs))
        i += 1

nporos = len(radios[len(radios)-1])

fo = open(SALIDA, "w")

for rs in radios:
    line = ""
    for r in rs:
        line += str(r) + " "
    for i in range(nporos - len(rs)):
        line += "0 "
    fo.write(line + "\n")

fo.close()
