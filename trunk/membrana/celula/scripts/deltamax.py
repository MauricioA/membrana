# muestra la maxima diferencia entre dos archivos

ARCHIVO_1 = "../salida/new/10_128/transporte/pH-4-0.019999.csv"
ARCHIVO_2 = "../salida/new/10_128/transporte/pH-0-0.000000.csv"
VALORES = 2
TOL = 1e-9

f1 = open(ARCHIVO_1)
f2 = open(ARCHIVO_2)

line1 = f1.readline()
line1 = f2.readline()
line1 = f1.readline()	#ignorar primera linea
line2 = f2.readline()	#ignorar primera linea

max = 0
nlin = 1
linmax = 0

while line1 != "":
	spl1 = line1.split(',')
	spl2 = line2.split(',')

	assert abs(float(spl1[0]) - float(spl2[0]) < TOL)
	assert abs(float(spl1[1]) - float(spl2[1]) < TOL)

	for i in range(2, VALORES+2):
		diff = abs(float(spl1[i]) - float(spl2[i]))
		if diff > max:
			max = diff
			linmax = nlin
			#print "linea %s: %s" % (nlin, line1)
			#print "linea %s: %s" % (nlin, line2)

	line1 = f1.readline()
	line2 = f2.readline()
	nlin += 1

print 'max diff: %s linea %s' % (max, linmax)

f1.close()
f2.close()
