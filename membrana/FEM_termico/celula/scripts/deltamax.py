# muestra la maxima diferencia entre dos archivos

ARCHIVO_1 = "../salida/new/25_128_cls_t1k/transporte/concent-3-0.015000.csv"
ARCHIVO_2 = "../salida/new/25_128_cls_t100/transporte/concent-3-0.015000.csv"
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

while line1 != "":
	spl1 = line1.split(',')
	spl2 = line2.split(',')

	assert abs(float(spl1[0]) - float(spl2[0]) < TOL)
	assert abs(float(spl1[1]) - float(spl2[1]) < TOL)

	for i in range(2, VALORES+2):
		diff = abs(float(spl1[i]) - float(spl2[i]))
		if diff > max:
			max = diff
			print "linea %s: %s" % (nlin, line1)
			print "linea %s: %s" % (nlin, line2)

	line1 = f1.readline()
	line2 = f2.readline()
	nlin += 1

print 'max diff: %s' % max

f1.close()
f2.close()
