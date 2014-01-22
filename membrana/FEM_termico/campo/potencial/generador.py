import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument('file' )
args = parser.parse_args()

f = open(args.file)
lines = f.readlines()
f.close()

v = []
t = []

for line in lines:
	if re.match("\*ELEMENT\_GROUPS", line):
		break
	if re.match("\d+\s\d+\.?\d*\s100\.?0*", line):
		v.append(re.match("\d+", line).group())
	elif re.match("\d+\s\d+\.?\d*\s0\.?0*\s", line):
		t.append(re.match("\d+", line).group())

f = open("input.in", "w")

f.write("Problema: celula\n")

f.write("\n")
f.write("data_entrada\n")
f.write("archivo1: " + args.file + "\n")
f.write("opcion: 1\n")
f.write("archivo3: sistema.dat\n")
f.write("""
sigint:  	0.00000015	#condutividad de la zona intraelular [S/um]
sigext:  	0.0000002    	#condutividad de la zona extra
sigmem: 	0.0000000000005		#condutividad de la membrana
permit :  	1.0       	#permitividad de la membrana
Potencial: 	1.0     	#voltaje 
Frecuencia:	0.0  		#freq del campo en MGhz
end_data

""")

f.write("dirichV: " + str(len(v)) + "\n")
for p in v:
	f.write(str(p) + "\n")

f.write("\n")

f.write("dirichT: " + str(len(t)) + "\n")
for p in t:
	f.write(str(p) + "\n")

f.write("\n")
f.write("end_dataCC")
f.write("\n")

f.close()
