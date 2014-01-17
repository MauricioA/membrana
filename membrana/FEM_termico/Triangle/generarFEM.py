import argparse, re, math

d = 5e-3
toler = 1e-3

parser = argparse.ArgumentParser()
parser.add_argument('node')
parser.add_argument('ele')
args = parser.parse_args()

node_f = open(args.node)
node_lines = node_f.readlines()
node_f.close()

ele_f = open(args.ele)
ele_lines = ele_f.readlines()
ele_f.close()

line = node_lines[0]
reg = re.match('\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+', line)
n_nodes = int(reg.groups()[0])

nodes = []

for i in range(n_nodes):
	line = node_lines[i+1]
	reg = re.match('\s*(\d+\.?\d*\e?\+?\-?\d*)\s+(\d+\.?\d*\e?\+?\-?\d*)\s+(\d+\.?\d*\e?\+?\-?\d*)\s+', line)
	if not reg: 
		print i 
		print line
		exit()
	nodes.append((
		float(reg.groups()[1]),
		float(reg.groups()[2])
	))

line = ele_lines[0]
reg = re.match('\s*(\d+)', line)
n_tri = int(reg.groups()[0])

elems = ([], [], [])

for i in range(n_tri):
	line = ele_lines[i+1]
	reg = re.match('\s*\d+\s+(\d+)\s+(\d+)\s+(\d+)', line)
	n1 = int(reg.groups()[0])
	n2 = int(reg.groups()[1])
	n3 = int(reg.groups()[2])
	nodos = [nodes[n1-1], nodes[n2-1], nodes[n3-1]]
	for nodo in nodos:
		dist = math.sqrt((nodo[0]-0)**2 + (nodo[1]-50)**2)
		if dist > (10 + d + toler):
			grupo = 0
			break
		elif dist < (10 - toler):
			grupo = 2
			break
		else:
			grupo = 1
	elems[grupo].append((n1, n2, n3))
	
print "*COORDINATES"
print n_nodes

i = 1
for node in nodes:
	print "%d %f %f" % (i, node[0], node[1])
	i += 1
	
print "*ELEMENT_GROUPS"
print "3"

for i in range(3):
	print "%d %d Tri3" % (i+1, len(elems[i]))

print "*INCIDENCES"

for g in range(3):
	i = 1
	for elem in elems[g]:
		print "%d %d %d %d" % (i, elem[0], elem[1], elem[2])
		i += 1
