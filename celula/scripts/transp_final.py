#archivo de transporte viejo -> concentracion en un instante

final_time = 20e-3
initial_pos = 100e6
input_folder = '../salida/25-64/160KVm/'
input_file = input_folder + 'transporte.csv'
output_file = input_folder + 'transp_final.csv'
epsilon = 1e-10

fo = open(output_file, 'w')

with open(input_file) as fi:
	fi.seek(initial_pos)
	fi.readline()
	fo.write('x, y, h, oh, na, cl\n')

	for line in fi:
		time, x, y, h, oh, na, cl = line.split(',')

		diff = float(time) - final_time
		if diff > -epsilon and diff < epsilon: 
			fo.write('%s, %s, %s, %s, %s, %s\n' % (x, y, h, oh, na, cl))
		elif diff > 1e-12: 
			break

fo.close()
