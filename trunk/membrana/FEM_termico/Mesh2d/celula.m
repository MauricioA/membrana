n = 32;
h = 100;
r = 10;
d = 0.1

dtheta = pi/n;
theta  = (-pi/2 : dtheta : pi/2)';
memb_i = [cos(theta)*r     sin(theta)*r+h/2];
memb_o = [cos(theta)*(r+d) sin(theta)*(r+d)+h/2];
nodes  = [0 0; 50 0; 50 100; 0 100];

nodes = vertcat(nodes, memb_i, memb_o);

n_6 = size(memb_i)(1) + 4;
n_7 = size(memb_i)(1) + 4 + 1;
n_8 = size(memb_i)(1) + 4 + size(memb_o)(1);

edges_r = [
	1 2
	2 3
	3 4
	4 n_8
	n_8 n_6
	n_6 5
	5 n_7
	n_7 1
];

for i = 1 : size(memb_i)(1)-1
	edges_i(i, 1) = i + 4;
	edges_i(i, 2) = i + 4 + 1;
end

for i = 1 : size(memb_o)(1)-1
	edges_o(i, 1) = i + 4 + size(memb_i)(1);
	edges_o(i, 2) = i + 4 + size(memb_i)(1) + 1;
end

edges = vertcat(edges_r, edges_i, edges_o);

faces{1} = [1, 2, 3, 4, 8, 8 + size(edges_i)(1) + 1: 8 + size(edges_i)(1) + size(edges_o)(1)];
faces{2} = [5, 7, 9 : 8 + size(edges_i)(1) + size(edges_o)(1)];
faces{3} = [6, 8 + 1 : 8 + size(edges_i)(1)];

[p, t, g] = meshfaces(nodes, edges, faces);

f = fopen('nueva.fem', 'w');

nNodos = size(p)(1);
nElems = size(t)(1);

fprintf(f, '*COORDINATES\n');
fprintf(f, ' %d\n', nNodos);

for i = 1 : nNodos
	fprintf(f, '%d %f %f\n', i, p(i, 1), p(i, 2));
end

fprintf(f, '*ELEMENT_GROUPS\n');
fprintf(f, ' 3\n');

elemsF = [];
elemsM = [];
elemsC = [];

for i = 1 : nElems
	if   g(i) == 1
		elemsF(end+1) = i;
	elseif g(i) == 2
		elemsM(end+1) = i;
	elseif g(i) == 3
		elemsC(end+1) = i;
	end
end

fprintf(f, '1 %d Tri3\n', size(elemsF)(2));
fprintf(f, '2 %d Tri3\n', size(elemsM)(2));
fprintf(f, '3 %d Tri3\n', size(elemsC)(2));

fprintf(f, '*INCIDENCES\n');

for i = 1 : size(elemsF)(2)
	fprintf(f, '%d %d %d %d\n', i, t(elemsF(1, i), 1), t(elemsF(1, i), 2), t(elemsF(1, i), 3))
end

for i = 1 : size(elemsM)(2)
	fprintf(f, '%d %d %d %d\n', i, t(elemsM(1, i), 1), t(elemsM(1, i), 2), t(elemsM(1, i), 3))
end

for i = 1 : size(elemsC)(2)
	fprintf(f, '%d %d %d %d\n', i, t(elemsC(1, i), 1), t(elemsC(1, i), 2), t(elemsC(1, i), 3))
end

fprintf(f, '*END\n');

fclose(f);
