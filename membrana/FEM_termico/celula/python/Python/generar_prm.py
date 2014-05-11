# Para generar archivos .prm para AutoMesh2D sin la malla 1d ni 2d

SALIDA = "generado.prm"
ALTO_orig = 500
ANCHO_orig = 250
RADIO_orig = 50
MEMBRANA_orig = 5e-3
ESCALA = 100

ALTO = ALTO_orig * ESCALA
ANCHO = ANCHO_orig * ESCALA
RADIO = RADIO_orig * ESCALA
MEMBRANA = MEMBRANA_orig * ESCALA

fo = open(SALIDA, "w")

fo.write("""
13

""")

fo.write("1  %s %s 0.00 0\n" % (0, ALTO))
fo.write("2  %s %s 0.00 0\n" % (ANCHO, ALTO))
fo.write("3  %s %s 0.00 0\n" % (0, 0))
fo.write("4  %s %s 0.00 0\n" % (ANCHO, 0))
fo.write("5  %s %s 0.00 0\n" % (0, ALTO/2 + RADIO + MEMBRANA))
fo.write("6  %s %s 0.00 0\n" % (0, ALTO/2 + RADIO))
fo.write("7  %s %s 0.00 0\n" % (0, ALTO/2 - RADIO))
fo.write("8  %s %s 0.00 0\n" % (0, ALTO/2 - RADIO - MEMBRANA))
fo.write("9  %s %s 0.00 0\n" % (RADIO, ALTO/2))
fo.write("10 %s %s 0.00 0\n" % (RADIO + MEMBRANA, ALTO/2))
fo.write("11 %s %s 0.00 0\n" % (0, ALTO/2 + RADIO + MEMBRANA/2))
fo.write("12 %s %s 0.00 0\n" % (0, ALTO/2 - RADIO - MEMBRANA/2))
fo.write("13 %s %s 0.00 0\n" % (RADIO + MEMBRANA/2, ALTO/2))

fo.write("""
13

0
3 4 0
0
4 2 0
0
2 1 0
0
1 5 0
0
6 7 0
0
8 3 0
2
5 8 10 0
2
11 12 13 0
2
6 7 9 0
0
5 11 0
0
11 6 0
0
7 12 0
0
12 8 0

4

16711680
1

6

     3          4          7          6          1          2          
16711680
1

4

     10         8          13         -7         
16711680
1

4

     11         9          12         -8         
16711680
1

2

     5          -9         

0

0

0

0

0

0

0

0

0

0

0

0

0
0
0
0
0
0
0
0
0
0
0

""")

fo.close()
