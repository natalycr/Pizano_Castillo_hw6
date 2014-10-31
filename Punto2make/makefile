FOO = $(kinetic)
FEE = $(alpha)



trayectoria_E_alpha.dat: particle_in_field.c
	cc particle_in_field.c
	./a.out FOO FEE
	python plotxyz.py plotxy.py #trayectoria_E_alpha.dat