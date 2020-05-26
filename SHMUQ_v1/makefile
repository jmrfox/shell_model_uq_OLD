all : anamil.x standardorderint.x

standardorderint.x : standardorderint.v2.f90
	gfortran -o standardorderint.x standardorderint.v2.f90

anamil.x : anamil.v3.f90 dlapackeig.f
	gfortran -o anamil.x anamil.v3.f90 dlapackeig.f

clean: 
	rm anamil.x standardorderint.x *.o
