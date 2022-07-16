echo ... build DLL from fortran files
gfortran -c -O3 -w -fimplicit-none f90\arrayprint.f90
gfortran -c -O3 -w -fimplicit-none f90\orthpoly.f90
gfortran -c -O3 -w -fimplicit-none f90\quad.f90
gfortran -shared -o quadlib.so quad.o arrayprint.o orthpoly.o -fPIC -Xlinker -Map=lib.map
del *.o
del *.mod
