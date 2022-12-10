include ../../Makefile.inc

all: example run

clean:
	rm -f example

example: example.f90 ../../lib/libadaptive_anderson_solver.so
	${FC} example.f90 -L../../lib -I../.. -ladaptive_anderson_solver $(LAPACK) -o example -g -lpthread -fbounds-check

run: example
	LD_LIBRARY_PATH=./../../lib:${LD_LIBRARY_PATH} ./example
