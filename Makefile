include Makefile.inc

all: lib/libadaptive_anderson_solver.so python test

clean:
	rm -rf lib/* test/test test/test_debug *.mod build

lib: lib/libadaptive_anderson_solver.so

lib/libadaptive_anderson_solver.so: adaptive_anderson_solver.f90
	${FC} adaptive_anderson_solver.f90 -shared -fPIC -o lib/libadaptive_anderson_solver.so ${OPTS}

lib/libadaptive_anderson_solver_debug.so: adaptive_anderson_solver.f90
	${FC} adaptive_anderson_solver.f90 -shared -fPIC -o lib/libadaptive_anderson_solver.so ${DEBUG_OPTS}

test/test: test/test.f90 lib/libadaptive_anderson_solver.so
	${FC} test/test.f90 -L./lib -ladaptive_anderson_solver $(LAPACK) -o test/test ${OPTS}

test/test_debug: test/test.f90 lib/libadaptive_anderson_solver_debug.so
	${FC} test/test.f90 -L./lib -ladaptive_anderson_solver $(LAPACK) -o test/test_debug ${DEBUG_OPTS}

run_test: test/test
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):./lib  test/test

run_test_debug: test/test_debug
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):./lib  test/test_debug

python: lib/libadaptive_anderson_solver.so
	${PYTHON} setup.py build_ext --inplace

python_test: python
	${PYTHON} test/test.py

test: run_test_debug run_test python_test

install: lib/libadaptive_anderson_solver.so
	cp lib/* ${PREFIX}/lib/

.PHONY: test lib
