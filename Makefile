include Makefile.inc

all: fortran python

fortran: lib fortran_test

clean:
	rm -rf dist *.mod lib/* test/test test/test_debug src/fortran/*.mod build example/*/example *.pyc src/*.cpython-* adaptive_anderson_solver.c src/*.egg-info

lib: lib/libadaptive_anderson_solver.so

lib/libadaptive_anderson_solver.so: src/fortran/adaptive_anderson_solver.f90
	${FC} src/fortran/adaptive_anderson_solver.f90 -shared -fPIC -o lib/libadaptive_anderson_solver.so ${OPTS}

lib/libadaptive_anderson_solver_debug.so: src/fortran/adaptive_anderson_solver.f90
	${FC} src/fortran/adaptive_anderson_solver.f90 -shared -fPIC -o lib/libadaptive_anderson_solver.so ${DEBUG_OPTS}

test/test: test/test.f90 lib/libadaptive_anderson_solver.so
	${FC} test/test.f90 -L./lib -Isrc/fortran -ladaptive_anderson_solver $(LAPACK) -o test/test ${OPTS}

test/test_debug: test/test.f90 lib/libadaptive_anderson_solver_debug.so
	${FC} test/test.f90 -L./lib -Isrc/fortran -ladaptive_anderson_solver $(LAPACK) -o test/test_debug ${DEBUG_OPTS}

run_test: test/test
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):./lib  test/test > /dev/null

run_test_debug: test/test_debug
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):./lib  test/test_debug > /dev/null

python: lib/libadaptive_anderson_solver.so
	pip install -r requirements.txt
	${PYTHON} setup.py build_ext --inplace

python_editable_install:
	pip install --user -e .

python_test: python
	${PYTHON} test/test.py > /dev/null

test: fortran_test python_test

fortran_test: lib test/test test/test_debug run_test_debug run_test python_test

install: fortran_install python_install

fortran_install: lib/libadaptive_anderson_solver.so
	cp lib/* ${PREFIX}/lib/

python_install: lib/libadaptive_anderson_solver.so
	pip install .

.PHONY: test lib python_install fortran_install fortran_test python_test python fortran clean python_edit_install
