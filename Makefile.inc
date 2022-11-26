FC=gfortran
LAPACK=-llapack -lblas -lpthread
PYTHON=python3
PREFIX=/usr/local
COMMON_OPTS=
DEBUG_OPTS=${COMMON_OPTS} -g -fbounds-check
OPTS=${COMMON_OPTS} -O3
