# Build Kedlaya's OpenMP enumeration driver (root-unitary/c-script/weilpoly.c)
# against FLINT 3.  By default we use the FLINT that ships with Sage; set
# FLINT_PREFIX to use another installation.
SAGE_LOCAL ?= $(shell sage -sh -c 'echo $$SAGE_LOCAL' 2>/dev/null)
FLINT_PREFIX ?= $(SAGE_LOCAL)
CC ?= gcc
CFLAGS ?= -O3

# c/flint_compat.h supplies a function that root-unitary needs and FLINT < 3.1 lacks.
bin/weilpoly: external/root-unitary/c-script/weilpoly.c external/root-unitary/power_sums.c external/root-unitary/power_sums.h c/flint_compat.h
	mkdir -p bin
	$(CC) $(CFLAGS) -fopenmp -I$(FLINT_PREFIX)/include -include c/flint_compat.h -L$(FLINT_PREFIX)/lib \
	    -Wl,-rpath,$(FLINT_PREFIX)/lib -o $@ $< -lflint -lgmp

external/root-unitary/c-script/weilpoly.c:
	git submodule update --init

.PHONY: test clean
test: bin/weilpoly
	sage -python -m pytest -q tests

clean:
	rm -rf bin
