OSKIDIR =  /global/homes/j/jjnjiang/project/poski-v1.0.0/build-1.0.0/build_oski
OSKIINC = $(OSKIDIR)/include
OSKILIB = $(OSKIDIR)/lib/oski

#Location of pOSKI library
POSKIDIR =  /global/homes/j/jjnjiang/project/poski-v1.0.0/build-1.0.0
POSKILIB = $(POSKIDIR)/lib
POSKIINC = $(POSKIDIR)/include/poski

PETSC_DIR=/global/u1/j/jjnjiang/Pet/petsc
PETSC_ARCH=arch-cray-xc40-knl-opt

#OSKI & pOSKI link flags
OSKILIB_SHARED = -I$(OSKIINC) -Wl,-rpath -Wl,$(OSKILIB) -L$(OSKILIB) `cat $(OSKILIB)/site-modules-shared.txt` -loski
POSKILIB_SHARED = -I$(POSKIINC) -Wl,-rpath -Wl,$(POSKILIB) -L$(POSKILIB) -lposki
LDFLAGS_SHARED = $(OSKILIB_SHARED) $(POSKILIB_SHARED) -lm

CC = mpiicc
CFLAGS = -g -O3 -pthread -qopenmp

SRC = test_p
all:$(SRC)-shared
$(SRC)-shared: $(SRC).o
	$(CC) $(CFLAGS) -o $(SRC) $(SRC).o $(LDFLAGS_SHARED) 
.c.o:
	$(CC) -c $(CFLAGS) -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include $(LDFLAGS_SHARED) -o $(SRC).o $(SRC).c
	#$(CC) -c $(CFLAGS) -I/opt/cray/pe/petsc/3.8.4.0/real/INTEL/16.0/haswell/include/ $(LDFLAGS_SHARED) -o $(SRC).o $(SRC).c
clean:
	rm -rf $(SRC)-shared $(SRC).o core*~

