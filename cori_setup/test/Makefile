OSKIDIR =  /global/homes/j/jjnjiang/project/poski-v1.0.0/build-1.0.0/build_oski
OSKIINC = $(OSKIDIR)/include
OSKILIB = $(OSKIDIR)/lib/oski

#Location of pOSKI library
POSKIDIR =  /global/homes/j/jjnjiang/project/poski-v1.0.0/build-1.0.0
POSKILIB = $(POSKIDIR)/lib
POSKIINC = $(POSKIDIR)/include/poski

#OSKI & pOSKI link flags
OSKILIB_SHARED = -I$(OSKIINC) -Wl,-rpath -Wl,$(OSKILIB) -L$(OSKILIB) `cat $(OSKILIB)/site-modules-shared.txt` -loski
POSKILIB_SHARED = -I$(POSKIINC) -Wl,-rpath -Wl,$(POSKILIB) -L$(POSKILIB) -lposki
LDFLAGS_SHARED = $(OSKILIB_SHARED) $(POSKILIB_SHARED) -lm

CC = icc
CFLAGS = -g -O3 -pthread -qopenmp

SRC = test_p
all:$(SRC)-shared
$(SRC)-shared: $(SRC).o
	$(CC) $(CFLAGS) -o $(SRC) $(SRC).o $(LDFLAGS_SHARED)
.c.o:
	$(CC) -c $(CFLAGS) $(LDFLAGS_SHARED) -o $(SRC).o $(SRC).c
clean:
	rm -rf $(SRC)-shared $(SRC).o core*~
