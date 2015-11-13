#############################################################
# TO BE CHANGED BY EACH USER TO POINT TO include/ AND lib/ 
# DIRS HOLDING CFITSIO *.h AND libcfitsio IF THEY ARE NOT IN 
# THE STANDARD PLACES
# 

CFITSIOINCDIR =  ../../cfitsio/include
LIBDIR        =  ../../cfitsio/lib

#
#
#############################################################
# COMPILATION OPTIONS BELOW
# 

# another good memory checker is valgrind : http://valgrind.kde.org/index.html
# valgrind --tool=memcheck hotpants

# for memory checking with libefence
# LIBS  = -L$(LIBDIR) -lm -lcfitsio -lefence

# for profiling with gprof
# COPTS = -pg -fprofile-arcs -funroll-loops -O3 -ansi -pedantic-errors -Wall -I$(CFITSIOINCDIR) 

# for gdbugging
#COPTS = -g3 -funroll-loops -O3 -ansi -pedantic-errors -Wall -I$(CFITSIOINCDIR) 

# standard usage
# recently added -std=c99 after a bug report
COPTS = -funroll-loops -O3 -ansi -std=c99 -pedantic-errors -Wall -I$(CFITSIOINCDIR) -D_GNU_SOURCE
LIBS  = -L$(LIBDIR) -lm -lcfitsio

# compiler
CC    = gcc 

#
#
############################################################# 
# BELOW SHOULD BE OK, UNLESS YOU WANT TO COPY THE EXECUTABLES
# SOMEPLACE AFTER THEY ARE BUILT eg. hotpants
#

STDH  = functions.h globals.h defaults.h
ALL   = main.o vargs.o alard.o functions.o 

all:	hotpants extractkern maskim

hotpants: $(ALL)
	$(CC) $(ALL) -o hotpants $(LIBS) $(COPTS)
#	cp hotpants ../../bin/$(ARCH)

main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

alard.o: $(STDH) alard.c
	$(CC) $(COPTS)  -c alard.c

functions.o: $(STDH) functions.c
	$(CC) $(COPTS)  -c functions.c

vargs.o: $(STDH) vargs.c
	$(CC) $(COPTS)  -c vargs.c

extractkern : extractkern.o 
	$(CC) extractkern.o -o extractkern $(LIBS) $(COPTS)

extractkern.o : $(STDH) extractkern.c
	$(CC) $(COPTS)  -c extractkern.c

maskim : maskim.o
	$(CC) maskim.o -o maskim $(LIBS) $(COPTS)

maskim.o: $(STDH) maskim.c
	$(CC) $(COPTS)  -c maskim.c

clean :
	rm -f *.o
	rm -f *~ .*~
	rm -f hotpants
	rm -f extractkern
	rm -f maskim
