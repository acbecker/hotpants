CFITSIOINCDIR = /usr/local/include

#use only with no debugging
COPTS = -funroll-loops -O3 -ansi -Wall -I$(CFITSIOINCDIR) -I/usr/include/malloc
LIBS  =  -lm -lcfitsio  

CC    = gcc 

STDH  = functions.h globals.h defaults.h
ALL   = main.o vargs.o alard.o functions.o 
SWIG  = alard.o functions.o 
ALLT  = main_test.o vargs.o alard.o functions.o 

all:	hotpants 

hotpants: $(ALL)
	$(CC) $(ALL) -o hotpants $(LIBS) $(COPTS)

hotpants_test: $(ALLT)
	$(CC) $(ALLT) -o hotpants_test $(LIBS) $(COPTS)

main_test.o: $(STDH) main_test.c
	$(CC) $(COPTS)  -c main_test.c

main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

alard.o: $(STDH) alard.c
	$(CC) $(COPTS)  -c alard.c

functions.o: $(STDH) functions.c
	$(CC) $(COPTS)  -c functions.c

vargs.o: $(STDH) vargs.c
	$(CC) $(COPTS)  -c vargs.c

NCOPTS = -funroll-loops -O3 -I$(CFITSIOINCDIR) -I/usr/include/malloc

extractkern : extractkern.c
	$(CC) $(NCOPTS) extractkern.c -o extractkern $(LIBS) 

maskim : maskim.c
	$(CC) $(NCOPTS) maskim.c -o maskim $(LIBS)

clean :
	rm -f *.o
	rm -f *~ .*~
	rm -f hotpants
	rm -f extractkern
	rm -f maskim
