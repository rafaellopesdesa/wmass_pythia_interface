# Use the following on IRIX 6.5 and Linux 2.4

LIBS = `root-config --libs `
OPTCOMP = `root-config --cflags` 

#SPECIALFLAGS= --exceptions
#### TO TURN ON PDF WEIGHTS UNCOMMENT
SPECIALFLAGS= --exceptions -D__USE_PDFS__
# SPECIALFLAGS= --exceptions -D__USE_PDFS_RESBOS__
CFLAGS = $(SPECIALFLAGS) -I- -I../ -I.
LFLAGS = $(SPECIALFLAGS) -L../../lib/$(SRT_SUBDIR)/

RCXX=$(CFLAGS) $(ROOTCFLAGS)

#CC = KCC -n32 --exceptions --thread_safe -O $(OPTCOMP)
#CC = KCC +K0 --exceptions --thread_safe -O $(OPTCOMP)
CC = g++ $(RCXX) $(OPTCOMP) 

all: tupleMaker

tupleMaker: tupleMaker.o Output.o
	$(CC) tupleMaker.o Output.o $(LIBS) -o tupleMaker -lEG

Output.o: Output.cpp Output.hpp
	$(CC) -c Output.cpp -o Output.o

tupleMaker.o: tupleMaker.cpp
	$(CC) -c tupleMaker.cpp -o tupleMaker.o

clean:
	\rm -f *.o
	\rm -f *~
	\rm -fr tupleMaker
