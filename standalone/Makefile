# Makefile for DRD

INSTALLDIR = /usr/local/bin

CC = g++
LD = g++
CFLAGS = -g -Wall -c
LDFLAGS = -pthread

all: drd

util.o: util.cpp drd.h
	$(CC) $(CFLAGS) util.cpp
	
orfs.o: orfs.cpp drd.h
	$(CC) $(CFLAGS) orfs.cpp
	
align.o: align.cpp alel.h drd.h
	$(CC) $(CFLAGS) align.cpp

ribodef.o: ribodef.cpp ribodef.h
	$(CC) $(CFLAGS) ribodef.cpp

heaviestpath.o: heaviestpath.cpp alel.h drd.h
	$(CC) $(CFLAGS) heaviestpath.cpp
	
maxalign.o: maxalign.cpp alel.h drd.h
	$(CC) $(CFLAGS) maxalign.cpp
	
finder.o: finder.cpp drd.h arg_struct.h
	$(CC) $(CFLAGS) finder.cpp
	
drd.o: drd.cpp drd.h arg_struct.h
	$(CC) $(CFLAGS) drd.cpp
	
drd: util.o orfs.o align.o ribodef.o heaviestpath.o maxalign.o finder.o drd.o
	$(LD) $(LDFLAGS) -o drd util.o orfs.o align.o ribodef.o heaviestpath.o maxalign.o finder.o drd.o

clean:
	rm *.o
	
install:
	cp drd $(INSTALLDIR)
	cp Ct2B.pl $(INSTALLDIR)
	cp -r defFiles $(INSTALLDIR)