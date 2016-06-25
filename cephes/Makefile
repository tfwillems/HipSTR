# Makefile for probability integrals.
# Be sure to set the type of computer and endianness in mconf.h.

CC = gcc
CFLAGS = -O2 -Wall
INCS = mconf.h

OBJS = bdtr.o btdtr.o chdtr.o drand.o expx2.o fdtr.o gamma.o gdtr.o \
igam.o igami.o incbet.o incbi.o mtherr.o nbdtr.o ndtr.o ndtri.o pdtr.o \
stdtr.o unity.o polevl.o const.o isnan.o

.PHONY: all
all: libprob.a

.PHONY: clean
clean:
	rm -f *.o libprob.a

libprob.a: $(OBJS) $(INCS)
	ar rv libprob.a $(OBJS)
	ranlib libprob.a
