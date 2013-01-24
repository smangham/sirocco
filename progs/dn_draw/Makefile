C = gcc
#CFLAGS = -xpg -xCC
# CFLAGS = -g -xCC -I$$HOME/include
CFLAGS = -g  -I$$HOME//include
LDFLAGS= -L$$HOME/lib -lm -lkpar -lcfitsio
BIN = $$HOME/bin


dn: dn.o vector.o
	${CC} ${CFLAGS} dn.o vector.o \
	${LDFLAGS} -o dn

tvector : tvector.o vector.o
	${CC} ${CFLAGS} tvector.o vector.o \
	${LDFLAGS} -o tvector

