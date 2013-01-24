CC = gcc
#CFLAGS = -xpg -xCC
# CFLAGS = -g -xCC -I$$HOME/include
CFLAGS = -g  -I$$HOME/progs/include
LDFLAGS= -L$$HOME/lib -lm -lkpar 
BIN = $$HOME/bin

dump_kurucz : dump_kurucz.o
	${CC} ${CFLAGS} dump_kurucz.o -o dump_kurucz
