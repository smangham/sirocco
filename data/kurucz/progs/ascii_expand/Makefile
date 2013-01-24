CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

ascii_expand: ascii_expand.o
		gcc ascii_expand.o  $(LDFLAGS) -o ascii_expand

