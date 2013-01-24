C = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

pyfit4: pyfit4.o get_models.o get_python_models.o ascii_read.o \
	get_data.o model_arith.o time.o
	gcc  pyfit4.o get_models.o  get_python_models.o  ascii_read.o \
		get_data.o  model_arith.o time.o $(LDFLAGS) -o pyfit4
	mv pyfit4 ../../bin


t_get_models: t_get_models.o get_models.o get_python_models.o ascii_read.o \
	get_data.o model_arith.o 
	gcc  t_get_models.o get_models.o  get_python_models.o  ascii_read.o \
		get_data.o  model_arith.o $(LDFLAGS) -o t_get_models

xtest: xtest.o
	gcc xtest.o  $(LDFLAGS) -o xtest


t_ascii: t_ascii.o ascii_read.o 
	gcc t_ascii.o ascii_read.o $(LDFLAGS) -o t_ascii


prototypes: 
	cproto -I$$HOME/include  get_data.c get_models.c get_python_models.c ascii_read.c model_arith.c >model_prototypes.h  2>proto.debug
