CC = gcc
CFLAGS = -g -Wall -I$$PYTHON/include
LDFLAGS= -L$$PYTHON/lib
INCLUDE=$$PYTHON/include
BIN = $$PYTHON/bin



py_compare_obj = py_compare.o  get_python_spec.o

py_compare_src = py_compare.c  get_python_spec.c

py_compare: $(py_compare_obj)
	$(CC) $(CFLAGS) $(LDFLAGS)  $(py_compare_obj) \
	-lkpar -lm -o   py_compare 


prototypes: 
	cp templates.h templates.h.old
	cproto -I$(INCLUDE)  $(py_compare_src) > foo.h
	mv foo.h templates.h
	
tread: tread.o  
	$(CC) $(CFLAGS) $(LDFLAGS) tread.o \
	-lkpar -lm -o  tread  

