CC = gcc
CFLAGS = -g -Wall -I$$PYTHON/include
LDFLAGS= -L$$PYTHON/lib
INCLUDE=$$PYTHON/include
BIN = $$PYTHON/bin



pydist: pydist.o  
	$(CC) $(CFLAGS) $(LDFLAGS) pydist.o \
	-lkpar -lm -o  pydist  
	mv $@ $(BIN)/pydist


pyfit3: pyfit3.o time.o 
	$(CC) $(CFLAGS) $(LDFLAGS) pyfit3.o time.o \
	-lkpar -lm -o  pyfit3  
	mv $@ $(BIN)/pyfit3

# Make templates for pyfit3
pyfit3_proto: 
	cp pyfit3_templates pyfit3_templates.old
	cproto -I$(INCLUDE) pyfit3.c > foo.h
	mv foo.h pyfit3_templates.h

py_compare_obj = py_compare.o  get_python_spec.o

py_compare_src = py_compare.c  get_python_spec.c

py_compare: $(py_compare_obj)
	$(CC) $(CFLAGS) $(LDFLAGS)  $(py_compare_obj) \
	-lkpar -lm -o   py_compare 
	mv $@ $(BIN)/py_compare


prototypes: 
	cp templates.h templates.h.old
	cproto -I$(INCLUDE)  $(py_compare_src) > foo.h
	mv foo.h templates.h
	
# Old
pyfit2: pyfit2.o  
	$(CC) $(CFLAGS) $(LDFLAGS) pyfit2.o \
	-lkpar -lm -o  pyfit2  
	mv $@ $(BIN)/pyfit2


pyfit: pyfit.o  
	$(CC) $(CFLAGS) $(LDFLAGS) pyfit.o \
	-lkpar -lm -o  pyfit  
	mv $@ $(BIN)/pyfit


