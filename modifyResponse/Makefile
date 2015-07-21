ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'

ifndef VERBOSE
  QUIET := @
endif

CC = g++
CMSROOT = ./
INCLUDE = $(shell root-config --cflags) -I$(CMSROOT) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(INCLUDE) $(EXTRACFLAGS)
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lTMVA #-lEG -lGenVector

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs) -lEG -lGenVector -lTMVA

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)

LIB = modifyResponse.so
EXE = modifyResponse
#
# the option "-Wl,-rpath,./" adds ./ to the runtime search path for libraries
$(EXE): $(LIB)
	$(QUIET) echo "Building $@"; \
	echo "$(CC) -o $@ $(LIB) $(ROOTLIBS) -Wl,-rpath,./"; \
	$(CC) -o $@ $+ $(ROOTLIBS) -Wl,-rpath,./ 

$(LIB): $(OBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $+ $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $+ $(LINKERFLAGS)

modifyResponse.o: 	modifyResponse.cc modifyResponse.h do_convolution.h do_scale.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

do_convolution.o: 	do_convolution.cc do_convolution.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

do_scale.o: 	do_scale.cc do_scale.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o:    %.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

b: build

clean:
	rm -f \
	*.o \
	*.so \
	$(LIB) \
	$(EXE) \

endif
