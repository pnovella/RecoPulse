# specific names for this package

GATEMAKEDIR = $(GATE_DIR)/Makefile
GATELIBDIR = $(GATE_DIR)/lib

DICT  = RecoPulseCint
SHLIB = libRecoPulse.so
SOURCES = $(wildcard *.cc)
#HEADERS = $(wildcard *.h)
HEADERS = $(filter-out RecoPulseCint.h,$(wildcard *.h))
#OBJECTS = $(SOURCES:.cc=.o)
OBJECTS = $(filter-out RunRecoPulse.o, $(SOURCES:.cc=.o))


# include options for this package
INCFLAGS = -I.
INCFLAGS += -I$(GATE_DIR)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(GATEMAKEDIR)/Makefile.${OSNAME}

# set compiler options for ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += '-fPIC'


# call the common GNUmakefile
include $(GATEMAKEDIR)/GNUmakefile.GATE

#all: bin

### binary compilation ###

LIBS += $(shell root-config --libs)  -lCore -lRIO -lHist
LIBS += -L../../lib  -lGATE -lGATEIO -lGATEUtils -lGateModule

bin: lib/$(SHLIB) RunRecoPulse.o
	@echo "<< compiling RunRecoPulse >>"
	@$(CXX) -g $(CXXFLAGS) -o ./bin/RunRecoPulse $^ $(LIBS) 

.PHONY: clean cleanapps

clean:  cleanapps

cleanapps:
	@rm -f ./bin/RunRecoPulse

