DICT  = SP_BaseDict
SHLIB = libSP_Base.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

INCFLAGS  = -I.
INCFLAGS += -I$(SP_TOPDIR)
INCFLAGS += $(shell python-config --includes)
INCFLAGS += -I$(shell python -c "import numpy; print numpy.get_include()")

LDFLAGS += -L$(shell python-config --prefix)/lib/
LDFLAGS += $(shell python-config --ldflags)

include $(SP_TOPDIR)/Makefile/Makefile
