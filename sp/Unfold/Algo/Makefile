DICT  = SP_UnfoldAlgoDict
SHLIB = libSP_UnfoldAlgo.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

INCFLAGS  = -I.
INCFLAGS += -I$(SP_TOPDIR)

LDFLAGS += -L$(SP_LIBDIR)
LDFLAGS += -lSP_UnfoldCore

include $(SP_TOPDIR)/Makefile/Makefile
