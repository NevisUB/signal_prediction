CC=clang++
SOMAKER=clang++
CFLAGS=-O3
INCFLAGS =-I.
INCFLAGS+=$(shell root-config --cflags)
INCFLAGS += $(shell python-config --includes)
INCFLAGS += -I$(shell python -c "import numpy; print numpy.get_include()")

LDFLAGS=-L.
LDFLAGS+=$(shell root-config --libs)
LDFLAGS += -L$(shell python-config --prefix)/lib/
LDFLAGS += $(shell python-config --ldflags)

SOFLAGS=-fPIC
SOFLAGS+=-g -shared
CINTFLAGS=-p

DICT=CombinedDict
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx)) 
FMWK_HEADERS = LinkDef.h $(DICT).h 
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS=$(SOURCES:.cxx=.o)

SHLIB=libCombined.so
ROOTMAP=$(SHLIB:.so=.rootmap)

all: $(DICT).cxx $(SHLIB) 

%.o: %.cxx
	@echo '<< compiling' $< '>>'
	@echo $(CC) $(CFLAGS) $(INCFLAGS) -c $<
	@$(CC) $(CFLAGS) $(INCFLAGS) -c $<

# root class dictionary
$(DICT).cxx: $(HEADERS) LinkDef.h
	@echo '<< generating dict >>'
	@echo rootcling -f $(DICT).cxx -s $(COMBINED_LIBDIR)/$(SHLIB) -rml $(SHLIB) -rmf $(COMBINED_LIBDIR)/$(ROOTMAP) $(CINTFLAGS) -c $(INCFLAGS) $(HEADERS) LinkDef.h
	@rootcling -f $(DICT).cxx -s $(COMBINED_LIBDIR)/$(SHLIB) -rml $(SHLIB) -rmf $(COMBINED_LIBDIR)/$(ROOTMAP) $(CINTFLAGS) -c $(INCFLAGS) $(HEADERS) LinkDef.h

# rule for making shareable library
$(SHLIB): $(OBJECTS) $(DICT).cxx $(DICT).o
	@echo '<< generating $(SHLIB) >>'
	@echo $(SOMAKER) -o $(COMBINED_LIBDIR)/$(SHLIB) $(filter-out $(DICT).o, $(OBJECTS)) $(DICT).o $(SOFLAGS) $(LDFLAGS)
	@$(SOMAKER) -o $(COMBINED_LIBDIR)/$(SHLIB) $(filter-out $(DICT).o, $(OBJECTS)) $(DICT).o $(SOFLAGS) $(LDFLAGS)


clean:
	@echo rm -f *.o $(COMBINED_LIBDIR)/*.so *Dict.cxx $(COMBINED_LIBDIR)/*.rootmap $(COMBINED_LIBDIR)/*.pcm
	@rm -f $(COMBINED_BASEDIR)/*.o $(COMBINED_LIBDIR)/*.so $(COMBINED_BASEDIR)/*Dict.cxx $(COMBINED_LIBDIR)/*.rootmap $(COMBINED_LIBDIR)/*.pcm $(COMBINED_BASEDIR)/*.pcm $(COMBINED_BASEDIR)/*Dict.cxx_tmp
