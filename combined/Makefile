CC=clang++
SOMAKER=clang++
CFLAGS=-g
INCFLAGS =-I.
INCFLAGS+=$(shell root-config --cflags)
LDFLAGS=-L.
LDFLAGS+=$(shell root-config --libs)

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
	@echo rootcling -f $(DICT).cxx -s $(SHLIB) -rml $(SHLIB) -rmf $(ROOTMAP) $(CINTFLAGS) -c $(INCFLAGS) $(HEADERS) LinkDef.h
	@rootcling -f $(DICT).cxx -s $(SHLIB) -rml $(SHLIB) -rmf $(ROOTMAP) $(CINTFLAGS) -c $(INCFLAGS) $(HEADERS) LinkDef.h

# rule for making shareable library
$(SHLIB): $(OBJECTS) $(DICT).cxx $(DICT).o
	@echo '<< generating $(SHLIB) >>'
	@echo $(SOMAKER) -o $(SHLIB) $(filter-out $(DICT).o, $(OBJECTS)) $(DICT).o $(SOFLAGS) $(LDFLAGS)
	@$(SOMAKER) -o $(SHLIB) $(filter-out $(DICT).o, $(OBJECTS)) $(DICT).o $(SOFLAGS) $(LDFLAGS)


clean:
	@echo rm -f *.o *.so *Dict.cxx *.rootmap *.pcm
	@rm -f *.o *.so *Dict.cxx *.rootmap *.pcm
