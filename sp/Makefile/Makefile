CC=g++
SOMAKER=g++
CFLAGS=-g -fPIC -std=c++11 -Wall -O0

INCFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --libs)

SOFLAGS = -fPIC
SOFLAGS += -g -shared -O0
CINTFLAGS = -p

ROOTMAP = $(SHLIB:.so=.rootmap)
PCM = $(SHLIB:.so=_rdict.pcm)

all: $(DICT).cxx $(SHLIB) 

%.o: %.cxx
	@echo '<< compiling' $< '>>'
	@$(CC) $(CFLAGS) $(INCFLAGS) -c $<

# root class dictionary
$(DICT).cxx: $(HEADERS) LinkDef.h
	@echo '<< generating $(DICT) >>'
	@rootcling -f $(DICT).cxx -s $(SP_LIBDIR)/$(SHLIB) -rml $(SHLIB) -rmf $(SP_LIBDIR)/$(ROOTMAP) $(CINTFLAGS) -c $(INCFLAGS) $(HEADERS) LinkDef.h

# rule for making shareable library
$(SHLIB): $(OBJECTS) $(DICT).cxx $(DICT).o
	@echo '<< generating $(SHLIB) >>'
	@$(SOMAKER) -o $(SP_LIBDIR)/$(SHLIB) $(filter-out $(DICT).o, $(OBJECTS)) $(DICT).o $(SOFLAGS) $(LDFLAGS)


clean:
	@echo '<< cleaning >>'
	@rm -f *.o *Dict.cxx *Dict.cxx_tmp *.pcm $(SP_LIBDIR)/$(PCM) $(SP_LIBDIR)/$(ROOTMAP) $(SP_LIBDIR)/$(SHLIB)
