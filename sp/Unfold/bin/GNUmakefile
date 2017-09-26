CXXFLAGS = -g -O0
CXXFLAGS += -I$(SP_TOPDIR)
CXXFLAGS += -I. $(shell root-config --cflags)

LDFLAGS = -L.
LDFLAGS += $(shell root-config --libs)
LDFLAGS += -L$(SP_LIBDIR) -lSP_Base -lSP_UnfoldCore -lSP_UnfoldAlgo -lSP_Combined

CXX = g++

PROGRAMS = multi_unfold stacked_backgrounds

all:	clean $(PROGRAMS)


$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
clean:	
	rm -f $(PROGRAMS) *.dYSM
