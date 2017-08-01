CXXFLAGS = -g -O0
CXXFLAGS += -I$(SP_TOPDIR)
CXXFLAGS += -I. $(shell root-config --cflags)

LDFLAGS = -L.
LDFLAGS += $(shell root-config --libs)
LDFLAGS += -L$(SP_LIBDIR) -lSP_Base -lSP_UnfoldCore -lSP_UnfoldAlgo

CXX = g++

PROGRAMS = unfold_svd NCunfold unfold lcurve

all:	clean $(PROGRAMS)


$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
clean:	
	rm -f $(PROGRAMS) *.dYSM
