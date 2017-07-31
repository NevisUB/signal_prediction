SUBDIRS := Core Algo

.phony: all clean

all:
	@for i in $(SUBDIRS); do ( echo "" && echo "\t> Compiling Unfold/$$i <" && cd $$i && $(MAKE) ) || exit $$?; done

clean:
	@for i in $(SUBDIRS); do ( echo "" && echo "\t> Cleaning Unfold/$$i <" && cd $$i && $(MAKE) clean && rm -f $(SP_LIBDIR)/$$i.* ) || exit $$?; done

