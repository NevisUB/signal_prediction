SUBDIRS := Core Algo

.phony: all clean

all:
	@for i in $(SUBDIRS); do ( echo "" && echo "\t>Compiling $$i<" && cd $$i && $(MAKE) ) || exit $$?; done

clean:
	@for i in $(SUBDIRS); do ( echo "" && echo "\t>Cleaning $$i<" && cd $$i && $(MAKE) clean && rm -f $(SP_LIBDIR)/$$i.* ) || exit $$?; done

