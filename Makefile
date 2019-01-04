REB = $(HOME)/rebound/src
include $(REB)/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I$(REB) -Wl,-rpath,./ $(OPT) $(PREDEF) Alexcode.c -L. -lrebound $(LIB) -o edisk
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REB)
	@-rm -f librebound.so
	@ln -s $(REB)/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(REB) clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
