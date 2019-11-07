# Make it slightly easier to compile external functions and build docs
#   The external functions all require the ASL, so set the ASL_BUILD environment
#   variable to point to the ASL build and this should be good to go.
include Makefile.in

ALL: iapws95 cubic_eos functions dist

.PHONY: iapws95
iapws95:
	$(MAKE) -C ./iapws95_lib

iapws95_clean:
	$(MAKE) -C ./iapws95_lib clean

.PHONY: cubic_eos
cubic_eos:
	$(MAKE) -C ./cubic_eos

cubic_eos_clean:
	$(MAKE) -C ./cubic_eos clean

.PHONY: functions
functions:
	$(MAKE) -C ./functions

functions_clean:
	$(MAKE) -C ./functions clean

.PHONY: dist
dist:
	$(RM) -rf dist
	$(MKDIR) dist/
	$(CP) iapws95_lib/*.so dist/
	$(CP) cubic_eos/*.so dist/
	$(CP) functions/*.so dist/

dist_clean:
	$(RM) -rf dist

clean: iapws95_clean cubic_eos_clean functions_clean dist_clean dist_clean
