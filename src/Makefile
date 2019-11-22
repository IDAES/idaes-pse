# Make it slightly easier to compile external functions and build docs
#   The external functions all require the ASL, so set the ASL_BUILD environment
#   variable to point to the ASL build and this should be good to go.
include Makefile.in

ALL: iapws95 cubic functions dist

.PHONY: iapws95
iapws95:
	$(MAKE) -C ./iapws95

iapws95_clean:
	$(MAKE) -C ./iapws95 clean

.PHONY: cubic
cubic:
	$(MAKE) -C ./cubic

cubic_clean:
	$(MAKE) -C ./cubic clean

.PHONY: functions
functions:
	$(MAKE) -C ./functions

functions_clean:
	$(MAKE) -C ./functions clean

.PHONY: dist
dist:
	$(RM) -rf dist
	$(MKDIR) dist/
	$(CP) iapws95/*.so dist/
	$(CP) cubic/*.so dist/
	$(CP) functions/*.so dist/

.PHONY: install
install:
	idaes lib-directory --create
	$(CP) dist/* `idaes lib-directory`

dist_clean:
	$(RM) -rf dist

clean: iapws95_clean cubic_clean functions_clean dist_clean dist_clean
