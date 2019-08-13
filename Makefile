# Make it slightly easier to compile external functions and build docs
#   The external functions all require the ASL, so set the ASL_BUILD environment
#   variable to point to the ASL build and this should be good to go.

ALL: iapws95 cubic_eos functions

iapws95:
	$(MAKE) -C ./idaes/property_models/iapws95_lib

iapws95_clean:
	$(MAKE) -C ./idaes/property_models/iapws95_lib clean

cubic_eos:
	$(MAKE) -C ./idaes/property_models/cubic_eos

cubic_eos_clean:
	$(MAKE) -C ./idaes/property_models/cubic_eos clean

functions:
	$(MAKE) -C ./idaes/functions

functions_clean:
	$(MAKE) -C ./idaes/functions clean

clean: iapws95_clean cubic_eos_clean functions_clean dist_clean

docs: docs_html

.PHONY: dist
dist:
	python setup.py bdist_wheel && /bin/ls -l dist/

alldocs:
	$(MAKE) -C ./docs all

docs_html:
	$(MAKE) -C ./docs apidoc html

docs_clean:
	$(MAKE) -C ./docs allclean

dist_clean:
	/bin/rm -rf dist
	/bin/rm -rf build