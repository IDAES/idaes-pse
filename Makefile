ALL: iapws95 cubic_eos

# can add contrib and core ... if needed

# for now iapws95 is the only one, but soon cubicEOS will be here too

iapws95:
	$(MAKE) -C ./idaes/property_models/iapws95

iapws95_clean:
	$(MAKE) -C ./idaes/property_models/iapws95 clean

cubic_eos:
	$(MAKE) -C ./idaes/property_models/cubic_eos

cubic_eos_clean:
	$(MAKE) -C ./idaes/property_models/cubic_eos clean

clean: iapws95_clean cubic_eos_clean


# Couldn't help throwing this in
docs: docs_html

docs_html:
	$(MAKE) -C ./docs_rst html

docs_clean:
	$(MAKE) -C ./docs_rst clean
