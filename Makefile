ALL: iapws95

# can add contrib and core ... if needed

# for now iapws95 is the only one, but soon cubicEOS will be here too

iapws95:
	$(MAKE) -C ./idaes/property_models/iapws95

iapws95_clean:
	$(MAKE) -C ./idaes/property_models/iapws95 clean

clean: iapws95_clean


# Couldn't help throwing this in
docs: docs_html

docs_html:
	$(MAKE) -C ./docs html

docs_clean:
	$(MAKE) -C ./docs clean
