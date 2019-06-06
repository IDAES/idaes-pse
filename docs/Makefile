# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
# DATE          = $(shell date +%Y-%m-%d_%H%M%S)
SPHINXOPTS    = -w sphinx-errors.txt -T
SPHINXBUILD   = sphinx-build
BUILDDIR      = ./build
SCRIPTDIR     = scripts

# Internal variables.
SPHINXDOCBUILDER  = doctest
SPHINXDOCDIR  	  = doctest
SPHINXDOCTEST     = $(SPHINXBUILD) -b $(SPHINXDOCBUILDER) .
ALLSPHINXOPTS     = -d $(BUILDDIR)/doctrees $(SPHINXOPTS) .
# For building API docs twice -- once indexed and once nonindexed --
# which allows us to avoid warnings in some cases and have cross-references
# work in the others. Fun!
SPHINXAPIDOC_X	  	= SPHINX_APIDOC_OPTIONS="members,ignore-module-all,noindex"  sphinx-apidoc
SPHINXAPIDOCOUT_X 	= apidoc_x
SPHINXAPIDOC_I	  	= SPHINX_APIDOC_OPTIONS="members,ignore-module-all"  sphinx-apidoc
SPHINXAPIDOCOUT_I 	= apidoc_i
SPHINXAPIDOCOUT   	= apidoc
SPHINXAPIDOCFLAGS 	= -fM 
SPHINXAPIDOCSRC 	= ../idaes
SPHINXAPIDOCEXCL  	= ../idaes/tests "../idaes/*/tests"
SPHINX_X 			= $(SPHINXAPIDOC_X) -o $(SPHINXAPIDOCOUT_X) $(SPHINXAPIDOCFLAGS) $(SPHINXAPIDOCSRC) $(SPHINXAPIDOCEXCL) 
SPHINX_I 			= $(SPHINXAPIDOC_I) -o $(SPHINXAPIDOCOUT_I) $(SPHINXAPIDOCFLAGS) $(SPHINXAPIDOCSRC) $(SPHINXAPIDOCEXCL) 
#ALLSPHINXAPIDOCOPTS_X = $(SPHINXAPIDOCFLAGS) $(SPHINXAPIDOCSRC) $(SPHINXAPIDOCEXCL)

# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(SPHINXOPTS) .

.PHONY: help
help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  apidoc     to make HTML API docs"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  applehelp  to make an Apple Help Book"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  epub3      to make an epub3"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  latexpdfja to make LaTeX files and run them through platex/dvipdfmx"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  xml        to make Docutils-native XML files"
	@echo "  pseudoxml  to make pseudoxml-XML files for display purposes"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"
	@echo "  coverage   to run coverage check of the documentation (if enabled)"
	@echo "  dummy      to check syntax errors of document sources"

clean:
	rm -f apidoc/*
	rm -f sphinx-build.log

allclean: clean
	rm -rf $(BUILDDIR)/*
	rm -f tutorials/Module_*.rst tutorials/*.html

apidoc:
	@ echo "** Build apidocs (twice)"
	$(SPHINX_I) >/dev/null
	$(SPHINX_X) >/dev/null
	@ echo "  combine indexed and non-indexed docs"
	@ mkdir -p apidoc
	@ /bin/cp $(SPHINXAPIDOCOUT_X)/* apidoc
	@for ix in idaes idaes.dmf idaes.vis ; do 	\
		/bin/cp $(SPHINXAPIDOCOUT_I)/$${ix}.rst apidoc;		\
	done
	@ /bin/rm -rf $(SPHINXAPIDOCOUT_I) $(SPHINXAPIDOCOUT_X)
	@ echo "** API docs are in ./apidoc"

all: apidoc tutorials html
	@echo "** Build finished"

html:
	SPHINX=1 $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)
	@echo "** HTML docs are in $(BUILDDIR)"

### 
# This section builds the workshop tutorials into the
# documentation (html only). It is not executed by default, so the
# process that builds for distribution should add, explicitly,
# a "tutorials" target.
###

tutorials: tutorial_setup tutorial_nb tutorial_imgs

tutorial_setup:
	@mkdir -p $(BUILDDIR)/tutorials
	@mkdir -p ./tutorials

tutorial_nb: ../examples/workshops/Module_*/*_Solution.ipynb
	@for nbfile in $^ ; do \
		echo "** Execute and convert: $$nbfile"; \
		name=$$(basename $$nbfile .ipynb); \
		title=$$(echo $$name | sed "s/_/ /g"); \
		cmd="jupyter nbconvert --to rst -y --output-dir=./tutorials $$nbfile"; \
		if test -n "$$READTHEDOCS"; then \
			echo "Generating static tutorials for ReadTheDocs..."; \
			echo "  $$cmd"; \
			$$cmd; \
		else \
			$$cmd --execute 2>/dev/null; \
		fi; \
		if [ $$? -ne 0 ]; then \
			echo "** Execute notebook failed: creating static"; \
			$${cmd}; \
		fi; \
	done

tutorial_imgs: ../examples/workshops/Module_*/*.png
	@echo "** Copy images from ../examples/workshops/* to $(BUILDDIR)/tutorials"
	@for f in $^; do \
		cp $$f $(BUILDDIR)/tutorials; \
	done

###
# End of tutorials section
###

.PHONY: dirhtml
dirhtml:
	mkdir -p $(BUILDDIR)/dirhtml
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

.PHONY: singlehtml
singlehtml:
	mkdir -p $(BUILDDIR)/singlehtml
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

.PHONY: pickle
pickle:
	mkdir -p $(BUILDDIR)/pickle
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

.PHONY: json
json:
	mkdir -p $(BUILDDIR)/json
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

.PHONY: htmlhelp
htmlhelp:
	mkdir -p $(BUILDDIR)/htmlhelp
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

.PHONY: qthelp
qthelp:
	mkdir -p $(BUILDDIR)/qthelp
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(BUILDDIR)/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/DMF.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/DMF.qhc"

.PHONY: applehelp
applehelp:
	mkdir -p $(BUILDDIR)/applehelp
	$(SPHINXBUILD) -b applehelp $(ALLSPHINXOPTS) $(BUILDDIR)/applehelp
	@echo
	@echo "Build finished. The help book is in $(BUILDDIR)/applehelp."
	@echo "N.B. You won't be able to view it unless you put it in" \
	      "~/Library/Documentation/Help or install it in your application" \
	      "bundle."

.PHONY: devhelp
devhelp:
	mkdir -p $(BUILDDIR)/devhelp
	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(BUILDDIR)/devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/DMF"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/DMF"
	@echo "# devhelp"

.PHONY: epub
epub:
	mkdir -p $(BUILDDIR)/epub
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

.PHONY: epub3
epub3:
	mkdir -p $(BUILDDIR)/epub3
	$(SPHINXBUILD) -b epub3 $(ALLSPHINXOPTS) $(BUILDDIR)/epub3
	@echo
	@echo "Build finished. The epub3 file is in $(BUILDDIR)/epub3."

.PHONY: latex
latex:
	mkdir -p $(BUILDDIR)/latex
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

.PHONY: latexpdf
latexpdf:
	mkdir -p $(BUILDDIR)/latex
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

.PHONY: latexpdfja
latexpdfja:
	mkdir -p $(BUILDDIR)/latex
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through platex and dvipdfmx..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf-ja
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

.PHONY: text
text:
	mkdir -p $(BUILDDIR)/text
	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

.PHONY: man
man:
	mkdir -p $(BUILDDIR)/man
	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(BUILDDIR)/man
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

.PHONY: texinfo
texinfo:
	mkdir -p $(BUILDDIR)/texinfo
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

.PHONY: info
info:
	mkdir -p $(BUILDDIR)/texinfo
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

.PHONY: gettext
gettext:
	mkdir -p $(BUILDDIR)/locale
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

.PHONY: changes
changes:
	mkdir -p $(BUILDDIR)/changes
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

.PHONY: linkcheck
linkcheck:
	mkdir -p $(BUILDDIR)/linkcheck
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

.PHONY: doctest
doctest:
	mkdir -p $(BUILDDIR)/$(SPHINXDOCDIR)
	$(SPHINXDOCTEST) $(BUILDDIR)/$(SPHINXDOCDIR)
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

.PHONY: coverage
coverage:
	mkdir -p $(BUILDDIR)/coverage
	$(SPHINXBUILD) -b coverage $(ALLSPHINXOPTS) $(BUILDDIR)/coverage
	@echo "Testing of coverage in the sources finished, look at the " \
	      "results in $(BUILDDIR)/coverage/python.txt."

xml:
	mkdir -p $(BUILDDIR)/xml
	$(SPHINXBUILD) -b xml $(ALLSPHINXOPTS) $(BUILDDIR)/xml
	@echo
	@echo "Build finished. The XML files are in $(BUILDDIR)/xml."

pseudoxml:
	mkdir -p $(BUILDDIR)/pseudoxml
	$(SPHINXBUILD) -b pseudoxml $(ALLSPHINXOPTS) $(BUILDDIR)/pseudoxml
	@echo
	@echo "Build finished. The pseudo-XML files are in $(BUILDDIR)/pseudoxml."

dummy:
	mkdir -p $(BUILDDIR)/dummy
	$(SPHINXBUILD) -b dummy $(ALLSPHINXOPTS) $(BUILDDIR)/dummy
	@echo
	@echo "Build finished. Dummy builder generates no files."

.PHONY: dummy pseudoxml xml clean allclean apidoc htmlnote \
		tutorial_setup tutorial_solns tutoral_probs tutorial_imgs
