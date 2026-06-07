# Makefile for the ABF--Fisher-Rao report.
#
#   make            -> build assets, then compile main.pdf (default)
#   make assets     -> regenerate figures/, tables/, and numbers from the CSVs
#   make check      -> verify the report numbers against the CSV-derived values
#   make pdf        -> compile main.pdf only (assumes assets exist)
#   make clean      -> remove LaTeX build artefacts (keeps main.pdf)
#   make distclean  -> also remove generated figures/tables and main.pdf
#
# Compilation prefers latexmk; if it is absent, falls back to pdflatex+bibtex.

PYTHON  ?= python3
MAIN     = main
SCRIPTS  = scripts
LATEXMK := $(shell command -v latexmk 2>/dev/null)

.PHONY: all assets check pdf clean distclean

all: assets check pdf

assets:
	$(PYTHON) $(SCRIPTS)/build_report_assets.py

check:
	$(PYTHON) $(SCRIPTS)/check_report_numbers.py

pdf:
ifneq ($(LATEXMK),)
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
else
	pdflatex -interaction=nonstopmode -halt-on-error $(MAIN).tex
	-bibtex $(MAIN)
	pdflatex -interaction=nonstopmode -halt-on-error $(MAIN).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(MAIN).tex
endif

clean:
ifneq ($(LATEXMK),)
	latexmk -c $(MAIN).tex
else
	$(RM) $(MAIN).aux $(MAIN).log $(MAIN).out $(MAIN).bbl $(MAIN).blg \
	      $(MAIN).toc $(MAIN).fls $(MAIN).fdb_latexmk $(MAIN).synctex.gz
endif

distclean: clean
	$(RM) $(MAIN).pdf
	$(RM) figures/fig0*.png figures/fig_*.png
	$(RM) tables/*.tex tables/*.csv tables/report_numbers.json
