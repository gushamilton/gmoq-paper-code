# Rebuild the paper's figures and supplementary tables from the aggregate results.
#
#   make figures   every figures/fig*.R    -> paper_assets/fig*.{png,pdf}, *_source.csv, *_caption.md
#   make tables    every tables/table*.R   -> paper_assets/table*.csv, *_caption.md
#   make all       both
#
# These read only the CSV files in results/ and run on a laptop in a few minutes. The
# analyses that produce results/ (from UK Biobank data, on a cluster) are in analysis/,
# one folder per paper section; see README.md.
SHELL := /bin/bash
RSCRIPT ?= Rscript

.PHONY: all figures tables

all: figures tables

figures:
	@for f in figures/fig*.R; do echo "== $$f"; $(RSCRIPT) "$$f" || exit 1; done

tables:
	@for f in tables/table*.R; do echo "== $$f"; $(RSCRIPT) "$$f" || exit 1; done
