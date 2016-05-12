# latex Makefile
ifndef texpath
texpath=/usr/texbin/
endif
PDFLATEX=${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
SKIPERR=${texpath}pdflatex --interaction=nonstopmode
LATEX=${PDFLATEX}
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf

all: referee w51_evla_h2co

.PHONY: w51_evla_h2co
w51_evla_h2co: 
	python make.py --texpath=${texpath}
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=cores_and_outflows_compressed.pdf cores_and_outflows.pdf

.PHONY: diff
diff:
	#python2 parse_macros.py w51_evla_h2co.tex w51_evla_h2co_fordiff.tex
	#python parse_macros.py w51_evla_h2co_aug13.tex w51_evla_h2co_aug13_fordiff.tex
	#python2 parse_macros.py w51_evla_h2co_feb16_2016_submitted.tex w51_evla_h2co_feb16_2016_submitted_fordiff.tex
	latexdiff original_for_diff.tex w51_evla_h2co_fordiff.tex > diff.tex
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=diff_compressed.pdf diff.pdf

.PHONY: referee
referee:
	python make.py --referee --texpath=${texpath}
