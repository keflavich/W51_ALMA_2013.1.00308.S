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

all: referee w51_alma_00308

.PHONY: w51_alma_00308
w51_alma_00308: 
	echo "texpath: ${texpath}"
	python make.py --texpath=${texpath}
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=w51_alma_00308_compressed.pdf w51_alma_00308.pdf

.PHONY: diff
diff:
	#python2 parse_macros.py w51_alma_00308.tex w51_alma_00308_fordiff.tex
	#python parse_macros.py w51_alma_00308_aug13.tex w51_alma_00308_aug13_fordiff.tex
	#python2 parse_macros.py w51_alma_00308_feb16_2016_submitted.tex w51_alma_00308_feb16_2016_submitted_fordiff.tex
	latexdiff original_for_diff.tex w51_alma_00308_fordiff.tex > diff.tex
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=diff_compressed.pdf diff.pdf

.PHONY: referee
referee:
	echo "texpath: ${texpath}"
	python make.py --referee --texpath=${texpath}
