
all: pp1.pdf clean

pp1.tex: pp1.Rnw
	$(R_HOME)/bin/R CMD Sweave pp1.Rnw

pp1.pdf: pp1.tex ref.bib
	pdflatex pp1
	bibtex pp1
	pdflatex pp1
	pdflatex pp1

clean: 
	rm -f *.aux *.log *.dvi *.blg *.bbl pp1-*.pdf Rplots.pdf
