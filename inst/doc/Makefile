all: pp1.pdf
pp1.pdf: pp1.tex pp1.bbl pp1.aux
	pdflatex pp1
pp1.tex: pp1.Rnw
	echo 'Sweave("pp1.Rnw")' | R --no-save --no-restore 
pp1.bbl: pp1.tex ref.bib pp1.aux
	bibtex pp1
	latex pp1
pp1.aux: pp1.tex 
	latex pp1
clean: 
	rm -f pp1.tex pp1.aux pp1.log pp1.dvi pp1-* Rplots.ps pp1.blg
