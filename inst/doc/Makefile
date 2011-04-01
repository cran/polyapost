
all: pp1.pdf clean

pp1.tex: pp1.Rnw
	$(R_HOME)/bin/R CMD Sweave pp1.Rnw

clean: 
	rm -f *.aux *.log *.dvi *.blg pp1-*.eps pp1-*.pdf Rplots.ps pp1.blg
