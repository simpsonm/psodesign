# LaTeX Makefile

## Define the file name (.tex)

DOC=paper
APP=appendix
DOCB=paperblind
APPB=appendixblind
FULL=paperfull
FULLB=paperblindfull

all:	$(DOC).pdf $(DOCB).pdf $(FULL).pdf $(FULLB).pdf

.PHONY: clean

clean:
	\rm *.aux *.blg *.out *.bbl *.log $(DOC).pdf $(DOCB).pdf \
 $(APP).pdf $(APPB).pdf $(FULL).pdf $(FULLB).pdf $(DOCB).tex $(APPB).tex

$(DOC).pdf: 	$(DOC).tex $(APP).tex
	pdflatex "\def\blind{0} \input{$(DOC).tex}"
	pdflatex "\def\blind{0} \input{$(APP).tex}"
	pdflatex "\def\blind{0} \input{$(DOC).tex}"
	pdflatex "\def\blind{0} \input{$(APP).tex}"
	bibtex $(DOC)
	bibtex $(APP)
	pdflatex "\def\blind{0} \input{$(DOC).tex}"
	pdflatex "\def\blind{0} \input{$(APP).tex}"
	pdflatex "\def\blind{0} \input{$(DOC).tex}"
	pdflatex "\def\blind{0} \input{$(APP).tex}"

$(DOCB).pdf: 	$(DOC).tex $(APP).tex
	cp $(DOC).tex $(DOCB).tex
	cp $(APP).tex $(APPB).tex
	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
	pdflatex "\def\blind{1} \input{$(APPB).tex}"
	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
	pdflatex "\def\blind{1} \input{$(APPB).tex}"
	bibtex $(DOCB)
	bibtex $(APPB)
	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
	pdflatex "\def\blind{1} \input{$(APPB).tex}"
	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
	pdflatex "\def\blind{1} \input{$(APPB).tex}"

$(FULL).pdf: $(DOC).pdf $(APP).pdf
	pdfunite $(DOC).pdf $(APP).pdf $(FULL).pdf

$(FULLB).pdf: $(DOCB).pdf $(APPB).pdf
	pdfunite $(DOCB).pdf $(APPB).pdf $(FULLB).pdf

