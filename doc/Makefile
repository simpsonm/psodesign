# LaTeX Makefile

## Define the file name (.tex)

DOC=psokriging
APP=psokrigingapp
# DOCB=psokrigingblind
# APPB=psokrigingappblind
# FULL=psokrigingfull
# FULLB=psokrigingblindfull
FIG=psokriginggraphs

all:	$(DOC).pdf $(FIG).pdf #$(DOCB).pdf $(FULL).pdf $(FULLB).pdf #

.PHONY: clean

clean:
	\rm *.aux *.blg *.out *.bbl *.log $(DOC).pdf $(APP).pdf $(FIG).pdf

 #$(DOCB).pdf  $(APPB).pdf $(FULL).pdf $(FULLB).pdf $(DOCB).tex $(APPB).tex 
 ##$(FIG).pdf $(FIG).tex

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

# $(DOCB).pdf: 	$(DOC).tex $(APP).tex
# 	cp $(DOC).tex $(DOCB).tex
# 	cp $(APP).tex $(APPB).tex
# 	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
# 	pdflatex "\def\blind{1} \input{$(APPB).tex}"
# 	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
# 	pdflatex "\def\blind{1} \input{$(APPB).tex}"
# 	bibtex $(DOCB)
# 	bibtex $(APPB)
# 	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
# 	pdflatex "\def\blind{1} \input{$(APPB).tex}"
# 	pdflatex "\def\blind{1} \input{$(DOCB).tex}"
# 	pdflatex "\def\blind{1} \input{$(APPB).tex}"

# $(FULL).pdf: $(DOC).pdf $(APP).pdf
# 	pdfunite $(DOC).pdf $(APP).pdf $(FULL).pdf

# $(FULLB).pdf: $(DOCB).pdf $(APPB).pdf
# 	pdfunite $(DOCB).pdf $(APPB).pdf $(FULLB).pdf

$(FIG).pdf: $(DOC).pdf
	cp $(DOC).tex $(FIG).tex
	pdflatex '\PassOptionsToPackage{active}{preview} \input{$(FIG).tex}'
# ## run this after in shell; can't get to work in makefile
# 	P=$(pdftk $(FIG).pdf dump_data|grep NumberOfPages| awk '{print $2}')	
# 	for i in $(seq 1 1 $P)
# 	do j=`expr $i - 1`
# 	echo $j
# 	convert -units PixelsPerInch -density 800 $(FIG).pdf[$j] -fuzz 1% -trim +repage $(FIG)-$i.tiff
# 	done
