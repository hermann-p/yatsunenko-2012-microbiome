RSCRIPT=R
PANDOC=~/.cabal/bin/pandoc

all: evaluation.md

evaluation.md: evaluation.Rmd
	$(RSCRIPT) -e "require(knitr); knit2html(\"evaluation.Rmd\")"

pdf: evaluation.md
	$(PANDOC) -s -S --biblio evaluation.bib --csl citationstyle.csl evaluation.md -t latex -o evaluation.pdf
