RSCRIPT=/psi/R-3.0.1/R-3.0.1/bin/Rscript

all: evaluation.Rmd
	$(RSCRIPT) -e "require(knitr); knit2html(\"evaluation.Rmd\")"

pdf: evaluation.md
	pandoc -s evaluation.md -t latex -o evaluation.pdf