#!/bin/bash


rm .DS_Store
rm figs/.DS_Store
rm codes/.DS_Store 
pdflatex notes1.tex
bibtex notes1.bib
bibtex notes1.bib
bibtex notes1.bib
pdflatex notes1.tex
bibtex notes1.bib
pdflatex notes1.tex
pdflatex notes1.tex
rm *.toc *.aux  *.log *.eps *Notes.bib *.blg *.out 
clear
open notes1.pdf
