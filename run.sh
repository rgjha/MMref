#!/bin/bash

rm figs/.DS_Store
rm codes/.DS_Store 
bibtex notes1.bib
pdflatex notes1.tex
bibtex notes1.bib
pdflatex notes1.tex
pdflatex notes1.tex
rm *.toc *.aux  *.log *.eps *Notes.bib *.blg *.out 
clear
open notes1.pdf
