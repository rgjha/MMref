#!/bin/bash

rm ../.DS_Store 
rm .DS_Store
rm ../full_ver/.DS_Store
rm figs/.DS_Store
rm figs/2MM/.DS_Store
rm codes/.DS_Store 
pdflatex v2.tex
bibtex v2.bib
bibtex v2.bib
pdflatex v2.tex
rm *.toc *.aux  *.log *.eps *Notes.bib *.blg *.out 
clear
open v2.pdf
