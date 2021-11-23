#!/bin/bash

rm ../.DS_Store 
rm .DS_Store
rm ../full_ver/.DS_Store
rm figs/.DS_Store
rm figs/2MM/.DS_Store
rm codes/.DS_Store 
pdflatex v1.tex
bibtex v1.bib
bibtex v1.bib
pdflatex v1.tex
rm *.toc *.aux  *.log *.eps *Notes.bib *.blg *.out 
clear
open v1.pdf
