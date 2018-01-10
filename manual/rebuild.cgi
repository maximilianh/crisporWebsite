#!/bin/sh

pandoc manual.md -o manual.epub
pandoc manual.md -o manual.pdf --variable=fontfamily:bera
pandoc manual.md -o manual.html --css pandoc.css --toc -s
