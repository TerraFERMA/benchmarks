execute_process(
     COMMAND pdflatex tfbenchmarks.tex && pdflatex tfbenchmarks.tex && bibtex tfbenchmarks.aux && bibtex tfbenchmarks.aux && pdflatex
tfbenchmarks.tex && pdflatex tfbenchmarks.tex
     )


