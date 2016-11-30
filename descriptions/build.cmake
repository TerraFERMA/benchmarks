execute_process(COMMAND rm tfbenchmarks.aux tfbenchmarks.toc tfbenchmarks.out tfbenchmarks.pdf tfbenchmarks.log tfbenchmarks.blg tfbenchmarks.bbl)
execute_process(COMMAND pdflatex tfbenchmarks.tex)
execute_process(COMMAND pdflatex tfbenchmarks.tex)
execute_process(COMMAND bibtex tfbenchmarks.aux)
execute_process(COMMAND bibtex tfbenchmarks.aux)
execute_process(COMMAND pdflatex tfbenchmarks.tex)
execute_process(COMMAND pdflatex tfbenchmarks.tex)


