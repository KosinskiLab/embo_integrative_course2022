.PHONY: clean cleaner

all:
	python3 ../../support/tutorial_tools/notebook/process_notebook.py analysis --branch master

clean:
	rm -f analysis.md Doxyfile
	rm -rf html

cleaner:
	rm -f analysis.ipynb analysis.sh analysis.md Doxyfile
	rm -rf html
