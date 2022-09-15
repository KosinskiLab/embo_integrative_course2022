.PHONY: clean cleaner

all:
	python3 ../../support/tutorial_tools/notebook/process_notebook.py Deposition --branch master

clean:
	rm -f Deposition.md Doxyfile
	rm -rf html

cleaner:
	rm -f Deposition.ipynb Deposition.py Deposition.md Doxyfile
	rm -rf html
