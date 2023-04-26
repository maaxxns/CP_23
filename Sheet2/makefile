
.PHONY: all
all: plots/plote.pdf


bin/02_03: src/02_03.cpp | bin
	g++ $^ -o $@

bin/euler.csv: bin/02_03 | bin 
	./bin/02_03

plots/plote.pdf: scripts/plot3.py bin/euler.csv | bin plots
	python scripts/plot3.py

plots:
	mkdir plots

bin:
	mkdir bin

clean:
	rm -rf bin
	rm -rf plots
