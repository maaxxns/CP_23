
.PHONY: all
all: plots/plot.pdf

bin/MD-template: src/MD-template.cpp | bin
	g++ $^ -o $@

bin/b)g.dat: bin/MD-template | bin
	./bin/MD-template

plots/plot.pdf: scripts/plots.py bin/b)g.dat | plots
	python scripts/plots.py
bin:
	mkdir bin

plots:
	mkdir plots

clean:
	rm -rf bin

