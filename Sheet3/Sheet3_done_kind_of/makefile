
.PHONY: all
all: plots/plot_E_2561.pdf

bin/MD-template: src/MD-template.cpp | bin
	g++ $^ -o $@

bin/n8.dat: bin/MD-template | bin
	./bin/MD-template

plots/plot_E_2561.pdf: scripts/plots.py bin/bN2561.dat | plots
	python scripts/plots.py
	
bin:
	mkdir bin

plots:
	mkdir plots

clean:
	rm -rf bin

