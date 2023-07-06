.PHONY: all 1 2

all: 1 2

2: build/ising.gif

build/ising.gif: bin/state.csv scripts/plot.py | build 
	python scripts/plot.py

bin/state.csv: bin/2 | bin 
	./bin/2 

bin/2: src/2.cpp | bin 
	g++ $^ -O3 -o $@

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin