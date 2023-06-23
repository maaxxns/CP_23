.PHONY: all 1

all: 1 2

1: build/newton.pdf

build/newton.pdf: bin/newton.csv scripts/plot.py | build 
	python scripts/plot.py

bin/newton.csv: bin/1 | bin 
	./bin/1 

bin/1: src/1.cpp | bin 
	g++ $^ -O3 -o $@

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin