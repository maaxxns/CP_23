.PHONY: all

all: build/d.pdf

build/d.pdf: bin/rng_number_d.csv scripts/histogram_1.py | build
	python scripts/histogram_1.py

bin/rng_number_d.csv: bin/1 | bin 
	./bin/1 

bin/1: src/1.cpp | bin 
	g++ $^ -O3 -o $@

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin
	rm -rf build