.PHONY: all 1 2

all: 1 2

1: build/d.pdf

2: build/rejection.pdf

build/rejection.pdf: bin/2_box_muller.csv scripts/histogram_2.py | build 
	python scripts/histogram_2.py

bin/2_box_muller.csv: bin/2
	./bin/2

bin/2: src/2.cpp | bin
	g++ $^ -O3 -o $@

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