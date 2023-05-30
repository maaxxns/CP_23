.PHONY: all 1

all: build/wave.gif build/schroedinger.gif
1: build/schroedinger.gif
2: build/wave.gif

build/wave.gif: bin/wiggle_that_membrane.csv | build
	python scripts/plot.py

bin/wiggle_that_membrane.csv: bin/2 | bin 
	./bin/2 

bin/2: src/2.cpp | bin 
	g++ $^ -O3 -o $@

bin/1: src/1.cpp | bin 
	g++ $^ -O3 -o $@

bin/schroedinger.csv: bin/1 | bin
	./bin/1 

build/schroedinger.gif: bin/schroedinger.csv scripts/plot1.py | build 
	python scripts/plot1.py

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin
	rm -rf build