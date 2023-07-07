.PHONY: all
 
Ex	= 2

all: bin/*.csv

plots: bin/*.csv | build
	python scripts/plot$(Ex).py

bin/*.csv: bin/Ex$(Ex) | build 
	./bin/Ex$(Ex) 

bin/Ex$(Ex): src/Ex$(Ex).cpp | bin 
	g++ $^ -O3 -o $@

 
bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin
	rm -rf build