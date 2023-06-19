.PHONY: all
 
Ex	= Ex1

all: bin/*.csv 

plots: bin/*.csv | build
	python scripts/plot.py

bin/*.csv: bin/$(Ex) | build 
	./bin/$(Ex) 

bin/$(Ex): src/$(Ex).cpp | bin 
	g++ $^ -O3 -o $@

 
bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin
	rm -rf build