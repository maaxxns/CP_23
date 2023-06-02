.PHONY: all
 

plots: bin/*.csv | build
	python scripts/plot.py

bin/*.csv: bin/Ex1 | build
	./bin/Ex1 

bin/Ex1: src/Ex1.cpp | bin 
	g++ $^ -O3 -o $@

 
bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin
	rm -rf build