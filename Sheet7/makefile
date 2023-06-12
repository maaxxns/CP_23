.PHONY: all

all: bin/Matrix_N_*.csv

bin/Matrix_N_*.csv: bin/3 | bin 
	./bin/3 

bin/3: src/3.cpp | bin 
	g++ $^ -O3 -o $@

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin