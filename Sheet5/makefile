.PHONY: all

all: build/*.mp4

build/*.mp4: bin/*.csv | build
	python scripts/plot.py

bin/*.csv: bin/2 | bin 
	./bin/2 

bin/2: src/2.cpp | bin 
	g++ $^ -O3 -o $@

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin