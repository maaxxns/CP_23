.PHONY: all

all: build/Diffusion_c_weird.mp4

build/Diffusion_c_weird.mp4: bin/*.csv | build
	python scripts/plot.py

bin/*.csv: bin/1 | bin 
	./bin/1 

bin/2: src/Ex1.cpp | bin 
	g++ $^ -O3 -o $@

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin