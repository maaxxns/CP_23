.PHONY: 2

2: build/Lorentz_r_20_3d.pdf

build/Lorentz_r_20_3d.pdf: bin/Lorentz_r_is_28.csv scripts/plot.py | build 
	python scripts/plot.py

bin/Lorentz_r_is_28.csv: bin/2 | bin 
	./bin/2 

bin/2: src/2.cpp | bin 
	g++ $^ -O3 -o $@

bin:
	mkdir bin

build:
	mkdir build

clean:
	rm -rf bin