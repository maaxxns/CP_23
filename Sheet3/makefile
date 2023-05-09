
.PHONY: all
all: bin/MD-template

bin/MD-template: src/MD-template.cpp | bin
	g++ $^ -o $@

bin:
	mkdir bin

clean:
	rm -rf bin

