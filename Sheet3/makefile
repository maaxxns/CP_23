
.PHONY: all
all: bin/MD-template

bin/MD-template: src/MD-template.cpp
	g++ $^ -o $@

bin:
	mkdir bin

clean:
	rm -rf bin

