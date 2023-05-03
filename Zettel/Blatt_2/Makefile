ODIR	= build

.PHONY: all
all: $(ODIR)/Ex3

plots: plot.py Ex3.cpp | $(ODIR)
	python plot.py 

$(ODIR)/Ex3: Ex3.cpp | $(ODIR)
	g++ $< -c -o $@



$(ODIR):
	mkdir $(ODIR)