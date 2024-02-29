tex=$(wildcard *.tex)
pdf=$(tex:.tex=.pdf)
jpeg=$(pdf:.pdf=.jpeg)
density=400
quality=100

.PHONY: all
all: $(jpeg)

early-int.jpeg: early-int.pdf
	convert -density $(density) -quality $(quality) $< $@
late-int.jpeg: late-int.pdf
	convert -density $(density) -quality $(quality) $< $@
test-train.jpeg: test-train.pdf
	convert -density $(density) -quality $(quality) $< $@

early-int.pdf: early-int.tex
	tectonic $<
late-int.pdf: late-int.tex
	tectonic $<
test-train.pdf: test-train.tex
	tectonic $<

.PHONY: clean
clean:
	rm -f $(pdf)
