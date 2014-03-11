
all: plot

plot: sine.dat plot.gpt
	gnuplot plot.gpt

sine.dat: writesin
	./writesin > $@

clean:
	$(RM) writesin
	$(RM) sine.dat sine.png
