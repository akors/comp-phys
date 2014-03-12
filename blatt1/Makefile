
CXXFLAGS += -std=c++0x

all: plot pendulum

plot: sine.dat plot.gpt
	gnuplot plot.gpt

sine.dat: writesin
	./writesin > $@

clean:
	$(RM) writesin pendulum
	$(RM) sine.dat sine.png
