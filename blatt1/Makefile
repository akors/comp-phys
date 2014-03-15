
CXXFLAGS += -std=c++0x

MODULES := 01_sine 02_pendulum

all: $(MODULES)

plots: $(addsuffix -plot, $(MODULES))

.PHONY: $(addsuffix -plot, $(MODULES))

$(addsuffix -plot, $(MODULES)) : %-plot : %.gpt %.dat
	gnuplot $<


$(addsuffix .dat, $(MODULES)) : %.dat : %.prm %
	./$* > $@


clean:
	$(RM) *.png
	$(RM) *.dat
	$(RM) $(MODULES) 

