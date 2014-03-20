
LIBDIR := ../lib
CXXFLAGS += -std=c++0x -I$(LIBDIR)
LDFLAGS += -L$(LIBDIR)
LIBS := -lcop

MODULES := 03_doublependulum

all: $(MODULES)

$(MODULES) : % : %.o $(LIBDIR)/libcop.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@


#.PHONY: $(LIBDIR)/libcop.a
$(LIBDIR)/libcop.a: $(wildcard $(LIBDIR)/*.cpp wildcard $(LIBDIR)/*.hpp)
	$(MAKE) -C $(LIBDIR)


plots: $(addsuffix -plot, $(MODULES))


.PHONY: $(addsuffix -plot, $(MODULES))
$(addsuffix -plot, $(MODULES)) : %-plot : %.gpt %.dat
	gnuplot $<


$(addsuffix .dat, $(MODULES)) : %.dat : %.prm %
	./$* > $@


clean:
	$(RM) *.png
	$(RM) *.dat
	$(RM) $(MODULES) *.o

