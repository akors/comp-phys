
LIBDIR := ../lib
CXXFLAGS += -std=c++0x -I$(LIBDIR)
LDFLAGS += -L$(LIBDIR)
LIBS := -lcop

REMOTEHOST := Calgary
REMOTEPATH := ~/comp-phys/blatt2

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

.PHONY: $(addprefix remotecalc-, $(MODULES))
$(addprefix remotecalc-, $(MODULES)) : remotecalc-% : %.prm %
	scp $^ $(REMOTEHOST):$(REMOTEPATH)
	ssh $(REMOTEHOST) 'cd $(REMOTEPATH) && ./$* | bzip2 ' | bunzip2 > $*.dat


clean:
	$(RM) *.png
	$(RM) *.dat
	$(RM) $(MODULES) *.o
	$(MAKE) -C $(LIBDIR) $@

