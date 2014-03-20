
CXXFLAGS += -std=c++0x

all: libcop.a

coplib.o: coplib.cpp coplib.hpp

libcop.a: coplib.o
	$(AR) r $@ $<

clean:
	$(RM) libcop.a *.o

