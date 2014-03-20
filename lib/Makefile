
CXXFLAGS += -std=c++0x

all: libcop.a

libcop.a: coplib.o
	$(AR) r $@ $<

clean:
	$(RM) libcop.a *.o

