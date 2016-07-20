CXX = g++

CXXFLAGS = -I. -I/opt/include -std=c++11
LDFLAGS  = -L/opt/lib
LFPIC = -fPIC
LIBS = -lmagma -lcublas -lcusparse -lpthread -lcudart -fopenmp -lopenblas

all: libNumRange numrangetest


numrange_test.o: numrange_test.cpp
	$(CXX) -c $(CXXFLAGS) numrange_test.cpp 

numrange.o: numrange.cpp numrange.h
	$(CXX) -c $(CXXFLAGS) $(LFPIC) numrange.cpp 

libNumRange: numrange.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LFPIC) numrange.o -shared -Wl,-soname,libNumRange.so -o libNumRange.so $(LIBS)

numrangetest: numrange_test.o numrange.o
	$(CXX) -o numrangetest $(CXXFLAGS) $(LDFLAGS) numrange_test.o numrange.o -lmagma -lcublas -lcusparse -std=c++11 -lpthread -lcudart -fopenmp -lopenblas 

clean:
	rm -f numrangetest libNumRange.so numrange_test.o numrange.o

 
