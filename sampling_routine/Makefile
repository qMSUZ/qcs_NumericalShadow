#NVCC = nvcc -ccbin g++-4.9
NVCC = nvcc
CXX = g++

CXXFLAGS = -I. -I/opt/include
#NVCCFLAGS = -g -G -m64 -I. --gpu-architecture=sm_20
NVCCFLAGS = -m64 -I. --gpu-architecture=sm_20 -Xcompiler '-fPIC'
#LDFLAGS = -L/usr/local/cuda/lib64/
LDFLAGS = 
LIBS = -lcuda -lcudart -lcudadevrt

all: libSamplingRoutine.so samplingroutinetest

samplingroutine.o: samplingroutine.cu
	$(NVCC) $(NVCCFLAGS) --device-c samplingroutine.cu -o samplingroutine.o

cplxNum.o: cplxNum.cu
	$(NVCC) $(NVCCFLAGS) --device-c cplxNum.cu -o cplxNum.o

samplingroutinetest.o: samplingroutinetest.cu
	$(NVCC) -c $(NVCCFLAGS) samplingroutinetest.cu -o samplingroutinetest.o

libSamplingRoutine.so:  samplingroutine.o cplxNum.o
	$(NVCC) -shared  samplingroutine.o cplxNum.o  -o libSamplingRoutine.so -I.

samplingroutinetest: samplingroutinetest.o samplingroutine.o cplxNum.o
	$(NVCC) $(NVCCFLAGS) -o samplingroutinetest samplingroutinetest.o samplingroutine.o cplxNum.o  $(LDFLAGS) $(LIBS)

clean:
	rm -f samplingroutinetest libSamplingRoutine.so samplingroutinetest.o samplingroutine.o cplxNum.o



# 	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) $(LIBS) --compiler-options '-fPIC' samplingroutine.cu -shared --compiler-options '-Wl,-soname,libSamplingRoutine.so' -o libSamplingRoutine.so 

