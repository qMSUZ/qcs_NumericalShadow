#NVCC = nvcc -ccbin g++-4.9
NVCC = nvcc
CXX = g++

CXXFLAGS = -I. -I/opt/include
#NVCCFLAGS = -g -G -m64 -I. --gpu-architecture=sm_20
NVCCFLAGS = -m64 -I. --gpu-architecture=sm_20
#LDFLAGS = -L/usr/local/cuda/lib64/
LDFLAGS = 
LIBS = -lcuda -lcudart -lcudadevrt

all: libSamplingRoutine.so samplingroutinetest


libSamplingRoutine.so: 
	$(NVCC) -shared -Xcompiler '-fPIC' samplingroutine.cu  -o libSamplingRoutine.so -I.

samplingroutinetest.o:
	$(NVCC) -c $(NVCCLAGS) samplingroutinetest.cu -o samplingroutinetest.o

samplingroutinetest: samplingroutinetest.o
	$(NVCC) $(NVCCFLAGS) -o samplingroutinetest samplingroutinetest.o $(LDFLAGS) $(LIBS)

clean:
	rm -f samplingroutinetest libSamplingRoutine.so samplingroutinetest.o samplingroutine.o



# 	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) $(LIBS) --compiler-options '-fPIC' samplingroutine.cu -shared --compiler-options '-Wl,-soname,libSamplingRoutine.so' -o libSamplingRoutine.so 

