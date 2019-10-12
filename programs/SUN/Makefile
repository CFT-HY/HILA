
OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer 
CC = g++ 
LD = g++

TRANSFOMER_GPU = ../../build/transformer -DCUDA  --CUDA
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA 
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

CXXFLAGS = -O3  $(OPTS) 
LDFLAGS = 



SUN_GPU: Makefile SUN.tr.gpu.o monte.gpu.o KennedyPendleton.gpu.o mersenne_inline.gpu.o lattice.gpu.o memory.gpu.o setup_layout_generic.gpu.o map_node_layout_trivial.gpu.o hila_cuda.gpu.o
	$(NVLD) $(LDFLAGS) -o SUN_GPU SUN.tr.gpu.o monte.gpu.o KennedyPendleton.gpu.o mersenne_inline.gpu.o lattice.gpu.o memory.gpu.o setup_layout_generic.gpu.o map_node_layout_trivial.gpu.o hila_cuda.gpu.o
	rm -f *.gpu.o

SUN_CPU: Makefile SUN.tr.cpu.o monte.cpu.o KennedyPendleton.cpu.o mersenne_inline.cpu.o lattice.cpu.o memory.cpu.o setup_layout_generic.cpu.o map_node_layout_trivial.cpu.o
	$(LD) $(LDFLAGS) -o SUN_CPU SUN.tr.cpu.o monte.cpu.o KennedyPendleton.cpu.o mersenne_inline.cpu.o lattice.cpu.o memory.cpu.o setup_layout_generic.cpu.o map_node_layout_trivial.cpu.o
	rm -f *.cpu.o



SUN.tr.cpu.cpp: Makefile SUN.cpp 
	$(TRANSFORMER) SUN.cpp -o=SUN.tr.cpu.cpp

%.cpu.o : %.cpp
	$(CC) $(CXXFLAGS) $< -c -o $@

%.cpu.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c -o $@


SUN.tr.gpu.cpp: Makefile SUN.cpp 
	$(TRANSFOMER_GPU) SUN.cpp -o=SUN.tr.gpu.cpp

%.gpu.o : %.cpp
	$(NVCC) $(CXXFLAGS) $< -c -o $@

%.gpu.o : ../plumbing/%.cpp
	$(NVCC) $(CXXFLAGS) $< -c -o $@


clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.cpt *.o SUN_GPU SUN_CPU

