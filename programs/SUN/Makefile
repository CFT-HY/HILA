
OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer 
CC = mpic++ -x c++ -DUSE_MPI
LD = mpic++ -DUSE_MPI

TRANSFORMER_GPU = ../../build/transformer -DGPU -target:CUDA
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA 
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

TRANSFORMER_openacc = ../../build/transformer -Dopenacc --target:openacc
PGCC = pgc++ -acc
PGLD = pgc++ -acc

CXXFLAGS = -O3  $(OPTS) 
LDFLAGS = 


SUN_openacc: Makefile SUN.tr.oac.o monte.oac.o KennedyPendleton.oac.o mersenne_inline.oac.o lattice.oac.o setup_layout_generic.oac.o map_node_layout_trivial.oac.o
	$(PGLD) $(LDFLAGS) -o SUN_openacc SUN.tr.oac.o monte.oac.o KennedyPendleton.oac.o mersenne_inline.oac.o lattice.oac.o setup_layout_generic.oac.o map_node_layout_trivial.oac.o
	rm -f *.oac.o

SUN_GPU: Makefile SUN.tr.gpu.o monte.gpu.o KennedyPendleton.gpu.o mersenne_inline.gpu.o lattice.gpu.o setup_layout_generic.gpu.o map_node_layout_trivial.gpu.o hila_cuda.gpu.o
	$(NVLD) $(LDFLAGS) -o SUN_GPU SUN.tr.gpu.o monte.gpu.o KennedyPendleton.gpu.o mersenne_inline.gpu.o lattice.gpu.o  setup_layout_generic.gpu.o map_node_layout_trivial.gpu.o hila_cuda.gpu.o
	rm -f *.gpu.o

SUN_CPU: Makefile SUN.tr.cpu.o monte.cpu.o KennedyPendleton.cpu.o mersenne_inline.cpu.o lattice.cpu.o setup_layout_generic.cpu.o map_node_layout_trivial.cpu.o com_mpi.cpu.o
	$(LD) $(LDFLAGS) -o SUN_CPU SUN.tr.cpu.o monte.cpu.o KennedyPendleton.cpu.o mersenne_inline.cpu.o lattice.cpu.o setup_layout_generic.cpu.o map_node_layout_trivial.cpu.o com_mpi.cpu.o
	rm -f *.cpu.o



SUN.tr.cpu.cpp: Makefile SUN.cpp 
	$(TRANSFORMER) SUN.cpp -o=SUN.tr.cpu.cpp

%.cpu.o : %.cpp
	$(CC) $(CXXFLAGS) $< -c -o $@

%.cpu.o : %.cpu.cpp
	$(CC) $(CXXFLAGS) $< -c -o $@

%.cpu.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c -o $@


SUN.tr.gpu.cpp: Makefile SUN.cpp 
	$(TRANSFORMER_GPU) SUN.cpp -o=SUN.tr.gpu.cpp

%.gpu.o : %.cpp
	$(NVCC) $(CXXFLAGS) $< -c -o $@

%.gpu.o : %.gpu.cpp
	$(NVCC) $(CXXFLAGS) $< -c -o $@

%.gpu.o : ../plumbing/%.cpp
	$(NVCC) $(CXXFLAGS) $< -c -o $@


SUN.tr.oac.cpp: Makefile SUN.cpp 
	$(TRANSFORMER_openacc) SUN.cpp -o=SUN.tr.oac.cpp

%.oac.o : %.cpp
	$(PGCC) $(CXXFLAGS) -Dopenacc $< -c -o $@

%.oac.o : %.oac.cpp
	$(PGCC) $(CXXFLAGS) -Dopenacc $< -c -o $@

%.oac.o : ../plumbing/%.cpp
	$(PGCC) $(CXXFLAGS) -Dopenacc $< -c -o $@


clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.tr.*.cpp *.o SUN_GPU SUN_CPU SUN_openacc

