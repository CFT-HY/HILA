ADDITION = -DMSIZE=2
OPTS = -DNDIM=4 $(ADDITION) -DSUBNODE_LAYOUT -I../../libraries/


#The next line work on puhti
HILAPP = ../../hilapp/build/hilapp -DPUHTI -DSUBNODE_LAYOUT -DHILAPP -I../../libraries/ -I../../../llvm/lib/clang/8.0.1/include/ 
#HILAPP = ../../hilapp/build/hilapp  -DUSE_MPI -I../../libraries/ -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ 

CC = mpic++ -DUSE_MPI -mavx2 -mfma  -I../../libraries/

CFLAGS = -Ofast -I../../libraries/ $(OPTS)
CXXFLAGS = -Ofast -ffast-math -x c++ --std=c++17 -I../../libraries/ $(OPTS)
# CXXFLAGS = -g -x c++ --std=c++17 $(OPTS)

# Don't remove cpt-files
.PRECIOUS: %.cpt

%.exe: %.tr.o Makefile inputs.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o timers.o memalloc.o test_gathers.o Makefile
	$(CC) -o $@ $< inputs.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o timers.o memalloc.o test_gathers.o

%.cpt: %.cpp Makefile
	$(HILAPP) $(ADDITION) $<

%.tr.o : %.cpt Makefile
	$(CC) $(CXXFLAGS) $< -c -o $@

%.tr.s : %.cpt Makefile
	$(CC) $(CXXFLAGS) $< -S -o $@

lattice.cpt : ../../libraries/plumbing/lattice.cpp Makefile
	$(HILAPP) $(ADDITION) $<

test_gathers.cpt : ../../libraries/plumbing/test_gathers.cpp Makefile
	$(HILAPP) $(ADDITION) $<

lattice.o : lattice.cpt Makefile
	$(CC) $(CXXFLAGS) $< -c

test_gathers.o : test_gathers.cpt Makefile
	$(CC) $(CXXFLAGS) $< -c 

%.o : ../../libraries/plumbing/%.cpp Makefile
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c Makefile
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.o *.cpt *.tr.s

cleanall:
	rm -rf *.o *.cpt *.exe *.tr.s
