# Interface Descriptions

Below you will find the necessary class methods that will have to be implemented for each backend. 
This document is also aimed at explaining the structure of the field class that is at the heart of 
the framework. 

## class field_struct

### Methods

### Class Members

* field_storage named payload : basically a class where raw data is stored and which is responsible for retrieving and scattering elements between MPI processes
* lattice_struct name lattice : contains information about the lattice geometry and memory alignment, which is needed by the above class to retrieve elements correctly. Currently this information is passed to payload through the field_struct class -> removes global dependency. 


## class field_storage 

### Methods

* inline T get(const int i, const int field_alloc_size) :  define how elements are retrieved from data buffer 

* inline void set(const T &value, const int i, const int field_alloc_size) : define how elements are set into data buffer.

* void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, parity par) const : gather the elements with parity "par" that are communicated with "to_node" into "buffer". 

* void place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, parity par) : place the elements in buffer, received from "from_node" into their correct places  

* void allocate_field( const int field_alloc_size ) : allocate field memory

* void free_field() : free memory 

Note: field_alloc_size is a variable offered by the lattice_struct class that tells how much memory to allocate in storage. This depends on a number of factors which are handled by the lattice_struct. 







