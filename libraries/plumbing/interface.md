# Interface Descriptions

The purpose of this document is to give all the information about the structure of the Field class and it's constituent interfaces that is necessary to implement a new backend. Below you will find a simplified description of the architecture of the Field class and a list of methods that need to be implemented for each new backend.  

## Architecture of the Field class 

To put it briefly, the Field entity is composed of 3 layers. The bottom layer field_storage contains the raw data buffer associated with the Field, and provides the routines for editing this data. The second layer, field_struct, is a wrapper on top of field_storage that adds various variables on MPI communications and a communication initialization routine, and a pointer to the lattice that the whole Field entity relies on. The final layer (class Field) provides the high level user interface, as well as the Field related routines that are written by the hilapp. 

## Methods that need to be implemented for each backend 

Each of these methods are template methods with a type parameter T. 

### field_storage::

* inline T get(const int i, const int field_alloc_size) :  define how elements are retrieved from data buffer 

* inline void set(const T &value, const int i, const int field_alloc_size) : define how elements are set into data buffer.

* void gather_comm_elements(char * buffer, lattice_struct::comm_node_struct to_node, Parity par) const : define how to gather all elements with Parity "par" into "buffer" that will be communicated with "to_node". The comm_node struct contains all necessary information about the destination node.  

* void place_comm_elements(char * buffer, lattice_struct::comm_node_struct from_node, Parity par) : place the elements of Parity par from buffer into their correct places in this field_storage. from_node contains information about the node that sent the information. 

* void allocate_field( const int field_alloc_size ) : allocate the data buffer.

* void free_field() : free the memory in the data buffer.  


## notes:

Different backends are selected and added to the code with preprocessor directives. After the methods above have been defined in a header file, add an include statement with this header file after a suitable preprocessor directive after the field_struct interface has been declared. Note that for now all backends are mutually exclusive. 







