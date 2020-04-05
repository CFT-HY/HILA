


#include "../plumbing/globals.h"
#include "../plumbing/lattice.h"
#include "../plumbing/field.h"


void test_std_gathers();


///***********************************************************
/// setup() lays out the lattice infrastruct, with neighbour arrays etc.

/// A list of all defined lattices
std::vector<lattice_struct*> lattices;


/// General lattice setup, including MPI setup
void lattice_struct::setup(int siz[NDIM], int &argc, char **argv) {
  // Add this lattice to the list 
  lattices.push_back( this );

  l_volume = 1;
  for (int i=0; i<NDIM; i++) {
    l_size[i] = siz[i];
    l_volume *= siz[i];
  }

#ifdef USE_MPI
  /* Initialize MPI */
  initialize_machine(argc, &argv);

  /* default comm is the world */
  mpi_comm_lat = MPI_COMM_WORLD;

  MPI_Comm_rank( lattices[0]->mpi_comm_lat, &this_node.rank );
  MPI_Comm_size( lattices[0]->mpi_comm_lat, &nodes.number );

#else 
  this_node.rank = 0;
  nodes.number = 1;
#endif

  setup_layout();

  setup_nodes();

  // set up the comm arrays 
  create_std_gathers();

  // Initialize wait_array structures - has to be after std gathers()
  initialize_wait_arrays();

  test_std_gathers();


}


#if NDIM==4
void lattice_struct::setup(int nx, int ny, int nz, int nt, int &argc, char **argv) {
  int s[NDIM] = {nx, ny, nz, nt};
  setup(s, argc, argv);
}
#elif NDIM==3
void lattice_struct::setup(int nx, int ny, int nz, int &argc, char **argv) {
  int s[NDIM] = {nx, ny, nz};
  setup(s, argc, argv);
}
#elif NDIM==2
void lattice_struct::setup(int nx, int ny, int &argc, char **argv) {
  int s[NDIM] = {nx, ny};
  setup(s, argc, argv);
}
#elif NDIM==1
void lattice_struct::setup(int nx, int &argc, char **argv) {
  int s[NDIM] = {nx};
  setup(s, argc, argv);
}
#endif

///////////////////////////////////////////////////////////////////////
/// The routines is_on_node(), node_rank(), site_index()
/// implement the "layout" of the nodes and sites of the lattice.
/// To be changed in different implementations!
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
/// Get the node rank for coordinates 
/// This is the fundamental routine which defines how the nodes
/// are mapped.  map_node_layout MUST BE compatible with this
/// algorithm!  So if you change this, change that too
/// 
/// Here the node number along one direction is calculated with
///    loc * n_divisions/l_size  (with integer division)
/// example: with l_size=14 and n_divisions=3,
/// the dividers are at 0, 5, 10, 14
/// 
///////////////////////////////////////////////////////////////////////

int lattice_struct::node_rank(const coordinate_vector & loc)
{
  int i;
  int dir;

  i = (loc[NDIM-1] * nodes.n_divisions[NDIM-1]) / l_size[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    i = i*nodes.n_divisions[dir] +
        ((loc[dir] * nodes.n_divisions[dir]) / l_size[dir]);
  }
  /* do we want to remap this?  YES PLEASE */
  i = nodes.remap(i);

  return( i );
}

///////////////////////////////////////////////////////////////////////
/// Is the coordinate on THIS node 
///////////////////////////////////////////////////////////////////////

bool lattice_struct::is_on_node(const coordinate_vector & loc)
{
  int d;

  for (int dir=0; dir<NDIM; dir++) {
    d = loc[dir] - this_node.min[dir];
    if (d < 0 || d >= this_node.size[dir] ) return false;
  }
  return true;
}


///////////////////////////////////////////////////////////////////////
/// give site index for ON NODE sites
/// Note: loc really has to be on this node
///////////////////////////////////////////////////////////////////////

unsigned lattice_struct::site_index(const coordinate_vector & loc)
{
  int dir,l,s;
  unsigned i;

  i = l = loc[NDIM-1] - this_node.min[NDIM-1];
  s = loc[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    l = loc[dir] - this_node.min[dir];
    i = i*this_node.size[dir] + l;
    s += loc[dir];
  }

  // now i contains the `running index' for site
#if defined(EVEN_SITES_FIRST)
  if (s%2 == 0) return( i/2 );    /* even site index */
  else return( i/2 + this_node.evensites );  /* odd site */
#else
  return( i );
#endif
}

///////////////////////////////////////////////////////////////////////
/// give site index for nodeid sites
/// compare to above
///////////////////////////////////////////////////////////////////////

unsigned lattice_struct::site_index(const coordinate_vector & loc, const unsigned nodeid)
{
  int dir,l,s;
  unsigned i;
  const node_info & ni = nodes.nodelist[nodeid];
  
  i = l = loc[NDIM-1] - ni.min[NDIM-1];
  s = loc[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    l = loc[dir] - ni.min[dir];
    i = i*ni.size[dir] + l;
    s += loc[dir];
  }

  // now i contains the `running index' for site
#if defined(EVEN_SITES_FIRST)
  if (s%2 == 0) return( i/2 );    /* even site index */
  else return( i/2 + ni.evensites );  /* odd site */
#else
  return( i );
#endif
}


///////////////////////////////////////////////////////////////////////
/// invert the this_node index -> location (only on this node)
///////////////////////////////////////////////////////////////////////

// this is defined in lattice.h

//const coordinate_vector & lattice_struct::coordinates(unsigned index)
//{
//  return this_node.coordinates[index];
//}


/////////////////////////////////////////////////////////////////////
/// Set up the basic list about all nodes
/// NOTE: in case of very large number of nodes this may be
/// inconveniently large.  Alternatives?
/////////////////////////////////////////////////////////////////////

void lattice_struct::setup_nodes() {

  nodes.number = numnodes();
  // Loop over all node "origins"
  nodes.nodelist.resize(nodes.number);

  // n keeps track of the node "root coordinates"
  coordinate_vector n(0);

  // use nodes.divisors - vectors to fill in stuff
  for (int i=0; i<nodes.number; i++) {
    coordinate_vector l;
    foralldir(d) l[d] = nodes.divisors[d][n[d]];

    int nn = node_rank(l);
    node_info & ni = nodes.nodelist[nn];
    int v = 1;
    foralldir(d) {
      ni.min[d]  = nodes.divisors[d][n[d]];
      ni.size[d] = nodes.divisors[d][n[d]+1] - nodes.divisors[d][n[d]];
      v *= ni.size[d];
    }
    if (v % 2 == 0)
      ni.evensites = ni.oddsites = v/2;
    else {
      // now node ni has odd number of sites
      if (l.coordinate_parity() == EVEN) {
        ni.evensites = v/2 + 1; 
        ni.oddsites = v/2;
      } else {
        ni.evensites = v/2; 
        ni.oddsites = v/2 + 1;
      }
    }
    
    // now to the next divisor - add to the lowest dir until all done
    foralldir(d) {
      n[d]++;
      if (n[d] < nodes.n_divisions[d]) break;
      n[d] = 0;
    }
    
    // use the opportunity to set up this_node when it is met
    if (nn == mynode()) this_node.setup(ni, *lattice);

    
  }
}


////////////////////////////////////////////////////////////////////////
/// Fill in this_node fields -- node_rank() must be set up OK
////////////////////////////////////////////////////////////////////////
void lattice_struct::node_struct::setup(node_info & ni, lattice_struct & lattice)
{
  
  rank = mynode();

  min  = ni.min;
  size = ni.size;

  evensites = ni.evensites;
  oddsites  = ni.oddsites;
  sites     = ni.evensites + ni.oddsites;

  first_site_even = (min.coordinate_parity() == EVEN);
   
  // neighbour node indices
  foralldir(d) {
    coordinate_vector l = min;   // this is here on purpose
    l[d] = (min[d] + size[d]) % lattice.l_size[d];
    nn[d] = lattice.node_rank(l);
    l[d] = (lattice.l_size[d] + min[d] - 1) % lattice.l_size[d];
    nn[opp_dir(d)] = lattice.node_rank(l);
  }

  // map site indexes to locations -- coordinates array
  coordinates.resize(sites);
  coordinate_vector l = min;
  for(unsigned i = 0; i<sites; i++){
    coordinates[ lattice.site_index(l) ] = l;
    // walk through the coordinates
    foralldir(d) {
      if (++l[d] < (min[d]+size[d])) break;
      l[d] = min[d];
    }
  }
}

/////////////////////////////////////////////////////////////////////
/// Create the neighbour index arrays 
/// This is for the index array neighbours
/// TODO: implement some other neighbour schemas!
/////////////////////////////////////////////////////////////////////

void lattice_struct::create_std_gathers()
{

  // allocate neighbour arrays - TODO: these should 
  // be allocated on "device" memory too!

  for (int d=0; d<NDIRS; d++) {
    neighb[d] = (unsigned *)memalloc(this_node.sites * sizeof(unsigned));
  }
  
  unsigned c_offset = this_node.sites;  // current offset in field-arrays

#ifdef SCHROED_FUN
  // special case for Schroedinger functional: fixed boundary
  // This requires we allocate one extra space for
  // CONSTANT schroedinger b.c.
  // Thus, add 1 item to latfield arrays to
  // take the boundary into account, and set the
  // neighb[d] to point to it.  This is for mpi or vanilla.
  // NOTE: IF NON-CONSTANT BOUNDARY NEED TO ADD c_offset FOR
  // EACH SITE, and set the corresponding neighb

  int sf_special_boundary = c_offset;
  /* add space for boundary only when needed */
  if (this_node.min[NDIM-1] + this_node.nodesize[NDIM-1] == size(NDIM-1])
    c_offset += 1;
  else sf_special_boundary = -(1<<30);

  output0 << "SPECIAL BOUNDARY LAYOUT for SF";
#endif


  // We set the communication and the neigbour-array here

  for (direction d=XUP; d<NDIRS; ++d) {

    nn_comminfo[d].index = neighb[d];    // this is not really used for nn gathers

    comm_node_struct & from_node = nn_comminfo[d].from_node;
    // we can do the opposite send during another pass of the sites. 
    // This is just the gather inverted
    // NOTE: this is not the send to direction d, but to -d!    
    comm_node_struct & to_node = nn_comminfo[-d].to_node;

    from_node.rank = to_node.rank = this_node.rank;    // invalidate from_node, for time being
    // if there are no communications the rank is left as is

    // counters to zero
    from_node.sites = from_node.evensites = from_node.oddsites = 0;

    // pass over sites
    int num = 0;  // number of sites off node
    for (int i=0; i<this_node.sites; i++) {
      coordinate_vector ln, l;
      l = coordinates(i);
      // set ln to be the neighbour of the site
      ln = mod(l + d, size());
      //ln = l;
      //if (is_up_dir(d)) ln[d] = (l[d] + 1) % size(d);
      //else ln[d] = (l[d] + size(-d) - 1) % size(-d);
 
#ifdef SCHROED_FUN
      if (d == NDIM-1 && l[NDIM-1] == size(NDIM-1)-1) {
	      // This is up-direction, give special site
    	  neighb[d][i] = sf_special_boundary;
      } else if (d == opp_dir(NDIM-1) && l[NDIM-1] == 0) {
	      //  We never should need the down direction, thus, short circuit!
	      neighb[d][i] = 1<<30;
      } else     // NOTE THIS UGLY else-if CONSTRUCT!
#endif
      if (is_on_node(ln)) {
        neighb[d][i] = site_index(ln);
      } else {
        // reset neighb array temporarily, as a flag
        neighb[d][i] = this_node.sites;

	      // Now site is off-node, this lead to fetching
        // check that there's really only 1 node to talk with
        unsigned rank = node_rank(ln);
        if ( from_node.rank == this_node.rank) {
          from_node.rank = rank;
        } else if (from_node.rank != rank) {
          hila::output << "Internal error in nn-communication setup\n";
          exit(-1);
        }

        from_node.sites++;
        if (l.coordinate_parity() == EVEN) from_node.evensites++;
        else from_node.oddsites++;

	      num++;
      }
    }

    // and set buffer indices
    from_node.buffer = c_offset;

    to_node.rank = from_node.rank;
    to_node.sites = from_node.sites;
    // note!  Parity is flipped, because parity is normalized to receieving node
    to_node.evensites = from_node.oddsites;
    to_node.oddsites = from_node.evensites;

    if (num > 0) {
      // sitelist tells us which sites to send
      to_node.sitelist = (unsigned *)memalloc(to_node.sites * sizeof(unsigned));
#ifndef VANILLA
      // non-vanilla code wants to have receive buffers, so we need mapping to field
      from_node.sitelist = (unsigned *)memalloc(from_node.sites * sizeof(unsigned));
#endif
    } else 
      to_node.sitelist = nullptr;

    if (num > 0) {
      // set the remaining neighbour array indices and sitelists in another go over sites.
      // temp counters
      // NOTE: ordering is automatically right: with a given parity,
      // neighbour node indices come in ascending order of host node index - no sorting needed
      int c_even, c_odd;
      c_even = c_odd = 0;

      for (int i=0; i<this_node.sites; i++) {
        if (neighb[d][i] == this_node.sites) {
          coordinate_vector ln, l;
          l = coordinates(i);
          ln = mod(l + d, size()); 

          if (l.coordinate_parity() == EVEN) {
            // THIS site is even
            neighb[d][i] = c_offset + c_even;

#ifndef VANILLA
            from_node.sitelist[c_even] = i;
#endif

            // flipped parity: this is for odd sends
            to_node.sitelist[c_even + to_node.evensites] = i;

            c_even++;

          } else {
            neighb[d][i] = c_offset + from_node.evensites + c_odd;

#ifndef VANILLA
            from_node.sitelist[c_odd + from_node.evensites] = i;
#endif

            // again flipped parity for setup
            to_node.sitelist[c_odd] = i;

            c_odd++;
          }
        }
      }
    }

    c_offset += from_node.sites;

  } /* directions */

  /* Finally, set the site to the final offset (better be right!) */
  this_node.field_alloc_size = c_offset;

  /* Setup backend-specific lattice info if necessary */
  backend_lattice = new backend_lattice_struct;
  backend_lattice->setup(*this);

}



///  Get message tags cyclically -- defined outside classes, so that it is global and unique

#define MSG_TAG_MIN  10000
#define MSG_TAG_MAX  (1<<28)     // a big int number

int get_next_msg_tag() {
  static int tag = MSG_TAG_MIN;
  ++tag;
  if (tag > MSG_TAG_MAX) tag = MSG_TAG_MIN;
  return tag;
}





/************************************************************************/

#ifdef USE_MPI
/* this formats the wait_array, used by forallsites_waitA()
 * should be made as fast as possible!
 *
 * wait_array[i] contains a bit at position 1<<dir if nb(dir,i) is out
 * of lattice.
 * Site OK if ((wait_arr ^ xor_mask ) & and_mask) == 0
 */

static_assert(NDIM <= 4 && 
  "Dimensions at most 4 in dir_mask_t = unsigned char!  Use larger type to circumvent");

void lattice_struct::initialize_wait_arrays()
{
  int i,dir;

  /* Allocate here the mask array needed for forallsites_wait
   * This will contain a bit at location dir if the neighbour
   * at that dir is out of the local volume
   */

  wait_arr_  = (dir_mask_t *)memalloc( this_node.sites * sizeof(unsigned char) );

  for (int i=0; i<this_node.sites; i++) {
    wait_arr_[i] = 0;    /* basic, no wait */
    foralldir(dir) {
      direction odir = -dir;
      if ( neighb[dir][i]  >= this_node.sites ) wait_arr_[i] = wait_arr_[i] | (1<<dir) ;
      if ( neighb[odir][i] >= this_node.sites ) wait_arr_[i] = wait_arr_[i] | (1<<odir) ;
    }
  }
}

#else

void lattice_struct::initialize_wait_arrays(){}

#endif


/////////////////////////////////////////////////////////////////////
/// Create the neighbour index arrays 
/// This is for the index array neighbours
/// TODO: implement some other neighbour schemas!
/////////////////////////////////////////////////////////////////////

#if 1

/// This is a helper routine, returning a vector of comm_node_structs for all nodes
/// involved with communication. 
/// If receive == true, this is "receive" end and index will be filled.
/// For receive == false the is "send" half is done.

std::vector<lattice_struct::comm_node_struct> 
lattice_struct::create_comm_node_vector( coordinate_vector offset, unsigned * index, bool receive) {

  // for send flip the offset
  if (!receive) offset = -offset;

  // temp work array: np = node points
  std::vector<unsigned> np_even(nodes.number);   // std::vector initializes to zero
  std::vector<unsigned> np_odd(nodes.number);

  // we'll go through the sites twice, in order to first resolve the size of the needed
  // buffers, then to fill them.  Trying to avoid memory fragmentation a bit

  // pass over sites
  int num = 0;  // number of sites off node
  for (int i=0; i<this_node.sites; i++) {
    coordinate_vector ln, l;
    l  = coordinates(i);
    ln = mod( l + offset, size() );
 
    if (is_on_node(ln)) {
      if (receive) index[i] = site_index(ln);
    } else {
	    // Now site is off-node, this will leads to fetching
      // use ci.index to store the node rank
      unsigned r = node_rank(ln);

      if (receive) {
        index[i] = this_node.sites + r;

        // using parity of THIS
        if (l.coordinate_parity() == EVEN) np_even[r]++;
        else np_odd[r]++;

      } else {
        // on the sending side - we use parity of target
        if (ln.coordinate_parity() == EVEN) np_even[r]++;  
        else np_odd[r]++;
      }

      num++;
    }
  }

  // count the number of nodes taking part
  unsigned nnodes = 0;
  for (int r=0; r<nodes.number; r++) {
    if (np_even[r] > 0 || np_odd[r] > 0) nnodes++;
  }

  // allocate the vector
  std::vector<comm_node_struct> node_v(nnodes);

  int n = 0;
  int c_buffer = 0;
  for (int r = 0; r<nodes.number; r++) {
    if (np_even[r] > 0 || np_odd[r] > 0) {
      // add the rank
      node_v[n].rank = r;
      node_v[n].evensites = np_even[r];
      node_v[n].oddsites  = np_odd[r];
      node_v[n].sites = np_even[r] + np_odd[r];

      // pre-allocate the sitelist for sufficient size
      if (!receive) node_v[n].sitelist = (unsigned *)memalloc(node_v[n].sites * sizeof(unsigned));

      node_v[n].buffer = c_buffer;        // running idx to comm buffer - used from receive
      c_buffer += node_v[n].sites;
      n++;
    }
  }

  // ci.receive_buf_size = c_buffer;  // total buf size

  // we'll reuse np_even and np_odd as counting arrays below
  for (int i=0; i<nnodes; i++) np_even[i] = np_odd[i] = 0;

  if (!receive) {
    // sending end -- create sitelists

    for (int i=0; i<this_node.sites; i++) {
      coordinate_vector ln, l;
      l  = coordinates(i);
      ln = mod( l + offset, size() );
 
      if (!is_on_node(ln)) {
        unsigned r = node_rank(ln);
        int n = 0;
        // find the node from the list 
        while (node_v[n].rank != r) n++;
        // we'll fill the buffers according to the parity of receieving node
        // first even, then odd sites in the buffer
        int k;
        if (ln.coordinate_parity() == EVEN) k = np_even[n]++;
        else k = node_v[n].evensites + np_odd[n]++;

        // and set the ptr to the site to be communicated
        node_v[n].sitelist[k] = i;        
      }
    }

  } else {
    // receive end
    // fill in the index pointers

    for (int i=0; i<this_node.sites; i++) {
      if (index[i] >= this_node.sites) {
        int r = index[i] - this_node.sites;
        int n = 0;
        // find the node which sends this 
        while (node_v[n].rank != r) n++;

        coordinate_vector l = coordinates(i);
        if (l.coordinate_parity() == EVEN)
          index[i] = node_v[n].buffer + (np_even[n]++);
        else 
          index[i] = node_v[n].buffer + node_v[n].evensites + (np_odd[n]++);
      }
    }
  }

  return node_v;
}




lattice_struct::gen_comminfo_struct lattice_struct::create_general_gather( const coordinate_vector & offset )
{
  // allocate neighbour arrays - TODO: these should 
  // be allocated on "device" memory too!
    
  gen_comminfo_struct ci;

  // communication buffer
  ci.index = (unsigned *)memalloc(this_node.sites*sizeof(unsigned));

  ci.from_node = create_comm_node_vector(offset, ci.index, true);  // create receive end
  ci.to_node   = create_comm_node_vector(offset, nullptr, false);  // create sending end

  // set the total receive buffer size from the last vector
  const comm_node_struct & r = ci.from_node[ci.from_node.size()-1];
  ci.receive_buf_size = r.buffer + r.sites;

  return ci;
}

#endif

//////////////////////////////////////////////////////////////////////////////
/// Test the standard gather here
//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct test_tt {
  T r[NDIM];

  using base_type = typename base_type_struct<T>::type;
};

using test_t = test_tt<int>;

void test_std_gathers()
{

  extern lattice_struct * lattice;
  field<test_t> t;
  field<double> f;
  
  onsites(ALL) {
    foralldir(d) {
      t[X].r[d] = X.coordinates()[d];
    }
  }

  for (parity p : {EVEN,ODD,ALL}) {

    foralldir(d) {
      direction d2;
      for (d2=d; is_up_dir(d2); d2=-d) {
      
        int diff = 0;
        int add;
        if (is_up_dir(d2)) add = 1; else add = -1;
        onsites(p) {
          element<int> i = (t[X].r[d] + add + lattice->size(d)) % lattice->size(d) - t[X+d2].r[d];
          diff += i;

        }

        if (diff != 0) {
          hila::output << "Std gather test error! Node " << mynode() 
                       << " Parity " << parity_name(p) << " direction " << (unsigned)d2 << '\n';
          exit(-1);
        }
      }
    }
  }
}



