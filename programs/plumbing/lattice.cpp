
#include "lattice.h"

///***********************************************************
/// setup() lays out the lattice infrastruct, with neighbour arrays etc.


void lattice_struct::setup(int siz[NDIM]) {
  l_volume = 1;
  for (int i=0; i<NDIM; i++) {
    l_size[i] = siz[i];
    l_volume *= siz[i];
  }
  
  this_node.index = mynode();
  nodes.number = numnodes();
  
  setup_layout();
  this_node.setup();

  /* then, set up the comm arrays */
  setup_communications();
      
#ifdef USE_MPI
  /* Initialize wait_array structures */
  initialize_wait_arrays();
#endif

}


#if NDIM==4
void lattice_struct::setup(int nx, int ny, int nz, int nt) {
  int s[NDIM] = {nx, ny, nz, nt};
  setup(s);
}
#elif NDIM==3
void lattice_struct::setup(int nx, int ny, int nz) {
  int s[NDIM] = {nx, ny, nz};
  setup(s);
}
#elif NDIM==2
void lattice_struct::setup(int nx, int ny) {
  int s[NDIM] = {nx, ny};
  setup(s);
}
#elif NDIM==1
void lattice_struct::setup(int nx) {
  int s[NDIM] = {nx};
  setup(s);
}
#endif

///////////////////////////////////////////////////////////////////////
/// The routines is_on_node(), node_number(), site_index()
/// implement the "layout" of the nodes and sites of the lattice.
/// To be changed in different implementations!
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
/// Get the node number for coordinates 
/// This is the fundamental routine which defines how the nodes
/// are mapped.  map_node_layout needs to be compatible with this
/// algorithm!
///////////////////////////////////////////////////////////////////////

unsigned lattice::node_number(const location & loc)
{
  unsigned i;
  int dir;

  i = (loc[NDIM-1] * nodes.ndir[NDIM-1]) / l_size[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    i = i*nodes.ndir[dir] +
        ((loc[dir] * nodes.ndir[dir]) / l_size[dir]);
  }
  /* do we want to remap this?  YES PLEASE */
  i = nodes.remap(i);

  return( i );
}

///////////////////////////////////////////////////////////////////////
/// Is the coordinate on THIS node 
///////////////////////////////////////////////////////////////////////

bool lattice::is_on_node(const location & loc)
{
  int d;

  for (direction dir=0; dir<NDIM; dir++) {
    d = loc[dir] - this_node.xmin[dir];
    if (d < 0 || d >= this_node.size[dir] ) return false;
  }
  return true;
}


///////////////////////////////////////////////////////////////////////
/// give site index for ON NODE sites
/// Note: loc really has to be on this node
///////////////////////////////////////////////////////////////////////

unsigned lattice::site_index(const location & loc)
{
  int dir,l,s;
  unsigned i;

  i = l = loc[NDIM-1] - this_node.xmin[NDIM-1];
  s = loc[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    l = loc[dir] - this_node.xmin[dir];
    i = i*this_node.size[dir] + l;
    s += loc[dir];
  }

  // now i contains the `running index' for site
#if defined(EVENFIRST)
  if (s%2 == 0) return( i/2 );    /* even site index */
  else return( i/2 + this_node.evensites );  /* odd site */
#else
  return( i );
#endif
}

///////////////////////////////////////////////////////////////////////
/// give site index for nodeid sites
/// Note: loc really has to be on this node
///////////////////////////////////////////////////////////////////////

unsigned lattice::site_index(const location & loc, const unsigned nodeid)
{
  int dir,l,s;
  unsigned i;
  node_info & ni = nodes.nodelist[nodeid];
  
  i = l = loc[NDIM-1] - ni.xmin[NDIM-1];
  s = loc[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    l = loc[dir] - ni.xmin[dir];
    i = i*ni.size[dir] + l;
    s += loc[dir];
  }

  // now i contains the `running index' for site
#if defined(EVENFIRST)
  if (s%2 == 0) return( i/2 );    /* even site index */
  else return( i/2 + ni.evensites );  /* odd site */
#else
  return( i );
#endif
}


///////////////////////////////////////////////////////////////////////
/// invert the this_node index -> location
///////////////////////////////////////////////////////////////////////

location lattice::site_location(unsigned index)
{
  // make the index lexicographic
#ifdef EVENFIRST
  if (index < this_node.evensites) {
    index *= 2;
    if (!this_node.first_site_even) index++;
  } else {
    index = (index - this_node.evensites)*2;
    if (this_node.first_site_even) index++;
  }
#endif
  
  location l;
  foralldir(d) {
    l[d] = index % this_node.size[d] + this_node.min[d];
    index /= this_node.size[d];
  }
  
  return l;
}






/////////////////////////////////////////////////////////////////////
/// Set up the basic list about all nodes
/// NOTE: in case of very large number of nodes this may be
/// inconveniently large.  Alternatives?
/////////////////////////////////////////////////////////////////////

void lattice::setup_nodes() {

  nodes.number = numnodes();
  // Loop over all node "origins"
  nodes.nodelist.resize(nodes.number);

  int n[NDIM];
  foralldir(d) n[d] = 0;

  // use nodes.divisors - vectors to fill in stuff
  for (int i=0; i<nodes.number; i++) {
    location l;
    foralldir(d) l[d] = nodes.divisors[n[d]];

    int nn = node_number(l);
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
      if (location_parity(l) == EVEN) {
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
      if (n[d] < nodes.ndir[d]) break;
      n[d] = 0;
    }
    
    // use the opportunity to set up this_node when it is met
    if (nn == mynode()) this_node.setup(ni);
    
  }
}


////////////////////////////////////////////////////////////////////////
/// Fill in this_node fields -- node_number() must be set up OK
////////////////////////////////////////////////////////////////////////

void lattice::this_node.setup(node_info & ni)
{

  this_node.index = mynode();

  foralldir(d) {
    this_node.min[d]  = ni.min[d];
    this_node.size[d] = ni.size[d];
  }
  this_node.evensites = ni.evensites;
  this_node.oddsites  = ni.oddsites;
  this_node.sites     = ni.evensites + ni.oddsites;

  int i = 0;
  foralldir(d) i += this_node.min[d];
  if (i % 2 == 0) this_node.first_site_even = true;
  else            this_node.first_site_even = false;
   
  // neighbour node indices
  foralldir(d) {
    location l = this_node.min;
    l[d] = (this_node.min[d] + this_node.size[d]) % l_size[d];
    this_node.nn[d] = node_number(l);
    l[d] = (l_size[d] + this_node.min[d] - 1) % l_size[d];
    this_node.nn[opp_dir(d)] = node_number(l);
  }
}


/////////////////////////////////////////////////////////////////////
/// Create the neighbour index arrays 
/// This is for the index array neighbours
/// TODO: implement some other neighbour schemas!
/////////////////////////////////////////////////////////////////////

void lattice::create_std_gathers()
{
  // allocate neighbour arrays - TODO: these should 
  // be allocated on "device" memory too!
  
  for (direction d=0; d<NDIRS; d++) {
    neigb[d] = (unsigned *)allocate_field_mem(this_node.sites * sizeof(unsigned));
  }
  
  comminfo.resize(MAX_GATHERS);
  for (auto & c : comminfo) c.label = -1;  // reset
  unsigned c_offset = this_node.sites;  // current offset in field-arrays

  // allocate work arrays, will be released later
  std::vector<unsigned> nnodes(this_node.sites); // node number
  std::vector<unsigned> index(this_node.sites);  // index on node
  std::vector<unsigned> parity(this_node.sites); // parity 
  std::vector<unsigned> here(this_node.sites);   // index of original site
  std::vector<unsigned> itmp(this_node.sites);   // temp sorting array

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
  if (this_node.xmin[NDIM-1] + this_node.nodesize[NDIM-1] == size[NDIM-1])
    c_offset += 1;
  else sf_special_boundary = -(1<<30);

  output0 << "SPECIAL BOUNDARY LAYOUT for SF");
#endif

  for (direction d=0; d<NDIRS; d++) {

    comminfo_struct & ci = comminfo[d];
    ci.index = neighb[d];  // the neighbour array - not needed really
    ci.label = d;          // and label it 

    // pass over sites
    int num = 0;  // number of sites off node
    for (int i=0; i<this_node.sites; i++) {
      location ln,l = site_location(i);
      // set ln to be the neighbour of the site
      if (is_up_dir(d)) {
        ln[d] = (l[d] + 1) % size[d];
      } else {
        direction k = opp_dir(d);
        ln[k] = (l[k] + size[k] - 1) % size[k];
      }
 
#ifdef SCHROED_FUN
      if (d == NDIM-1 && l[NDIM-1] == size[NDIM-1]-1) {
	// This is up-direction, give special site
	neighb[d][i] = sf_special_boundary;
      } else if (d == opp_dir(NDIM-1) && l[NDIM-1] == 0) {
	// We never should need the down direction, thus, short circuit!
	neighb[d][i] = 1<<30;
      } else     // NOTE THIS UGLY else-if CONSTRUCT!
#endif
      if (is_on_node(ln)) {
        neighb[d][i] = site_index(ln);
      } else {
	// Now site is off-node, this lead to fetching
	nodes[num] = node_number(ln);
	index[num] = node_index(x, allnodes + nodes[num] );
	parity[num] = location_parity(l);  // parity of THIS
	here[num]  = i;
	num++;
      }
    }

    // construct nodes to be gathered from
    ci.from_node = {};    

    if (num > 0) {
      // now, get the number of nodes to be gathered from
      for (int i=0; i<num; i++) {
	// chase the list until the node found
        for (int j=0; j<ci.from_node.size() && nodes[i] != ci.from_node[j].node; j++);
	if (j == ci.from_node.size()) {
	  // NEW NODE to receive from
          comm_node_struct s;
          s.node = nodes[i];
          s.n = s.n_even = s.n_odd = 0;
          ci.from_node.push_back(s);
	}
	// add to OLD NODE or just added node
        ci.from_node[j].n++;
        if ( parity[i] == EVEN ) ci.from_node[j].n_even ++;  
        else ci.from_node[j].n_odd++;
      }

      // Calculate the buffer offsets for the gathers
      for (comm_node_struct & fn : ci.from_node) {
        fn.buffer = c_offset;
        c_offset += fn.n;  // and increase the offset
      }

      // and NOW, finish the NEIGHBOR array 

      for (comm_node_struct & fn : ci.from_node) {
        // Now, accumulate the locations to itmp-array, and sort the
        // array according to the index of the sending node .
        // First even neighbours

	for (parity par=EVEN; par<=ODD; par++) {
          int n,i;
	  for (n=i=0; i<num; i++) if (nodes[i] == fn.node && parity[i] == par) {
	    itmp[n++] = i;
	    // bubble sort the tmp-array according to the index on the neighbour
            // node
	    for (int k=n-1; k > 0 && index[itmp[k]] < index[itmp[k-1]]; k--)
	      swap( itmp[k], itmp[k-1] );
	  }
	  off = fn.buffer;
	  if (par == ODD) off += fn.n_even;
	  // finally, root indices according to offset
	  for (int k=0; k<n; k++) neighb[d][here[itmp[k]]] = off + k;
	}
      }
    } // num > 0 

#if 0
    //Check that we are not dividing odd lattices, breaks communication
    if (num > n_sublattices) if((nx%2==1)||(ny%2==1)||(nz%2==1)||(nt%2==1))
    halt("Odd number of sites cannot be devided into nodes.");
#endif

    // receive done, now opposite send. This is just the gather
    // inverted

    od = opp_dir(d);
    comminfo[od].n_send = comminfo[d].n_receive;

    if (num > 0) {
      comm_node_struct & fn = comminfo[d].from_node;
      for (j=0, s=&(comlist[od].to_node); j<comlist[od].n_send;
	   j++, s = &((*s)->next), p = p->next) {
	(*s) = q = (send_struct *)memalloc( sizeof(send_struct) );
	q->node   = p->node;
	q->n      = p->n;
	q->n_even = p->n_odd;     /* Note the swap !  even/odd refers to type of gather */
	q->n_odd  = p->n_even;
	q->next   = NULL;
	q->sitelist = (int *)memalloc(q->n * sizeof(int));

	/* now, initialize sitelist -- Now, we first want ODD parity, since
	 * this is what even gather asks for!
	 */

	for (n=0,par=ODD; par>=EVEN; par--) {
	  for (i=0; i<num; i++) if (nodes[i] == q->node && parity[i] == par) {
	    (q->sitelist)[n++] = here[i];
	  }
	  if (par == ODD && n != q->n_even) halt("Parity odd error 3");
	  if (par == EVEN && n != q->n) halt("Parity even error 3");
	}

	// ignore return value
	(void)parallel_initDevSitelist(q);

      }
    }
  } /* directions */

  free(nodes);
  free(index);
  free(parity);
  free(here);
  free(itmp);

  /* Finally, set the site to the final offset (better be right!) */
  node.latfield_size = c_offset;

}

