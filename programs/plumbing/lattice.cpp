
#include "lattice.h"

///***********************************************************
/// setup() lays out the lattice infrastruct, with neighbour arrays etc.


void lattice_struct::setup(int siz[NDIM]) {
  volume = 1;
  for (int i=0; i<NDIM; i++) {
    size[i] = siz[i];
    volume *= siz[i];
  }
  
  /* Do now the basic lattice layout - array map_node_list contains the
   * non-trivial ordering of the nodes
   */
  #ifdef USE_MPI
  map_node_list = (int *)memalloc(numnodes()*sizeof(int));
  #endif
  setup_layout( siz, nsquares, squaresize, map_node_list );

  /* then, set up the comm arrays */
  make_lattice_arrays( &lattice );

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
/// Is the coordinate on THIS node 
///////////////////////////////////////////////////////////////////////

bool lattice_struct::is_on_node(const location & loc)
{
  int d,dir;

  foralldir(dir) {
    d = loc[dir] - node.xmin[dir];
    if (d < 0 || d >= node.nodesize[dir] ) return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
/// Get the node number for (BLOCKED) coordinates 
///////////////////////////////////////////////////////////////////////

int lattice_struct::node_number(const location & loc)
{
  unsigned i;
  int dir;

  // 
  i = (loc[NDIM-1] * nodes.number_to_dir[NDIM-1]) / size[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    i = i*nodes.number_to_dir[dir] +
        ((loc[dir] * nodes.number_to_dir[dir]) / size[dir]);
  }
#ifdef USE_MPI
  /* do we want to remap this?  YES PLEASE */
  i = remap_node(i);
#endif

  return( i );
}


/************** give site index for ON NODE sites */

int node_index(int loc[NDIM], node_struct *node)
{
  int dir,l,i,s;
  i = l = loc[NDIM-1] - node->xmin[NDIM-1];
  s = loc[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    l = loc[dir] - node->xmin[dir];
    i = i*node->nodesize[dir] + l;
    s += loc[dir];
  }


#if USE_MORTON_CODE || LOCAL_SITES_FIRST
  //static int mapdone = 0;
  //static int* map = NULL;
  if (!node->map){
    int nodesize[NDIM];
    foralldir(dir) nodesize[dir] = node->nodesize[dir];
#if defined(EVENFIRST) && !LOCAL_SITES_FIRST
    for (dir = NDIM-1; dir >= 0; dir--) {
      if ((nodesize[dir] % 2) == 0){
        nodesize[dir] /= 2;
        break;
      }
    }
#endif
    node->map = domap(nodesize, node);
    //mapdone = 1;
  }

#if LOCAL_SITES_FIRST
    #if 0
      if (node->checkIndex){
          printf("Site [ %2d", loc[0]);
          for (dir = 1; dir < NDIM; dir++) printf(", %2d", loc[dir]);
          printf(" ] = %6d = %d\n", i, node->map[i]);
      }
    #endif
  i = node->map[i];
#else // LOCAL_SITES_FIRST
    #if defined(EVENFIRST)
      i = 2*node->map[i/2];
    #else
      i = node->map[i];
    #endif
#endif  // LOCAL_SITES_FIRST
#endif // USE_MORTON_CODE || LOCAL_SITES_FIRST
  /* now i contains the `running index' for site */
// NOTE: LOCAL_SITES_FIRST compute explicitly the even-odd mapping!
#if defined(EVENFIRST) && !(LOCAL_SITES_FIRST)
  if (s%2 == 0) return( i/2 );    /* even site index */
  else return( i/2 + node->evensites );  /* odd site */
#else
  return( i );
#endif
}



void lattice_struct::loop_begin( parity p )
{
  #ifndef USE_MPI
    #if defined(EVENFIRST)
      if(p==ODD){
        return mynode->evensites;
      } else {
        return 0; // EVEN or ALL
      }
    #endif
  #else //USE_MPI
    not implemented
  #endif
}


void lattice_struct::loop_begin( parity p )
{
  #ifndef USE_MPI
    #if defined(EVENFIRST)
      if(p==EVEN){
        return mynode->evensites; 
      } else {
        return mynode->sites; // ODD or ALL
      }
    #endif
  #else //USE_MPI
    not implemented
  #endif
}
