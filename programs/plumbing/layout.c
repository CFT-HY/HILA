/****************************************************************
 *                                                              *
 *    Hypercubic lattice layout routines                        *
 *    Based on MILC lattice QCD code, pretty much modified      *
 *                                                              *
 *    These determine the distribution of sites on nodes,       *
 *    and do the necessary setting up.                          *
 *                                                              *
 ***************************************************************/
#ifdef ALIGN
#define _XOPEN_SOURCE 600
#endif


#include "../include/comdefs.h"
#include "../include/generic.h"
#include "../wilson/lattice.h"
#include "../parallellib/parallel_defs.h"
#include "../parallellib/parallel_util.h"

//Allow for odd number of sites when not using LU
//This causes problems if trying to devide the lattice to nodes
#ifndef LU
#define IGNORE_PARITY
#endif

/* static variables for node calculations */
int squaresize[NDIM];	/* dimensions of hypercubes */
int nsquares[NDIM];	/* number of hypercubes in each direction */

/* GLOBALS for communications; needed by com_XXX.c and block_lattice.c */
node_struct *allnodes;          /* structure for all nodes on this run */
comlist_struct *comlist;        /* gather pointer for all gathers */

static int* map_node_list = NULL;

#define swap(a,b) {register int t; t=a; a=b; b=t; }

void test_std_gathers( lattice_struct *lat );
void make_std_gathers( lattice_struct *lat );
void test_gather(int gat,
		 void (*fmap)(int x[NDIM], int *arg, int mode, int xp[NDIM]),
		 int *arg, int mode);


static int coordinate_parity(int x[NDIM]);

/***************************************************************
 *  BASIC CALL FOR SETUP
 *
 *  setup_lattice(int size[NDIM]);
 */

void setup_lattice(int siz[NDIM])
{
  int dir;

  /* static global */
  this_node = mynode();

  /* store the baseline */
  base_lattice = lattice;

  /* reset the blocking level (just in case) */
  foralldir(dir) current_blocking_level[dir] = 0;

  /* set lattice volume and dimensions */
  lattice.volume = 1;
  foralldir(dir) {
    lattice.size[dir] = siz[dir];
    lattice.volume *= lattice.size[dir];
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

/**************** Get the node number for (BLOCKED) coordinates */
int node_number(int loc[NDIM])
{
  register int i,dir;

  i = (loc[NDIM-1] << current_blocking_level[NDIM-1]) / squaresize[NDIM-1];
  for (dir=NDIM-2; dir>=0; dir--) {
    i = i*nsquares[dir] +
      ((loc[dir] << current_blocking_level[dir]) / squaresize[dir]);
  }
#ifdef USE_MPI
  /* do we want to remap this?  YES PLEASE */
  i = map_node_list[i];
#endif

#ifdef GRAYCODE
  return( i ^ (i>>1) );	/* Gray code of i */
#else
  return( i );
#endif
}

/************** fast routine for clarifying if we're on THIS node */

int is_on_node(int loc[NDIM])
{
  register int d,dir;

  foralldir(dir) {
    d = loc[dir] - node.xmin[dir];
    if (d < 0 || d >= node.nodesize[dir] ) return(0);
  }
  return(1);
}

#define USE_MORTON_CODE    0
#if USE_MORTON_CODE

static int complog2(unsigned int x)
{
  int r;
  int shift;
  // get log2 - http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
  r =     (x > 0xFFFF) << 4; x >>= r;
  shift = (x > 0xFF  ) << 3; x >>= shift; r |= shift;
  shift = (x > 0xF   ) << 2; x >>= shift; r |= shift;
  shift = (x > 0x3   ) << 1; x >>= shift; r |= shift;
  r |= (x >> 1);
  return r;
}

static int computeMorton(int* xs, int* sizes){
 // Iterate over every log2 - size (ie. level)
 int level = 0;
 int dir;
 int sizeslog2[NDIM];
 int maxlog2 = 0;
 int bitIndex = 0;
 unsigned int code = 0; // This will contain the resulting code
 // First find out max dimension (to compute loop-limit)
 foralldir(dir){
   sizeslog2[dir] = complog2(sizes[dir]);
   if (sizeslog2[dir] > maxlog2) maxlog2 = sizeslog2[dir];
 }

 // Ok now we can start iterating over all log2-sizes
 for (level = 0; level < maxlog2; level++)
 {
   // And iterate over every direction
   foralldir(dir)
   {
     // Only add a bit if this dimension is larger or equal than the current level-size
     if (level < sizeslog2[dir])
     {
       int coordinate = xs[dir];
       // extract bit
       int result = (coordinate & (1 << level)) >> level;
       // And put it to its place in the Morton-code
       code |= (result << bitIndex);
       // Finally increment the bitIndex to the next position:
       bitIndex++;
     }
   }
 }
 return code;
}

static int* domap_Mortoncode(int* result, int* nodesize){
  int* resptr = result;
  int xs[NDIM] = { 0 };
  int done = 0;
  while (!done)
  {
    int dir;
    // First compute Morton code for current coordinates:
    *resptr++ = computeMorton(xs, nodesize);
    // Then inc coordinate
    for (dir = 0; dir < NDIM; dir++)
    {
      xs[dir]++;
      if (xs[dir] < nodesize[dir]) break;
      // else
      xs[dir] = 0;
      if (dir == NDIM - 1)
        done = 1;
    }
  }
  return result;
}


static int* domap(int* nodesize, node_struct* node)
{
  int* result;
  size_t size = 1;
  int dir;
  int pow2 = 1;
  foralldir(dir) size *= nodesize[dir];
  result = (int*)malloc(size * sizeof(int));
  // Check if every direction is of power-of-two size
  foralldir(dir){
    unsigned int dimsize = (unsigned int)nodesize[dir];
    // Check power of two - simple to check: If x is power of two, then
    // it has only one bit up, therefore if we subtract one from x then
    // for power-of-two numbers it will have all bits filled from zero to
    // one less than the position of the highest bit of the original x
    // Therefore if x is pow2, then x = 0010 0000 and (x-1) = 0001 1111
    // No two bits intersect. On the other hand if x is not pow2, then
    // the highest bit of x will not flip in x-1 and at least this bit
    // will intersect with x and (x-1)
    // Therefore x is pow 2 iff. x > 0 and (x & (x-1)) == 0,
    // (if 0 is not considered 2^-infinity)
    if (dimsize == 0 || (dimsize & (dimsize - 1)) != 0){
      pow2 = 0;
      break;
    }
  }
  if (pow2) return domap_Mortoncode(result, nodesize);
  // Otherwise do a simple blocking-theme. This should help with normal caches, but
  // might not help as much with different levels of cache-structure
  #define TMP_BLOCK_SIZE    4
  // Sites go like: (0,0,0,0), (0,0,0,1), (0,0,0,2), ... (0,0,1,0), and so on
  {
    int x0s[NDIM] = { 0 };
    int xs[NDIM] = { 0 };
    int limits[NDIM];
    int done = 0;
    int i = 0;
    while (!done)
    {
      int subdone = 0;
      foralldir(dir) limits[dir] =
          TMP_BLOCK_SIZE + x0s[dir] < nodesize[dir] ? TMP_BLOCK_SIZE + x0s[dir] : nodesize[dir];
      foralldir(dir) xs[dir] = x0s[dir];
      while (!subdone)
      {
        int index = 0;
        int subsize = 1;
        foralldir(dir) { index += xs[dir] * subsize; subsize *= nodesize[dir]; }
        result[i++] = index;
        // inc indices and move forwards - TODO:Check logic
        foralldir(dir) { xs[dir]++; if (xs[dir] < limits[dir]) break; xs[dir] = x0s[dir]; if (dir == NDIM-1) subdone = 1; }
      }
      foralldir(dir) {
        x0s[dir] = limits[dir];
        if (x0s[dir] < nodesize[dir]) break;
        // else
        x0s[dir] = 0;
        if (dir == NDIM-1) done = 1;
      }
    }
  }
  return result;
  #undef BLOCK_SIZE
}
#endif



#if LOCAL_SITES_FIRST

// Check whether any nearest neighbor site is across node borders
static int is_on_node_border(int coords[NDIM])
{
    int dir;
    for (dir = 0; dir < NDIM; dir++){
        int x0 = coords[dir];
        int x1 = x0 - 1;
        int x2 = x0 + 1;
#ifdef SCHROED_FUN
        if (dir == TUP){
            if (x0 == lattice.size[dir] - 1)
                x2 = x0;
            if (x0 == 0)
                x1 = x0;
        }
#endif
        if (x1 < 0) x1 = lattice.size[dir] - 1;
        if (x2 >= lattice.size[dir]) x2 = 0;
        if (is_across_node_boundary(x0, x1, dir))
            return 1;
        if (is_across_node_boundary(x0, x2, dir))
            return 1;
    }
    return 0;
}

static int* domap(int* nodesize, node_struct* node)
{
    int* result;
    int* types;
    size_t size = 1;
    int dir;

    // Note - we don't know yet how many sites there are in the inner volume so
    // first write them as negative numbers and after passing the whole thing, adjust them
    int evenNInner = 0;
    int oddNInner = 0;
    int evenNBorder = 0;
    int oddNBorder = 0;
    int nIndices = 0;

    foralldir(dir) size *= nodesize[dir];
    result = (int*)malloc(size * sizeof(int));
    types = (int*)malloc(size * sizeof(int));
    if (!result || !types)
        halt("Malloc failure in layout.c:domap()");

    int xs[NDIM] = { 0 };
    int done = 0;
    while (!done)
    {
      int dir;
      int tmpxs[NDIM];
      int parity = EVEN;
      for (dir = 0; dir < NDIM; dir++)
          tmpxs[dir] = xs[dir] +  node->xmin[dir];
#ifdef EVEN_SITES_FIRST
      parity = coordinate_parity(tmpxs);
#endif

      if (is_on_node_border(tmpxs)){
          int index;
          if (parity == EVEN){
              index = evenNBorder++;
              types[nIndices] = 1;
          } else {
              index = oddNBorder++;
              types[nIndices] = 2;
          }
          result[nIndices] = index;
      } else {
          int index;
          if (parity == EVEN){
              index = evenNInner++;
              types[nIndices] = 3;
          } else {
              index = oddNInner++;
              types[nIndices] = 4;
          }
          result[nIndices] = index;
      }
      // Then inc coordinate
      for (dir = 0; dir < NDIM; dir++)
      {
        xs[dir]++;
        if (xs[dir] < nodesize[dir]) break;
        // else
        xs[dir] = 0;
        if (dir == NDIM - 1)
          done = 1;
      }
      nIndices++;
    }
    // Finally walk through the list and put in the correct offsets according to type:
    {
        int i;
        for (i = 0; i < nIndices; i++){
            if (types[i] == 1) // Even border
                result[i] += evenNInner;
            else if (types[i] == 2) // Odd border
                result[i] += node->evensites + oddNInner;
            else if (types[i] == 4) // Odd inner
                result[i] += node->evensites;
        }
    }
    node->nEvenInner = evenNInner;
    node->nEvenBorder = evenNBorder;
    node->nOddInner = oddNInner;
    node->nOddBorder = oddNBorder;
#if 0
    node->checkIndex = 0;
#endif
    free(types);
    return result;
}
#endif


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
#if defined(EVEN_SITES_FIRST) && !LOCAL_SITES_FIRST
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
    #if defined(EVEN_SITES_FIRST)
      i = 2*node->map[i/2];
    #else
      i = node->map[i];
    #endif
#endif  // LOCAL_SITES_FIRST
#endif // USE_MORTON_CODE || LOCAL_SITES_FIRST
  /* now i contains the `running index' for site */
// NOTE: LOCAL_SITES_FIRST compute explicitly the even-odd mapping!
#if defined(EVEN_SITES_FIRST) && !(LOCAL_SITES_FIRST)
  if (s%2 == 0) return( i/2 );    /* even site index */
  else return( i/2 + node->evensites );  /* odd site */
#else
  return( i );
#endif
}


/*   Helper routine for determining if
 *   given coordinate pair to direction dir
 *   is across node boudary
 *
 *   For safety, we're using the node_number routine,
 *   even though this could be done slightly quicker
 *
 *   ASSUMPTION: node boundaries cut through the volume,
 *   thus, this can be done by setting other coords = 0
 */

int is_across_node_boundary( int x1, int x2, int dir )
{
  int d[NDIM],i,n1,n2;

  if (dir < 0 || dir >= NDIM ||
      x1 < 0 || x1 > lattice.size[dir] ||
      x2 < 0 || x2 > lattice.size[dir] ) {
    fprintf0(outf," Invalid coordinates to is_across_node_boundary\n");
    fprintf0(outf," x1 %d  x2 %d  dir %d\n",x1,x2,dir);
    finishrun();
  }

  for (i=0; i<NDIM; i++) d[i] = 0;
  d[dir] = x1;
  n1 = node_number( d );
  d[dir] = x2;
  n2 = node_number( d );
  return( n1 != n2 );
}

/************** give parity of site */

int coordinate_parity(int x[NDIM])
{
  int d,p;

  p = 0;
  foralldir(d) p += x[d];
  if (p % 2 == 0) return(EVEN);
  else return(ODD);
}


/******************************************************
 * routines for stepping through the lattice in
 * coordinates, as in
 * #define forallcoordinates(x) \
 *  for(zero_arr(x); is_coord(x,&lattice); step_coord(x,&lattice) )
 */

void zero_arr(int x[NDIM]) { register int d; foralldir(d) x[d] = 0; }

int is_allowed_coord(int x[NDIM],lattice_struct *l)
{
  int d,i;
  i = 1;
  foralldir(d) i = (i && (x[d] >= 0) && (x[d] < l->size[d]));
  return(i);
}

void step_coord(int x[NDIM],lattice_struct *l)
{
  int d;

  for(d=0; d<NDIM && ++(x[d]) >= l->size[d]; x[d++] = 0) ;

  /* check if the lattice is 'full' */
  if (d >= NDIM) x[NDIM-1] = l->size[NDIM-1];
}


/**************
 * set up the node structure for all of the nodes in
 * the run (for BLOCKED coordinates).
 */

void setup_node( int loc[NDIM], node_struct *n )
{
  register int offset,dir,blev,l,c0,c1,s;

  n->sites = 1;
#if USE_MORTON_CODE || LOCAL_SITES_FIRST
  n->map = (int*)0;
#endif
  s = 0;
  foralldir(dir) {
    blev = 1 << current_blocking_level[dir];
    l = loc[dir] << current_blocking_level[dir];  /* normalized coord */
    offset = l % squaresize[dir];         /* normalized coord from node 'origin' */
    c0 = l - offset;                      /* coordinate of the origin */
    c1 = c0 + squaresize[dir] - 1;        /* coordinate of the last point */

    /* calculate the coordinate of the first blocked point on the
     * node.  If the origin is divisible by the blocking factor,
     * then it belongs to the blocked lattice and the coordinate is
     * just l0 / blev.  However, if not, then the first blocked point
     * is l0 / blev + 1.
     */
    if (c0 % blev == 0) n->xmin[dir] = c0 >> current_blocking_level[dir];
    else n->xmin[dir] = (c0 >> current_blocking_level[dir]) + 1;

    /* Now the coordinate of the last blocked point.  This is
     * always c1/blev, regardless if it is divisible or not.
     */

    c1 = c1 >> current_blocking_level[dir];

    /* now the length of the blocked lattice */
    n->nodesize[dir] = c1 - n->xmin[dir] + 1;

    /* need to accumulate size */
    n->sites *= n->nodesize[dir];
    /* and parity of the origin */
    s += n->xmin[dir];
  }

  if ( n->sites % 2 ) {
    /* now odd sized node */
    if ( s % 2 == 0) n->evensites = n->sites/2 + 1;
    else n->evensites = n->sites/2;
    n->oddsites = n->sites - n->evensites;
  } else {
    n->evensites = n->oddsites = n->sites/2;
  }
}

/************************************************************
 * set up the node struct for all nodes
 */

node_struct * setup_nodes(lattice_struct *lat)
{
  int i,l,d,n,x[NDIM];
  node_struct *p;

  /* allocate the node array */
  p = (node_struct *)memalloc( (l=numnodes()) * sizeof(node_struct) );
  for (i=0; i<l; i++) p[i].sites = -1;   /* flag as not initialized */

  forallcoordinates(x) {
    /* now x contains the coordinates of the site */
    n = node_number(x);  /* which node? */
    /* set up the node structure, if not done before */
    if (p[n].sites < 0) setup_node(x,p+n);
  }

  for (i=0; i<l; i++) {
    /* not initialized, reset */
    if (p[i].sites == -1) {
      p[i].sites = p[i].evensites = p[i].oddsites = 0;
      foralldir(d) {
	p[i].xmin[d] = p[i].nodesize[0] = 0;
      }
    }

    /* find now the neighbour node indices */
    if (p[i].sites) {
      int xn[NDIM];
      /* bottom corner coordinate of node */
      foralldir(d) xn[d] = x[d] = p[i].xmin[d];
      foralldir(d) {
	/* directions down */
	xn[d] = (x[d] - 1 + lattice.size[d]) % lattice.size[d];
	p[i].down_node[d] = node_number(xn);
	xn[d] = x[d];
      }
      /* upper corner coordinate of node */
      foralldir(d) xn[d] = x[d] = p[i].xmin[d] + p[i].nodesize[d] - 1;
      foralldir(d) {
	/* directions up */
	xn[d] = (x[d] + 1) % lattice.size[d];
	p[i].up_node[d] = node_number(xn);;
	xn[d] = x[d];
      }
    }
  }
  return (p);
}

/************************************************************/


void make_lattice_arrays(lattice_struct * l)
{
  int x[NDIM],i,d,j;

  /* First, set up the node structure */

  allnodes = setup_nodes(l);
  node = allnodes[ this_node ];
  sites = node.sites;

  /* Setup the SITE array */

  site = (site_struct *)memalloc((node.sites+3) * sizeof(site_struct));

  i = 0;  /* index to sites */
  forallcoordinates(x) {
    /* now x contains the coordinates of the site */

    if (is_on_node(x)) {
      j = node_index(x,&node);    /* array index to site */

      /* set the site */
      foralldir(d) site[j].x[d] = x[d];
      site[j].parity = coordinate_parity(x);
      site[j].index = i;
    } /* else fprintf(outf,"Offsite!\n"); */
    i++;
  }

  /* and then, set up the neighbour arrays and bookkeeping */

  make_std_gathers(l);

#ifdef SCHROED_FUN
  special_sf_site = getSpecialSiteIdx();
#endif // SCHROED_FUN

#if 0
      int tmpnodei = 0;
      for (tmpnodei = 0; tmpnodei < numnodes(); tmpnodei++){
          fflush(stdout);
          g_sync();
          if (0 == mynode() && tmpnodei == 0){
              printf("\n\n Node = %d\n\n", tmpnodei);
              node.checkIndex = 1;
              forallsites(i){
                  j = node_index(site[i].x, &node);
              }
              node.checkIndex = 0;
          }
          fflush(stdout);
          g_sync();
      }
#endif

#if LOCAL_SITES_FIRST
  {
      int i, dir;
      int border = 0;
#ifdef EVEN_SITES_FIRST
      forevensites(i){
          int isOnBorder = is_on_node_border(site[i].x);
          if (isOnBorder)
              border = 1;
          else if (border == 1)
              halt("Found an inner-volume site after border-sites ");
          if (!border){
              for(dir = 0; dir < NDIRS; dir++)
                  if_sf_site(i)
                  if (!h_is_on_node(nb(dir,i))){
                      printf("i, dir = (%d, %d), nb(dir, i) = %d\n", i, dir, nb(dir,i));
                      fflush(stdout);
                      halt("Found an inner-volume site on the border!");
                  }
          }
          else {
              int out = 0;
              for(dir = 0; dir < NDIRS; dir++)
                  if (!h_is_on_node(nb(dir,i)))
                      out = 1;
              if (!out){
                  int dir2;
                  printf("i = %d\n", i);
                  printf("FAULT = [ %d", site[i].x[0]);
                  for (dir2 = 1; dir2 < NDIM; dir2++){
                      printf(", %d", site[i].x[dir2]);
                  }
                  printf("\n");
                  fflush(stdout);
                  halt("Found a border site in the inner-volume!");
              }
          }
      }
      border = 0;
      foroddsites(i){
          int isOnBorder = is_on_node_border(site[i].x);
          if (isOnBorder)
              border = 1;
          else if (border == 1)
              halt("Found an inner-volume site after border-sites ");
          if (!border){
              for(dir = 0; dir < NDIRS; dir++)
                  if_sf_site(i)
                  if (!h_is_on_node(nb(dir,i))){
                      printf("i, dir = (%d, %d), nb(dir, i) = %d\n", i, dir, nb(dir,i));
                      fflush(stdout);
                      halt("Found an inner-volume site on the border!");
                  }
          }
          else {
              int out = 0;
              for(dir = 0; dir < NDIRS; dir++)
                  if (!h_is_on_node(nb(dir,i)))
                      out = 1;
              if (!out)
                  halt("Found a border site in the inner-volume!");
          }
      }
#else // EVEN_SITES_FIRST
      forallsites(i){
          int isOnBorder = is_on_node_border(site[i].x);
          if (isOnBorder)
              border = 1;
          else if (border == 1)
              halt("Found an inner-volume site after border-sites ");
          if (!border){
              for(dir = 0; dir < NDIRS; dir++)
                  if (!h_is_on_node(nb(dir,i)))
                      halt("Found an inner-volume site on the border!");
          }
          else {
              int out = 0;
              for(dir = 0; dir < NDIRS; dir++)
                  if (!h_is_on_node(nb(dir,i)))
                      out = 1;
              if (!out)
                  halt("Found a border site in the inner-volume!");
          }
      }
#endif // EVEN_SITES_FIRST
  }
#endif // LOCAL_SITES_FIRST
  /* and test it - Note: This test has been moved back, since CUDA-arrays have not
   * been allocated yet. Therefore this test is run after alloc-fields from control.c */
  //test_std_gathers(l);

}

/***************************************************************************
 * A long routine to make the gather datastructs (comlists)
 * and neighbour arrays
 */

void make_std_gathers( lattice_struct *lat )
{
  int d,i,j,k,c_offset,x[NDIM];
  int *nodes, *index, *parity, *here, *itmp, *tmpbuf, num;
  int par,n,off,od;
  receive_struct **r,*p;
  send_struct **s,*q;

  /* First, allocate neighbour arrays */
  foralldir(d) {
    neighb[d]          = (int *)memalloc(node.sites * sizeof(int));
    neighb[opp_dir(d)] = (int *)memalloc(node.sites * sizeof(int));
  }

  /* allocate comlist structure - needed as ptr for blocking! */
  comlist = (comlist_struct *)memalloc(MAX_GATHERS * sizeof(comlist_struct));
  for (i=0; i<MAX_GATHERS; i++) comlist[i].label = -1;  /* reset */

  /* work arrays, will be released later */
  /* allocate first itmp-array, release it after allocation,
     in order to ensure low memory */
  tmpbuf = (int *)memalloc(3*node.sites * sizeof(int));   /* temporary space */
  nodes  = (int *)memalloc(node.sites * sizeof(int));  /* node number */
  index  = (int *)memalloc(node.sites * sizeof(int));  /* index on node */
  parity = (int *)memalloc(node.sites * sizeof(int));  /* and parity */
  here   = (int *)memalloc(node.sites * sizeof(int));  /* index of original site */
  itmp   = (int *)memalloc(node.sites * sizeof(int));  /* temporary sorting array */
  free(tmpbuf);   /* release here to give "low memory", avoid fragmentation */

  c_offset = node.sites;  /* current offset in arrays */

  /* now, construct nn-gather to direction d */

#ifdef SCHROED_FUN
  /* special case for Schroedinger functional: fixed boundary
   * This requires we allocate one extra space for
   * CONSTANT schroedinger b.c.
   * Thus, add 1 item to latfield arrays to
   * take the boundary into account, and set the
   * neighb[d] to point to it.  This is for mpi or vanilla.
   * NOTE: IF NON-CONSTANT BOUNDARY NEED TO ADD c_offset FOR
   * EACH SITE, and set the corresponding neighb
   */
  int sf_special_boundary = c_offset;
  /* add space for boundary only when needed */
  if (node.xmin[NDIM-1] + node.nodesize[NDIM-1] == lat->size[NDIM-1])
    c_offset += 1;
  else sf_special_boundary = -(1<<30);

  fprintf0(outf," SPECIAL BOUNDARY LAYOUT for SF\n");
#endif

  for (d=0; d<NDIRS; d++) {

    comlist[d].index = neighb[d];  /* the neighbour array - not needed, but consistency*/
    comlist[d].label = d;          /* and label it */

    /* first pass over sites */
    num = 0;  /* number of sites off node */
    forallsites(i) {
      foralldir(j) x[j] = site[i].x[j];
      if (is_up_dir(d)) {
	x[d] = (x[d] + 1) % lat->size[d];   /* neighbour of site */
      } else {
	k = opp_dir(d);
	x[k] = (x[k] - 1 + lat->size[k]) % lat->size[k]; /* neighbour of site */
      }
#ifdef SCHROED_FUN
      if (d == NDIM-1 && x[NDIM-1] == 0) {
	/* This is up-direction, give special site */
	neighb[d][i] = sf_special_boundary;
      } else if (d == opp_dir(NDIM-1) && x[NDIM-1] == lat->size[NDIM-1]-1) {
	/* We never should need the down direction, thus, short circuit! */
	neighb[d][i] = 1<<30;
      } else     /* NOTE THIS UGLY else-if CONSTRUCT! */
#endif
      if (is_on_node(x)) {
	neighb[d][i] = node_index(x,&node);
      } else {
	/* Now site is off-node, this lead to fetching */
	nodes[num] = node_number(x);
	index[num] = node_index(x, allnodes + nodes[num] );
	parity[num] = site[i].parity;  /* parity of THIS */
	here[num] = i;
	num++;
      }
    }

    comlist[d].n_receive = 0;
    if (num > 0) {
      /* now, get the number of nodes to be gathered from */
      for (i=0; i<num; i++) {

	/* chase the list until the node found */
	for (j=0, r=&(comlist[d].from_node); j<comlist[d].n_receive &&
	       nodes[i] != (*r)->node; j++) r = &((*r)->next);
	if (j == comlist[d].n_receive) {
	  /* NEW NODE to receive from */
	  comlist[d].n_receive++;
	  (*r) = p = (receive_struct *)memalloc( sizeof(receive_struct) );
	  /* and fill in the node structure */
	  p->node = nodes[i];
	  p->n = 1; /* first site */
	  p->n_even = p->n_odd = 0;
	  if ( parity[i] == EVEN ) p->n_even = 1;  else p->n_odd = 1;
	  p->next = NULL;
	} else {
	  /* add to OLD NODE */
	  p = *r;
	  p->n ++;
	  if ( parity[i] == EVEN ) p->n_even ++;  else p->n_odd ++;
	}
      }

      /* Calculate the offsets for the gathers */
      for (j=0, p=comlist[d].from_node; j<comlist[d].n_receive; j++, p = p->next) {
	p->offset = c_offset;
	c_offset += p->n;  /* and increase the offset */
      }

      /* and NOW, finish the NEIGHBOR array */

      for (j=0, p=comlist[d].from_node; j<comlist[d].n_receive; j++, p = p->next) {
	/* Now, accumulate the locations to itmp-array, and sort the
	 * array according to the index of the sending node .
	 * First even neighbours
	 */
	for (par=EVEN; par<=ODD; par++) {
	  for (n=i=0; i<num; i++) if (nodes[i] == p->node && parity[i] == par) {
	    itmp[n++] = i;
	    /* bubble sort the tmp-array */
	    for (k=n-1; k > 0 && index[itmp[k]] < index[itmp[k-1]]; k--)
	      swap( itmp[k], itmp[k-1] );
	  }
	  off = p->offset;
	  if (par == ODD) off += p->n_even;
	  /* finally, root indices according to offset */
	  for (k=0; k<n; k++) neighb[d][here[itmp[k]]] = off + k;
	}
      }
    } /* num > 0 */

#ifndef LU
    //Check that we are not dividing odd lattices, breaks communication
    if (num > n_sublattices) if((nx%2==1)||(ny%2==1)||(nz%2==1)||(nt%2==1))
    halt("Odd number of sites cannot be devided into nodes.");
#endif

    /* receive done, now opposite send. This is just the gather
     * inverted
     */

    od = opp_dir(d);
    comlist[od].n_send = comlist[d].n_receive;

    if (num > 0) {
      p = comlist[d].from_node;
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


/***************************************************************************
 * A long routine to make general gather datastructs (comlists)
 * and neighbour/index arrays
 *
 * fmap is a function mapping
 *    x->xp; if mode FORWARDS
 *    xp->x; if mode INVERSE
 * arg is an int-array of passthrough arguments to fmap
 *
 * example: mapping (x,y,z) -> (x+1,y+1,z+1)   FORWARDS
 *                  (x,y,z) -> (x-1,y-1,z-1)   INVERSE
 *
 * Return value of make_gather is a label for the gather.
 *
 * After making the gather it can be used with
 *  start_gather( U, label, parity, ptr_arr )
 * defined in com_mpi.c
 */

int make_gather
( void (*fmap)(int x[NDIM], int *arg, int mode, int xp[NDIM]),
  int *arg,
  int mode        /* FORWARDS, INVERSE */
  )
{
  comlist_struct *cp;
  int gather_index,i,j,k,c_offset,x[NDIM];
  int *nodes, *index, *parity, *here, *idx, *itmp, *tbuf, num;
  int par,n,off;
  for (i=0; i<MAX_GATHERS && comlist[i].label >= 0; i++) ;
  gather_index = i;

  if (gather_index >= MAX_GATHERS) {
    fprintf0(outf," Too many gathers: max %d\n",MAX_GATHERS-NDIRS);
    terminate(1212);
  }

  cp = comlist + gather_index;
  cp->label = gather_index;      /* Mark it in use */
  
  /* First, allocate neighbour/index arrays */
  idx = cp->index = (int *)memalloc(node.sites * sizeof(int));
  #ifdef CUDA
    cp->d_index = (int *)parallel_alloc(MemType_DEV , sizeof(int) * node.sites);
  #endif

  /* work arrays, will be released later */

  tbuf   = (int *)memalloc(node.sites * sizeof(int)); /* space holder */
  nodes  = (int *)memalloc(node.sites * sizeof(int)); /* node number */
  index  = (int *)memalloc(node.sites * sizeof(int)); /* index on node */
  parity = (int *)memalloc(node.sites * sizeof(int)); /* and parity */
  here   = (int *)memalloc(node.sites * sizeof(int)); /* index of original site */
  itmp   = (int *)memalloc(node.sites * sizeof(int)); /* temporary index array */
  free(tbuf); /* and release space-try to avoid fragmentation */

  c_offset = node.sites;  /* current offset in arrays */

  /* now, construct gather */

  /* first pass over sites */
  num = 0;  /* number of sites off node */
  forallsites(i) {
    fmap( site[i].x, arg, mode, x );

    if (is_on_node(x)) idx[i] = node_index(x,&node);
    else {
      nodes[num] = node_number(x);
      index[num] = node_index(x, allnodes + nodes[num] );
      parity[num] = site[i].parity;  /* parity of THIS */
      here[num] = i;
      num++;
    }
  }

  cp->n_receive = 0;
  if (num > 0) {
    receive_struct **r,*p;
    /* now, get the number of nodes to be gathered from */
    for (i=0; i<num; i++) {

      /* chase the list until the node found */
      for (j=0, r=&(cp->from_node); j<cp->n_receive &&
	     nodes[i] != (*r)->node; j++) r = &((*r)->next);
      if (j == cp->n_receive) {
	/* NEW NODE to receive from */
	cp->n_receive++;
	(*r) = p = (receive_struct *)memalloc(sizeof(receive_struct));
	/* and fill in the node structure */
	p->node = nodes[i];
	p->n = 1; /* first site */
	p->n_even = p->n_odd = 0;
	if ( parity[i] == EVEN ) p->n_even = 1;  else p->n_odd = 1;
	p->next = NULL;
      } else {
	/* add to OLD NODE */
	p = *r;
	p->n ++;
	if ( parity[i] == EVEN ) p->n_even ++;  else p->n_odd ++;
      }
    }
    /* Calculate the offsets for the gathers */
    for (j=0, p=cp->from_node; j<cp->n_receive; j++, p = p->next) {
      p->offset = c_offset;
      c_offset += p->n;  /* and increase the offset */
    }

    /* and NOW, finish the index array */

    for (j=0, p=cp->from_node; j<cp->n_receive; j++, p = p->next) {
      /* Now, accumulate the locations to itmp-array, and sort the
       * array according to the index of the sending node .
       * First even neighbours
       */
      for (par=EVEN; par<=ODD; par++) {
	for (n=i=0; i<num; i++) if (nodes[i] == p->node && parity[i] == par) {
	  itmp[n++] = i;
	  /* bubble sort the tmp-array */
	  for (k=n-1; k > 0 && index[itmp[k]] < index[itmp[k-1]]; k--)
	    swap( itmp[k], itmp[k-1] );
	}
	off = p->offset;
	if (par == ODD) off += p->n_even;
	/* finally, root indices according to offset - GIVES OFFSETS ABOVE THE
	 * ARRAY LENGTH
	 */
	for (k=0; k<n; k++) idx[here[itmp[k]]] = off + k;
      }
    }
  } /* num > 0 */

  /* total size of receive buffer */
  cp->bufsiz = c_offset-node.sites;
  if (cp->bufsiz != num) fprintf(outf,"bufsiz error!\n");

  /* receive done, now the send.  Again pass over the sites, collecting
   * nodes.  Use inverse function
   */
  num = 0;  /* number of sites off node */
  forallsites(i) {

    fmap( site[i].x, arg, OPP_MODE(mode), x );

    if (!is_on_node(x)) {
      nodes[num]  = node_number(x);
      index[num]  = node_index( x, allnodes + nodes[num] );
      parity[num] = coordinate_parity(x);  /* parity of the coordinate to sent to! */
      here[num]   = i;
      num++;
    }
  }

  cp->n_send = 0;
  if (num > 0) {
    send_struct **s,*q;
    /* now, get the number of nodes to be gathered from */
    for (i=0; i<num; i++) {

      /* chase the list until the node found */
      for (j=0, s=&(cp->to_node); j<cp->n_send &&
	     nodes[i] != (*s)->node; j++) s = &((*s)->next);
      if (j == cp->n_send) {
	/* NEW NODE to send to */
	cp->n_send++;
	(*s) = q = (send_struct *)memalloc(sizeof(send_struct));
	/* and fill in the node structure */
	q->node = nodes[i];
	q->n = 1; /* first site */
	q->n_even = q->n_odd = 0;
	if ( parity[i] == EVEN ) q->n_even = 1;  else q->n_odd = 1;
	q->next = NULL;
      } else {
	/* add to OLD NODE */
	q = *s;
	q->n ++;
	if ( parity[i] == EVEN ) q->n_even ++;  else q->n_odd ++;
      }
    }

    /* allocate & initialize sitelist-buffers
     * - must be done outside the loop above
     */
    for (j=0, q=cp->to_node; j<cp->n_send; j++, q = q->next) {
      q->sitelist = (int *)memalloc(q->n * sizeof(int));
      /* now, initialize sitelist -- Now, we first want EVEN parity, since
       * we label the parity by the node to which we send to.
       * Irrelevant if parity is not needed
       */
      for (n=0,par=EVEN; par<=ODD; par++) {
        for (i=0; i<num; i++) if (nodes[i] == q->node && parity[i] == par) {
          (q->sitelist)[n++] = here[i];
        }
      }
      // Handle parallel-stuff
      {
        int err = parallel_initDevSitelist(q);
        if (err < 0) printf("parallel_initDevSitelist(q) failed! %s:%d\n", __FILE__, __LINE__);
      }
    }
  } /* num > 0 */

  // Indices
  #ifdef CUDA
  parallel_copy(cp->d_index, MemType_DEV, cp->index, MemType_HOST, sizeof(int) * node.sites);
  #endif
  /* send done
   */

  free(nodes);
  free(index);
  free(parity);
  free(here);
  free(itmp);


  /* put in a test too - can be taken away */
  // TODO: Does not work correctly in structure of arrays layout-case. This should be fixed
  //#if 0// !STRUCT_OF_ARRAYS_LAYOUT
  #ifndef CUDA
  #ifndef HREP // for some reason hrep+test_gather does not work
  test_gather( gather_index, fmap, arg, mode );
  #endif
  #endif
  //#endif
  /*
   * if (numnodes() > 1) {
   *   fprintf0(outf," Gather %d initialized: node0 ->(%d nodes, %d sites),",
   *            gather_index, cp->n_send, cp->bufsiz);
   * fprintf0(outf," <-(%d nodes, %d sites)\n",
   *          cp->n_receive, num);
   * }
   */

  return( gather_index );

}

/************************************************************************
 * delete gather permanently (need make-gather again)
 * NOT YET DONE
 */

/************************************************************************
 * Do gathers where offset is given by an offset array
 * This MUST be cleared with cleanup_gather( tag )!!!
 */

/* Helper function to make gathers */
/* max number of offset gathers */
/* Set SF boundary condition gathers here now too! */

#define N_G_MAX 120

void offset_fmap(int x[NDIM], int *offset, int mode, int xp[NDIM])
{
  int d;

  if (mode == FORWARDS) {
    foralldir(d) {
      xp[d] = ( x[d] + offset[d] + 6*lattice.size[d] ) % lattice.size[d];
    }
  } else {
    foralldir(d) {
      xp[d] = ( x[d] - offset[d] + 6*lattice.size[d] ) % lattice.size[d];
    }
  }
}

msg_tag * start_gather_offset_func(
       /* arguments */
       char * field,  /* pointer to some latfield */
       int size,      /* size in bytes of the field (eg sizeof(su3_vector))*/
       int offset[NDIM],  /* offsets to gather from */
       int parity,    /* parity of sites whose neighbors we gather.
			 one of EVEN, ODD or EVENODD (EVEN+ODD). */
       fieldptr_array **out_fieldptrs
                       // So this is a pointer to a variable, which receives the
                       // pointer to an array of pointers (some indirection)
                       // But it's actually - we just want to return an array of
                       // pointers, and since the output of this function is
                       // "already used" - we just write it into a location provided
                       // by the user of the function
       )
{
  static int ngathers = 0;
  static unsigned int hash_arr[N_G_MAX];
  static int dir_arr[N_G_MAX];
  int i,d;
  unsigned int hash;

  /* calculate a compressed "hash" */
  hash = 0;
  foralldir(d) {
    if (abs(offset[d]) > 120) {
      fprintf0(outf,"Too large offset for start_gather_offset, max 120, is %d\n",
	      offset[d]);
      finishrun();
    }
    hash = (hash<<8)|(offset[d]+127);
  }
  /* search for the hash */
  for (i=0; i<ngathers && hash_arr[i] != hash; i++) ;
   
  if (i == ngathers) {
    /* now did not find, make new */
    if (i >= N_G_MAX) {
      fprintf0(outf," Too many offset gathers: max %d\n",N_G_MAX);
      finishrun();
    }
    hash_arr[i] = hash;
    dir_arr[i]  = make_gather( offset_fmap, offset, FORWARDS );
    ngathers++;
  }
  /* and call general gather function */
  return( start_gather_func( field, size, dir_arr[i], parity, out_fieldptrs ) );
}


/************************************************************************
 * Do some test to validate the correctness of the gather
 */




void gather_test_error( char *abuse, int dir, tst_struct *a,
			tst_struct *n, int par )
{
  int l;

  fprintf(outf," *** %s, parity %d, from dir %d: ( ",abuse,par,dir);
  foralldir(l) fprintf(outf,"%d ",a->x[l]);
  fprintf(outf,") -> ( ");
  foralldir(l) fprintf(outf,"%d ",n->x[l]);
  fprintf(outf,"), parity %d -> %d\n",a->parity,n->parity);
}


void run_std_gather_test(void)
{
    test_std_gathers(&lattice);
    test_simple_histogram();
}


void test_std_gathers( lattice_struct *lat )
{
  int i,d,k,j,n,off,dir,n_err,par,checkparity;
  tst_struct *a;
  tst_struct *d_a;
  msg_tag *tag[NDIM];
#ifdef CUDA
  const memHeader* header;
#endif
  d_a = new_latfield( tst_struct );
  #ifdef CUDA
    header = getFieldHeader_h(d_a);
      a = (tst_struct*)malloc(sizeof(tst_struct)*node.latfield_size);
  #else
    a = d_a;
  #endif
  /* ignore parity if blocked lattice - usually OK */
  checkparity = 1;
  foralldir(d) if (current_blocking_level[d]) checkparity = 0;

  n_err = 0;
  for (k=0; k<2; k++) {
    for (par=EVEN; par<=EVENODD; par++) {

      forallsites(i) {
	foralldir(d) a[i].x[d] = site[i].x[d];
	a[i].parity = site[i].parity;
      }
#ifdef CUDA
      // Copy the data to actual field buffer
      parallel_copy(d_a, MemType_DEV, a, MemType_HOST, sizeof(tst_struct) * node.sites);
#endif

      foralldir(d) {
	if (k) dir = opp_dir(d); else dir = d;
	tag[d]  = start_get( d_a, dir, par );
      }

      foralldir(d) {
	if (k) dir = opp_dir(d); else dir = d;

	wait_get(tag[d]);
#ifdef CUDA
      // Copy the data from the actual field to host-side buffer - note: Copy all sites (on and off-node)
      parallel_copy(a, MemType_HOST, d_a, MemType_DEV, sizeof(tst_struct) * node.latfield_size);
#endif

	if (is_up_dir(dir)) off = 1; else off = lat->size[d] - 1;

	forparity(i,par) {
#ifdef SCHROED_FUN
	  if ((dir != NDIM-1 || a[i].x[d] != lat->size[d]-1) &&
	      (dir != opp_dir(NDIM-1) && a[i].x[d] != 0))
#endif
	    foralldir(j) {
	      n = nb(dir,i);
	      if (( j != d && a[n].x[j] != a[i].x[j]) ||
		  ( j == d && a[n].x[j] != ((a[i].x[j] + off) % lat->size[d]))
#ifndef IGNORE_PARITY
		  || (( a[i].parity != opp_parity(a[n].parity))&&checkparity )
#endif
		  ) {
		if (n_err < 10)
		  gather_test_error("Hello! NN-map error",dir,a+i,a+n,par);
		n_err ++;
	      }
	    }
	}
      }
    }
  }
  // TODO-CUDA: Reimplement rest of the test
   #ifdef CUDA
    free_latfield(d_a);
    return;
  #endif

  /* test scatter too - inverse.  Sensible only for EVEN or ODD */
  /* can be up or down */
  for (dir=0; dir<NDIRS; dir++) {
    int odir = opp_dir(dir);
    for (par=EVEN; par<=EVENODD; par++) {
      int opar = opp_parity(par);
      int error = 0;

      forallsites(i) foralldir(d) a[i].x[d] = 0;

      forparity(i,par) {
	/* put other parity (neighb) sites to data values */
#ifdef SCHROED_FUN
	if ((site[i].x[NDIM-1] != 0 || dir != opp_dir(NDIM-1)) &&
	    (site[i].x[NDIM-1] != lat->size[NDIM-1]-1 ||
	     dir != NDIM-1))
#endif
	  {
	    n = nb(dir,i);
	    foralldir(d) {
	      if (d == dir)
		a[n].x[d] = ((site[i].x[d] + 1)%lat->size[d]);
	      else if (d == odir)
		a[n].x[d] = ((site[i].x[d] - 1 + lat->size[d])%lat->size[d]);
	      else
		a[n].x[d] = site[i].x[d];
	    }
	    a[n].parity = opp_parity(site[i].parity);
	  }
      }

      wait_put( start_put( a, dir, par ) );
#ifdef LU
      forparity(i,opar) {
#ifdef SCHROED_FUN
	if ((site[i].x[NDIM-1] != 0 || dir != (NDIM-1)) &&
	    (site[i].x[NDIM-1] != lat->size[NDIM-1]-1 ||
	     dir != opp_dir(NDIM-1)))
#endif
	  {
	    error = 0;
#ifndef IGNORE_PARITY
	    if (checkparity && a[i].parity != site[i].parity) error = 1;
#endif
	    foralldir(d) if (a[i].x[d] != site[i].x[d]) error = 1;
	    if (error) {
              if (n_err < 10) gather_test_error("HALOO! Scatter error",
						dir,a+i,a+nb(odir,i),par);
	      n_err ++;
	    }
	  }
      }
#endif
    }
  }

  //if (n_err > 0) halt(" Lattice layout error (BUG in com_mpi.c or layout.c)");
  /* else fprintf0(outf," Gather/Scatter tests passed\n"); */

  free_latfield(a);
}

/****************************************************************
 * Test general gather
 */

void error_gather(int gat, int i, int d, int x[NDIM],
                  tst_struct *tp, int parity, const char* msg)
{
  if (msg)
    fprintf(outf," Gather (%s) %d maps ( ",msg, gat);
  else
    fprintf(outf," Gather (normal) %d maps ( ",gat);
  foralldir(d) fprintf(outf,"%d ",site[i].x[d]);
  fprintf(outf,") -> ( ");
  foralldir(d) fprintf(outf,"%d ",tp[i].x[d]);
  fprintf(outf,") should be ( ");
  foralldir(d) fprintf(outf,"%d ",x[d]);
  fprintf(outf,") parity %s\n",
          parity == EVEN ? "EVEN" : (parity == ODD ? "ODD" : "EVENODD") );

  fprintf(outf," Node %d, site %d\n", this_node, i);

  terminate(2121);
}

#if 1 //!STRUCT_OF_ARRAYS_LAYOUT


typedef struct copyTstToAOS_in_s
{
    fieldptr_array* tst_struct_ptrs;
    tst_struct* dst;
} copyTstToAOS_in;

PARALLEL_CALL_KERNEL_BEGIN(copyTstToAOS, copyTstToAOS_in, input, i)
{
#if STRUCT_OF_ARRAYS_LAYOUT
  tst_struct* in = (tst_struct*)input.tst_struct_ptrs[i].ptrToFirst;
#else
  tst_struct* in = (tst_struct*)input.tst_struct_ptrs[i];
#endif
    input.dst[i] = *in;
}
PARALLEL_CALL_KERNEL_END();

typedef struct tstToFmatrix_in_s
{
  tst_struct* src;
  fmatrix_array* dest;
} tstToFmatrix_in;

PARALLEL_CALL_KERNEL_BEGIN(tstToFmatrix, tstToFmatrix_in, input, i)
{
  // Store tst-struct in fmatrices (at least as much as we can)
  radix* psrc = (radix*)(&input.src[i]);
  fmatrix mat;
  radix* res = (radix*)&mat;
  int j;
  const int size1 = sizeof(tst_struct) / sizeof(radix);
  const int size2 = sizeof(fmatrix) / sizeof(radix);
  const int size = size1 < size2 ? size1 : size2;
  fm_zero(&mat);
  for (j = 0; j < size; j++)
  {
    res[j] = psrc[j];
  }
  set_fmatrix_element(input.dest, &mat, i);
}
PARALLEL_CALL_KERNEL_END();

typedef struct fmatToTst_in_s
{
  fmatrixptr_array * sources;
  tst_struct * dst;
} fmatToTst_in;

PARALLEL_CALL_KERNEL_BEGIN(fmatToTst, fmatToTst_in, input, i)
{
  fmatrix mat = get_fmatptr_element(input.sources, i);
  tst_struct res = input.dst[i];
  radix* psrc = (radix*)(&mat);
  radix* pdst = (radix*)(&res);
  int j;
  const int size1 = sizeof(tst_struct) / sizeof(radix);
  const int size2 = sizeof(fmatrix) / sizeof(radix);
  const int size = size1 < size2 ? size1 : size2;
  for (j = 0; j < size; j++)
  {
    pdst[j] = psrc[j];
  }
  // Store result!
  input.dst[i] = res;
}
PARALLEL_CALL_KERNEL_END();


void test_gather(int gat,
		 void (*fmap)(int x[NDIM], int *arg, int mode, int xp[NDIM]),
		 int *arg, int mode)
{
  tst_struct *a, *d_a;
  tst_struct *result;
  fieldptr_array* ptrs;
  fmatrixptr_array* mptrs;
  fmatrix_array * mat_a;
  int i,d,x[NDIM];
  msg_tag * tag;

  tst_struct* tmp = (tst_struct *)parallel_alloc(MemType_DEV, node.sites * sizeof(tst_struct) );

  a = (tst_struct *)memalloc( node.sites * sizeof(tst_struct) );
  result = (tst_struct *)memalloc( node.sites * sizeof(tst_struct) );
  d_a = new_latfield( tst_struct );
  mat_a = new_latfield( fmatrix_array );

  //printf("test_gather, gat = %d\n", gat);

  forallsites(i) {
    foralldir(d) a[i].x[d] = site[i].x[d];
    a[i].parity = site[i].parity;
  }

  // Copy the data to actual field buffer
  parallel_copy(d_a, MemType_DEV, a, MemType_HOST, sizeof(tst_struct) * node.sites);
  // And to fmatrices:
  {
    tstToFmatrix_in input;
    input.dest = mat_a;
    input.src = d_a;
    FOR_ALL_SITES_KERNEL_CALL(tstToFmatrix, input);
  }

  tag = start_gather( d_a, gat, EVEN, ptrs );
  wait_gather( tag );

  {
    copyTstToAOS_in input;
    input.tst_struct_ptrs = ptrs;
    input.dst = tmp;
    FOR_EVEN_SITES_KERNEL_CALL(copyTstToAOS, input);
    parallel_copy(result, MemType_HOST, tmp, MemType_DEV, node.sites * sizeof(tst_struct) );
  }
  forevensites(i) {
    fmap( site[i].x, arg, mode, x );
    foralldir(d) if ( x[d] != result[i].x[d] ) error_gather(gat,i,d,x, result,EVEN, NULL);
  }
  cleanup_gather( tag );


  // Now do the same using fmatrices as transport vehicle
  tag = start_gather( mat_a, gat, EVEN, mptrs);
  wait_gather( tag );
  {
    fmatToTst_in input;
    input.dst = tmp;
    input.sources = mptrs;
    FOR_EVEN_SITES_KERNEL_CALL(fmatToTst, input);
    parallel_copy(result, MemType_HOST, tmp, MemType_DEV, node.sites * sizeof(tst_struct) );
  }
  forevensites(i) {
    fmap( site[i].x, arg, mode, x );
    foralldir(d) if ( x[d] != result[i].x[d] ) error_gather(gat,i,d,x, result,EVEN, "fmatrix transport 1");
  }
  cleanup_gather( tag );




  tag = start_gather( d_a, gat, ODD, ptrs );
  wait_gather( tag );

  {
    copyTstToAOS_in input;
    input.tst_struct_ptrs = ptrs;
    input.dst = tmp;
    FOR_ODD_SITES_KERNEL_CALL(copyTstToAOS, input);
    parallel_copy(result, MemType_HOST, tmp, MemType_DEV, node.sites * sizeof(tst_struct) );
  }
  foroddsites(i) {
    fmap( site[i].x, arg, mode, x );
    foralldir(d) if ( x[d] != result[i].x[d] ) error_gather(gat,i,d,x,result,ODD, NULL);
  }
  cleanup_gather( tag );

  // Now do the same using fmatrices as transport vehicle
  tag = start_gather( mat_a, gat, ODD, mptrs);
  wait_gather( tag );
  {
    fmatToTst_in input;
    input.dst = tmp;
    input.sources = mptrs;
    FOR_ODD_SITES_KERNEL_CALL(fmatToTst, input);
    parallel_copy(result, MemType_HOST, tmp, MemType_DEV, node.sites * sizeof(tst_struct) );
  }
  forevensites(i) {
    fmap( site[i].x, arg, mode, x );
    foralldir(d) if ( x[d] != result[i].x[d] ) error_gather(gat,i,d,x, result,ODD, "fmatrix transport 2");
  }
  cleanup_gather( tag );







  tag = start_gather( d_a, gat, EVENODD, ptrs );
  wait_gather( tag );

  {
    copyTstToAOS_in input;
    input.tst_struct_ptrs = ptrs;
    input.dst = tmp;
    FOR_ALL_SITES_KERNEL_CALL(copyTstToAOS, input);
    parallel_copy(result, MemType_HOST, tmp, MemType_DEV, node.sites * sizeof(tst_struct) );
  }
  forallsites(i) {
    fmap( site[i].x, arg, mode, x );
    foralldir(d) if ( x[d] != result[i].x[d] ) error_gather(gat,i,d,x,result,EVENODD, NULL);
  }
  cleanup_gather( tag );


  // Now do the same using fmatrices as transport vehicle
  tag = start_gather( mat_a, gat, EVENODD, mptrs);
  wait_gather( tag );
  {
    fmatToTst_in input;
    input.dst = tmp;
    input.sources = mptrs;
    FOR_ALL_SITES_KERNEL_CALL(fmatToTst, input);
    parallel_copy(result, MemType_HOST, tmp, MemType_DEV, node.sites * sizeof(tst_struct) );
  }
  forevensites(i) {
    fmap( site[i].x, arg, mode, x );
    foralldir(d) if ( x[d] != result[i].x[d] ) error_gather(gat,i,d,x, result,EVENODD, "fmatrix transport 3");
  }
  cleanup_gather( tag );





  free(a);
  free(result);
  parallel_free(tmp, MemType_DEV);
  free_latfield(d_a);
  free_latfield(mat_a);
}
#endif

/****************************************************************/

char *copy_latfield_func( char *f )
{
  char *t;
  const memHeader* header = getFieldHeader_h(f);
  size_t typesize = size_of_fieldtype( header->fieldtype, MemType_DEV);
  t = latfield_alloc( header->fieldtype );
  //memcpy( t, f, (siz * node.latfield_size) );
  // TODO: Fix storage location issues (if any)
  parallel_copy_field(t, MemType_DEV, f, MemType_DEV, (typesize * node.latfield_size));
  return( t );
}

/****************************************************************/


void latfield_free(char * field)
{
  //parallel_free_field(field, MemType_DEV);
  pool_free_field(field, MemType_DEV);
}

void hostfield_free(char * field)
{
  //parallel_free_field(field, MemType_DEV);
  pool_free_field(field, MemType_HOST);
}



static char* fieldAllocImpl(int field_type, MemType allocType)
{
  char *t;
  char f[150];
  size_t size = size_of_fieldtype(field_type, allocType);
  size_t allocsize = node.latfield_size * size + GATHER_STATUS_SIZE;

#ifndef ALIGN
  //t = (char *)parallel_alloc_field(MemType_DEV, allocsize, node.latfield_size, field_type, node.latfield_size);
  t = (char *)pool_alloc_field(allocType, allocsize, node.latfield_size, field_type, node.latfield_size);
  //t = (char *)malloc();
  if (t == NULL) {
    sprintf(f,"Could not allocate a latfield of %d bytes", (int)size);
    halt(f);
  }

#else
  int e;
  /* Align all to 32-byte boundary */
  if ((e=posix_memalign( (void **)&t, (size_t)32,
			 node.latfield_size * size
			 + GATHER_STATUS_SIZE)) != 0) {
    fprintf(outf,"Could not allocate a latfield of %d bytes\n",size);
    sprintf(f,"posix_memalign error code %d\n",e);
    halt(f);
  }
  if (((long)t) & 31L) {
    sprintf(f,"MEMORY NOT ALIGNED; alignment %d\n,",(int)(((long)t) & 31L));
    halt(f);
  }
#endif

#ifdef USE_MPI
  gather_status_reset( t, size );
#endif
  return( t );
}

char* hostfield_alloc(int field_type)
{
  return fieldAllocImpl(field_type, MemType_HOST);
}


char* latfield_alloc(int field_type)
{
  return fieldAllocImpl(field_type, MemType_DEV);
}

#if !defined(MEMALLOC_ALIAS) || defined(ALIGN)
// MEMALLOC_ALIAS defined in layout.h, ALIGN in Makefile

char *memalloc(int size)
{
  char *t;

#ifndef ALIGN
  t = (char *)malloc( size );
  if (t == NULL) {
    char f[150];
    sprintf(f,"Memalloc: could not allocate %d bytes",size);
    halt(f);
  }
#else
  int e;
  /* Align all to 32-byte boundary */
  if ((e=posix_memalign( (void **)&t, (size_t)32, size )) != 0) {
    char f[150];
    fprintf(outf,"Could not allocate %d bytes\n",size);
    sprintf(f,"posix_memalign error code %d\n",e);
    halt(f);
  }

  if (((long)t) & 31L) {
    char f[150];
    sprintf(f,"MEMORY NOT ALIGNED; alignment %d\n,",(int)( ((long) t) & 31L));
    halt(f);
  }

#endif

  return( t );
}

#endif

/*********************************************************/

void *halt(char *s)
{
  fprintf0(outf,"*** %s\n",s);
  terminate(0);
  return((void *)NULL);
}


