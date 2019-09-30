/* This routine reorders the nodes in 3*2^n blocks, which
 * can fit within one multicore node.  This
 * potentially reduces n.n. communications.
 *
 * Suitable for some Intel cpu's with 12 (24) cores per die
 *
 * NOTE: THIS ASSUMES MPI GIVES THE CORES IN ORDER; I.E. 
 * THE NODES 0-11 ARE THE CORES OF 1st 12- CORE CPU,
 * 12-23 2nd ETC.  IF NOT, THIS PROCESS IS SUBOPTIMAL!
 *
 * This takes the blocklat divisible by 3 first.
 *
 * New: if system not divisible by 3, use just block2 -layout
 *
 * Kari Rummukainen 8/2014
 */

#include "../plumbing/lattice.h"

void map_node_layout()
{
  int nodes,i,d,c[NDIM],bdim[NDIM],bsize[NDIM],block_c[NDIM],na;
  int blocksize,a,b,j;
  int order[NDIM],di;

  nodes = numnodes();

  foralldir(d) order[d] = d;
  
  /* Find direction div by 3 */
  j = -1;
  for (d=0; j<0 && d<NDIM; d++) if (blocklat[d] % 3 == 0) j=d;
  if (j >= 0) {
    fprintf0(outf," map_node_layout: BLOCK 3*2^(dim-1)\n");
    
    // set the order of the directions
    order[0] = j;
    order[j] = 0;
  
    /* recognize here the dimensions where we do blocks */ 
    bsize[j] = 3;  // this was 3
    bdim[j] = blocklat[j] / 3;
    blocksize = 3;
    di = 1;

  } else {

    fprintf0(outf," No mod 3 dimensions found, attempt map_node_layout: BLOCK 2^dim\n");

    blocksize = 1;
    di = 0;
  }
  
  na=0;
  // di set above
  for ( ; di<NDIM; di++) {
    d = order[di];
    bsize[d] = (blocklat[d] % 2 == 0) ? 2 : 1;  // these 2 or 1
    bdim[d] = blocklat[d]/bsize[d];
    blocksize *= bsize[d];
    na += (bsize[d]-1);
  }

  if (this_node == 0) {
    if (j>=0) {
      fprintf(outf," map_node_layout: logical node blocks 3x2^%d, %d nodes/block\n",
              na, blocksize);
    } else {
      fprintf(outf," map_node_layout: logical node blocks 2^%d, %d nodes/block\n",
              na, blocksize);
    }
    fprintf(outf,"   In coordinate directions: ");
    for(d=0; d<NDIM-1; d++) fprintf(outf,"%d x ",bsize[d]);
    fprintf(outf,"%d\n",bsize[NDIM-1]);
    fprintf(outf,"   Taking directions in order: ");
    for(d=0; d<NDIM-1; d++) fprintf(outf,"%d ",order[d]);
    fprintf(outf,"%d\n",order[NDIM-1]);
 
  }

  /* i is the MPI index of the node */
  for (i=0; i<nodes; i++) {
    /* now what is the reshuffled coordinate? */
    /* block number and index */
    b = i/blocksize;
    j = i-b*blocksize;

    /* and get cart. coordinates */
    /* block_c contans block cart. coordinate */
    a = 1;
    foralldir(di) {
      d = order[di];
      block_c[d] = (b/a) % bdim[d];
      a = a*bdim[d];
    }
    /* and c[d] the inside block coordinate */
    a = 1;
    foralldir(di) {
      d = order[di];
      c[d] = (j/a) % bsize[d];
      a = a*bsize[d];
    }

    /* now the cartesian coord of this node
     * is block_c[d]*bsize[d] + c[d]
     * Thus, now the natural cart. index is */
    a = 1; j = 0;
    foralldir(d) {
      j += a*(block_c[d]*bsize[d] + c[d]);
      a = a*blocklat[d];
    }
    
    /* now, map the array to nodes 
     * when i is the natural "coordinate" index of
     * the node, map_arr[i] gives the true MPI index
     * of it.  Thus,
     */
    map_arr[j] = i;

    /* and print debug info here */

    // fprintf0(outf," map %d -> %d, ( ", i, j);
    // foralldir(d) fprintf0(outf,"%d ",block_c[d]*bsize[d] + c[d]);
    // fprintf0(outf,")\n");
  }
}
