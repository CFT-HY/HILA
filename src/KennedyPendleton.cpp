#ifndef KENNEDYPENDLETON_INL_
#define KENNEDYPENDLETON_INL_



#ifdef SU2


#ifndef _SU2_H
#   include "../include/su2.h"
#endif

#ifndef _RANDOM_H
#   include "../include/random.h"
#endif

#ifndef _GENERIC_H
//#   include "../include/generic.h"
#endif

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

/*#ifndef PARALLEL_DEFS_H_
#	include "../parallellib/parallel_defs.h"
#endif*/

/*************************************************************
 *                                                           *
 *     Kennedy-Pendleton su2 heatbath code                   *
 *     p contains the 'staple'                               *
 *     act contains the resulting action                     *
 *                                                           *
 *     act = u*p', weight exp[-u*p']                         *
 *                                                           *
 ************************************************************/


#define looplim 200

#if USE_HOSTSIDE_RND && !defined(__CUDA_ARCH__)
#define TMPRND(SITE_IDX)    dran()
#else
#define TMPRND(SITE_IDX) taus_rnd( SITE_IDX )
#endif

//this is might be too much to have in a single kernel
int KennedyPendleton( su2_matrix *u, su2_matrix *p,
            int siteindex)
{
  int nloops;
  su2_matrix a,ua;
  radix e,f,r1,r2,r3,rd,b,theta;

  /* first, normalize and move to u [-] because want exp[up] */

  b = 1.0/sqrt(su2_sqr((*p)));
  su2_scalar_mult((*p),-b,ua);

  /*     generate random su2 variable according to Pendleton
   *     and Kennedy, PL 156B(1985)393
   */

  nloops = 0;
  do {
    r1 = log(1.0-TMPRND( siteindex ));
    r2 = log(1.0-TMPRND( siteindex ));
    r3 = cos(pi2*TMPRND( siteindex ));
    r3 *= r3;
    e  = 1.0 + (r1+r2*r3)*b;
    r3 = TMPRND( siteindex );
    f  = e - 2.0*sqr(r3) + 1.0;
  } while (nloops++ < looplim && f <= 0.0);

  if (nloops >= looplim) {
    //fprintf(outf," ---> K-P loop limit %d exceeded!\n",looplim);
    //fprintf(outf," ---> staple magnitude %g\n",1/b);
    /* exit(0); */
    e = 1e-9;
  }

  /* generate random vector on s2, radius sqrt(1.0-a.d**2) */
  a.d   = e;
  rd    = 1.0 - sqr(e);
  if (rd <= 0.0) {
    // TODO: Is this ok way to handle warnings? We are wasting 31 of 32 bits
    // Could we perhaps do better?
//    warning=1;
    //fprintf(outf," ---> neg. K-P sqrt, val %g\n",rd);
    rd = 0.0;
  }
  r3    = 2.0*TMPRND( siteindex ) - 1.0;
  a.c   = sqrt(rd)*r3;
  theta = pi2*TMPRND( siteindex );
  rd    = sqrt(rd*(1.0 - sqr(r3)));
  a.b   = rd*cos(theta);
  a.a   = rd*sin(theta);

  /* a = -u*p' => u = -a*p = a*ua */

  su2_mult_nn(a,ua,(*u));

  return(nloops);
}

#undef TMPRND

#endif

#endif //KENNEDYPENDLETON_INL_
