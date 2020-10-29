
/************************** monte.cpp *******************************/
/* Kennedy-Pendleton quasi heat bath on SU(2) subgroups */
/* MIMD version 3 */
/* T. DeGrand March 1991 */
/* modified by K.R 97 & 2002 & 2004 */

/* This does the heatbath with
 *       beta ( 1 - 1/Ncol Tr U S^+ )  S is the sum of staples.
 *  ->  -beta/Ncol Tr u (US+) = -beta/Ncol u_ab (US+)_ba
 *  here u is a SU(2) embedded in SU(N); u = su2 * 1.
 */

#include "SUN.h"

#pragma hila loop_function
void mult_su2_mat_vec_elem_n(
  matrix<2,2,cmplx<double>> &u,
  cmplx<double> &x0,
  cmplx<double> &x1)
{
  /* Multiplies the complex column spinor (x0, x1) by the SU(2) matrix u */
  /* and puts the result in (x0,x1).  */
  /* Thus x <- u * x          */
  /* C. DeTar 3 Oct 1990 */

  cmplx<double> z0, z1, t0, t1;

  t0 = x0; t1 = x1;
  
  z0 = u.c[0][0]*t0;
  z1 = u.c[0][1]*t1;
  x0 = z0 + z1;
  z0 = u.c[1][0] * t0;
  z1 = u.c[1][1] * t1;
  x1 = z0 + z1;

} /* mult_su2_mat_vec_elem_n */


#pragma hila loop_function
void left_su2_hit_n(
  matrix<2,2,cmplx<double>> &u,
  int p,
  int q,
  matrix<N,N,cmplx<double>> &link)
{
  /* link <- u * link */
  /* The 0 row of the SU(2) matrix u matches row p of the SU(3) matrix */
  /* The 1 row of the SU(2) matrix u matches row q of the SU(3) matrix */
  /* C. DeTar 18 Oct 1990 */

  for (int m = 0; m < 3; m++)
    mult_su2_mat_vec_elem_n(u, link.c[p][m], link.c[q][m]);

} /* left_su2_hit_n */


#pragma hila loop_function
double monte(
  matrix<N,N,cmplx<double>> &U, 
  matrix<N,N,cmplx<double>> &staple,
  double beta)
{
  /* Do K-P quasi-heat bath by SU(2) subgroups */
  int ina, inb;
  double utry,uhit;
  double xr1,xr2,xr3,xr4;
  double a0,a1,a2,a3;
  double v0,v1,v2,v3, vsq;
  double h0,h1,h2,h3;
  double r,r2,rho,z;
  double al,d, xl,xd;
  int  k, test;
  double b3;
  matrix<N,N,cmplx<double>> action;
  matrix<2,2,cmplx<double>> h;
  double pi2 = PI*2.0;

  b3=beta/N;

  utry = uhit = 0.0;

  /* now for the qhb updating */

  action = U * staple;

  /*  loop SU(2) subgroups */
  for (ina=0; ina<N-1; ina++) for (inb=ina+1; inb<N; inb++) {

      /* decompose the action into SU(2) subgroups using
       * Pauli matrix expansion
       * The SU(2) hit matrix is represented as
       * a0 + i * Sum j (sigma j * aj)
       */
      v0 =  action.c[ina][ina].re + action.c[inb][inb].re;
      v3 =  action.c[ina][ina].im - action.c[inb][inb].im;
      v1 =  action.c[ina][inb].im + action.c[inb][ina].im;
      v2 =  action.c[ina][inb].re - action.c[inb][ina].re;

      vsq = v0*v0 + v1*v1 + v2*v2 + v3*v3;

      z = sqrt(vsq );
      /* Normalize   u */
      v0 = v0/z; v1 = v1/z; v2 = v2/z; v3 = v3/z;
      /* end norm check--trial SU(2) matrix is a0 + i a(j)sigma(j)*/
      /* test
	 if(this_node == 0)printf("v= %e %e %e %e\n",v0,v1,v2,v3);
	 if(this_node == 0)printf("z= %e\n",z);
      */
      /* now begin qhb */
      /* get four random numbers */

      xr1 = log(1.0 - hila_random());
      xr2 = log(1.0 - hila_random());
      xr3 = hila_random();
      xr4 = hila_random();

      xr3 = cos(pi2*xr3);

      /*
      if(this_node == 0)printf("rand= %e %e %e %e\n",xr1,xr2,xr3,xr4);
      */

      /*
	generate a0 component of suN matrix

	first consider generating an su(2) matrix h
	according to exp(bg/Ncol * re tr(h*s))
	rewrite re tr(h*s) as re tr(h*v)z where v is
	an su(2) matrix and z is a real normalization constant
	let v = z*v. (z is 2*xi in k-p notation)
	v is represented in the form v(0) + i*sig*v (sig are pauli)
	v(0) and vector v are real

	let a = h*v and now generate a
	rewrite beta/3 * re tr(h*v) * z as al*a0
	a0 has prob(a0) = n0 * sqrt(1 - a0**2) * exp(al * a0)
      */
      al = b3*z;
      /*if(this_node == 0)printf("al= %e\n",al);*/

      /*
	let a0 = 1 - del**2
	get d = del**2
	such that prob2(del) = n1 * del**2 * exp(-al*del**2)
      */

      d = -(xr2  + xr1*xr3*xr3)/al;

      /* monte carlo prob1(del) = n2 * sqrt(1 - 0.5*del**2)
       * then prob(a0) = n3 * prob1(a0)*prob2(a0)
       */

      /* now  beat each  site into submission */
      if ((1.00 - 0.5*d) <= xr4*xr4) {
	if(al > 2.0) { /* k-p algorithm */
	  test=0;
	  for(k=0; k<20 && !test;k++) {
	    /*  get four random numbers */
	    xr1 = log(1.0 - hila_random());
	    xr2 = log(1.0 - hila_random());
	    xr3 = hila_random();
	    xr4 = hila_random();

	    xr3 = cos(pi2*xr3);

	    d = -(xr2 + xr1*xr3*xr3)/al;
	    if ((1.00 - 0.5*d) > xr4*xr4) test = 1;
	  }
	  utry += k;
	  uhit++;
#ifndef __CUDA_ARCH__
	  // TODO: Make parallel-warning out of this
	  //output0 << "site took 20 kp hits\n";
#endif
	} else { /* now al <= 2.0 */

	  /* creutz algorithm */
	  xl=exp((double)(-2.0*al));
	  xd= 1.0 - xl;
	  test=0;
	  for(k=0;k<20 && test == 0  ;k++) {
	    /*        get two random numbers */
	    xr1=hila_random();
	    xr2=hila_random();

	    r = xl + xd*xr1;
	    a0 = 1.00 + log((double)r)/al;
	    if((1.0 -a0*a0) > xr2*xr2) test = 1;
	  }
	  d = 1.0 - a0;
	  utry += k;
	  uhit++;
#ifndef __CUDA_ARCH__
	  //output0 << "site  took 20 creutz hits\n";
#endif
	} /* endif al */

      } else {
	/* direct hit */
	utry += 1.0;
	uhit += 1.0;
      }

      /*  generate full su(2) matrix and update link matrix*/

      /* find a0  = 1 - d*/
      a0 = 1.0 - d;
      /* compute r */
      r2 = 1.0 - a0*a0;
      r2 = fabs(r2);
      r  = sqrt(r2);

      /* compute a3 */
      a3=(2.0*hila_random() - 1.0)*r;

      /* compute a1 and a2 */
      rho = r2 - a3*a3;
      rho = fabs(rho);
      rho = sqrt(rho);

      /*xr2 is a random number between 0 and 2*pi */
      xr2 = pi2*hila_random();
      a1 = rho*cos((double)xr2);
      a2 = rho*sin((double)xr2);

      /* now do the updating.  h = a*v^dagger, new u = h*u */
      h0 = a0*v0 + a1*v1 + a2*v2 + a3*v3;
      h1 = a1*v0 - a0*v1 + a2*v3 - a3*v2;
      h2 = a2*v0 - a0*v2 + a3*v1 - a1*v3;
      h3 = a3*v0 - a0*v3 + a1*v2 - a2*v1;

      /* Elements of SU(2) matrix */
      h.c[0][0].re =  h0; h.c[0][0].im =  h3;
      h.c[0][1].re =  h2; h.c[0][1].im =  h1;
      h.c[1][0].re = -h2; h.c[1][0].im =  h1;
      h.c[1][1].re =  h0; h.c[1][1].im = -h3;

      /* update the link and 'action' */

      left_su2_hit_n(h,ina,inb,U);
      left_su2_hit_n(h,ina,inb,action);

    } /*  hits */

  /* check_unitarity( U );
   */
  return( uhit/utry );

} /* site */


