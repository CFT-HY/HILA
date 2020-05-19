




template<typename vector>
void dirac_staggered_diag_gNJL_force(
  const field<vector> &psi,
  const field<vector> &chi,
  field<double> &sforce,
  field<double> &pforce,
  int sign, parity par )
{
  sforce[ALL] = 0;
  pforce[ALL] = 0;
  if(par!=ODD)
  onsites(EVEN){
    cmplx<double> dot = chi[X].dot(psi[X]);
    sforce[X] = sforce[X] + dot.re;
    pforce[X] = pforce[X] + sign * dot.im;
  }
  if(par!=EVEN)
  onsites(ODD){
    cmplx<double> dot = chi[X].dot(psi[X]);
    sforce[X] = sforce[X] + dot.re;
    pforce[X] = pforce[X] - sign * dot.im;
  }
}



// Apply the diagonally
template<typename vector>
void dirac_staggered_diag_gNJL(
  const double mass,
  const field<vector> &v_in,
  field<vector> &v_out,
  field<double> &sigma,
  field<double> &pi,
  parity par, int sign)
{
  if(par!=ODD)
    v_out[EVEN] = v_out[X] + cmplx(mass+sigma[X], sign*pi[X])*v_in[X];
  if(par!=EVEN)
    v_out[ODD]  = v_out[X] + cmplx(mass+sigma[X],-sign*pi[X])*v_in[X];
}


// Apply the inverse of the diagonal
template<typename vector>
void dirac_staggered_diag_gNJL_inverse(
  const double mass,
  field<vector> &v_out,
  const field<double> &sigma,
  const field<double> &pi,
  parity par, int sign)
{
  if(par!=ODD)
    onsites(EVEN){
      element<double> sqr = (sigma[X]+mass)*(sigma[X]+mass) + pi[X]*pi[X];
      v_out[X] = cmplx(sigma[X]+mass,-sign*pi[X])/sqr * v_out[X];
  }
  if(par!=EVEN)
    onsites(ODD){
      element<double> sqr = (sigma[X]+mass)*(sigma[X]+mass) + pi[X]*pi[X];
      v_out[X] = cmplx(sigma[X]+mass, sign*pi[X])/sqr * v_out[X];
  }
}




template<typename vector, typename matrix>
class dirac_staggered_gNJL {
  private:
    double mass;
    field<double> staggered_eta[NDIM];
    field<double> &sigma, &pi;

    // Note array of fields, changes with the field
    field<matrix> (&gauge)[NDIM];
  public:

    using vector_type = vector;

    // Constructor: initialize mass, gauge and eta
    dirac_staggered_gNJL(dirac_staggered_gNJL &d) : gauge(d.gauge), mass(d.mass), sigma(d.sigma), pi(d.pi) {
      init_staggered_eta(staggered_eta);
    }
    dirac_staggered_gNJL(double m, field<matrix> (&g)[NDIM], field<double> &s, field<double> &p) : gauge(g), mass(m), sigma(s), pi(p) {
      init_staggered_eta(staggered_eta);
    }

    // Applies the operator to in
    void apply( const field<vector> & in, field<vector> & out){
      out[ALL] = 0;
      dirac_staggered_diag_gNJL(mass, in, out, sigma, pi, ALL, 1);
      dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, 1);
    }

    // Applies the conjugate of the operator
    void dagger( const field<vector> & in, field<vector> & out){
      out[ALL] = 0;
      dirac_staggered_diag_gNJL(mass, in, out, sigma, pi, ALL, -1);
      dirac_staggered_hop(gauge, in, out, staggered_eta, ALL, -1);
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<vector> & chi, const field<vector> & psi, field<matrix> (&gforce)[NDIM], field<double> &sforce, field<double> &pforce, int sign){
      dirac_staggered_calc_force(gauge, chi, psi, gforce, staggered_eta, sign, ALL);
      dirac_staggered_diag_gNJL_force(chi, psi, sforce, pforce, sign, ALL);
    }
};



template<typename vector, typename matrix>
class dirac_staggered_gNJL_evenodd {
  private:
    double mass;
    field<double> staggered_eta[NDIM];

    field<double> &sigma, &pi;
    field<matrix> (&gauge)[NDIM];
  public:

    using vector_type = vector;
    using matrix_type = matrix;

    dirac_staggered_gNJL_evenodd(dirac_staggered_gNJL_evenodd &d) : gauge(d.gauge), mass(d.mass), sigma(d.sigma), pi(d.pi) {
      init_staggered_eta(staggered_eta);
    }
    dirac_staggered_gNJL_evenodd(double m, field<matrix> (&g)[NDIM], field<double> &s, field<double> &p) : gauge(g), mass(m), sigma(s), pi(p) {
      init_staggered_eta(staggered_eta);
    }

    // Applies the operator to in
    inline void apply( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag_gNJL(mass, in, out, sigma, pi, EVEN, 1);

      dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, 1);
      dirac_staggered_diag_gNJL_inverse(mass, out, sigma, pi, ODD, 1);
      dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, 1);
      out[ODD] = 0;
    }

    // Applies the conjugate of the operator
    inline void dagger( const field<vector_type> & in, field<vector_type> & out){
      out[ALL] = 0;
      dirac_staggered_diag_gNJL(mass, in, out, sigma, pi, EVEN, -1);

      dirac_staggered_hop(gauge, in, out, staggered_eta, ODD, -1);
      dirac_staggered_diag_gNJL_inverse(mass, out, sigma, pi, ODD, -1);
      dirac_staggered_hop(gauge, out, out, staggered_eta, EVEN, -1);
      out[ODD] = 0;
    }

    // Applies the derivative of the Dirac operator with respect
    // to the gauge field
    void force( const field<vector> & psi, const field<vector> & chi, field<matrix> (&gforce)[NDIM], field<double> &sforce, field<double> &pforce, int sign){
      field<matrix_type> gforce2[NDIM];
      field<double> sforce2, pforce2;
      field<vector_type> tmp, tmp2;
      tmp.copy_boundary_condition(chi);
      tmp2.copy_boundary_condition(chi);
      dirac_staggered_diag_gNJL_force(psi, chi, sforce, pforce, sign, EVEN);

      tmp[ALL] = 0;
      dirac_staggered_hop(gauge, chi, tmp, staggered_eta, ODD, -sign);
      dirac_staggered_diag_gNJL_inverse(mass, tmp, sigma, pi, ODD, sign);
      dirac_staggered_calc_force(gauge, tmp, psi, gforce, staggered_eta, sign, EVEN);

      tmp2[ALL] = 0;
      dirac_staggered_hop(gauge, psi, tmp2, staggered_eta, ODD, sign);
      dirac_staggered_diag_gNJL_inverse(mass, tmp2, sigma, pi, ODD, -sign);
      dirac_staggered_calc_force(gauge, chi, tmp2, gforce2, staggered_eta, sign, ODD);

      dirac_staggered_diag_gNJL_force(tmp2, tmp, sforce2, pforce2, sign, ODD);

      foralldir(dir){
        gforce[dir][ALL] = gforce[dir][X] + gforce2[dir][X];
      }
      // Note negative sign since the second part is the derivative
      // of an inverse
      sforce[ALL] = sforce[X] - sforce2[X];
      pforce[ALL] = pforce[X] - pforce2[X];

    }
};

